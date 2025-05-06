##'
#' @title Pathway Enrichment Analysis
#' @description Performs over‐representation testing via DOSE/clusterProfiler against KEGG or custom gene sets.
#' @param gene Character vector of gene or protein IDs to test for enrichment.
#' @param organism KEGG organism code (e.g. `"ppu"`).
#' @param keyType ID type (default `"kegg"`).
#' @param pvalueCutoff Numeric; p‐value cutoff (default `0.05`).
#' @param pAdjustMethod Character; adjustment method (default `"BH"`).
#' @param universe Character vector of background IDs.
#' @param minGSSize Integer; minimum gene set size.
#' @param maxGSSize Integer; maximum gene set size.
#' @param qvalueCutoff Numeric; q‐value cutoff (default `0.2`).
#' @param use_internal_kegg Logical; use internal KEGG DB if `TRUE`.
#' @param custom_gene_sets Logical; use custom gene sets if `TRUE`.
#' @param custom_pathways Data frame with columns `"Pathway"` and `"Accession"`, or a GMT file path.
#' @return An `enrichResult` object (or `NULL` if no pathways pass the thresholds).
#' @export
#' @importFrom clusterProfiler read.gmt
#' @importFrom Biobase reverseSplit
#' @importFrom stats p.adjust phyper
#' @importFrom qvalue qvalue
#' @importFrom tibble deframe
#' @importFrom dplyr group_by summarise mutate rename
#' @importFrom magrittr %>%
pathway_enrich <- function (gene, organism = "ppu", keyType = "kegg",
                            pvalueCutoff = 0.05, pAdjustMethod = "BH", universe,
                            minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_kegg = FALSE, custom_gene_sets = FALSE, custom_pathways = NULL)
{
  # — allow custom_pathways to be a GMT filepath —
  if (custom_gene_sets &&
      is.character(custom_pathways) &&
      length(custom_pathways) == 1 &&
      grepl("\\.gmt$", custom_pathways) &&
      file.exists(custom_pathways))
  {
    # read.gmt returns a data.frame with columns `term` and `gene`
    gmt_df <- clusterProfiler::read.gmt(custom_pathways)
    # collapse genes into comma-space lists per pathway
    custom_pathways <- gmt_df %>%
      dplyr::group_by(term) %>%
      dplyr::summarise(Accession = paste(gene, collapse = ", "), .groups = "drop") %>%
      dplyr::rename(Pathway = term)
  }

  if (custom_gene_sets) {
    custom_pathways$Accession <-
      custom_pathways$Accession %>% str_replace_all(., " // ", ", ") %>%
      str_replace_all(., ", ", ",") %>% str_replace_all(., ",", ", ")
    custom_vec <-
      custom_pathways[, match(c("Pathway", "Accession"), colnames(custom_pathways))] %>% tibble::deframe()
    NAME2EXTID <- strsplit(as.character(custom_vec), ", ")
    names(NAME2EXTID) <- names(custom_vec)
    EXTID2NAME <- Biobase::reverseSplit(NAME2EXTID)
    DATA <- new.env()
    DATA$NAME2EXTID  <- NAME2EXTID
    DATA$EXTID2NAME  <- EXTID2NAME
    DATA$PATHID2EXTID <- NAME2EXTID
    res <- enricher_custom(gene, pvalueCutoff = pvalueCutoff,
                           pAdjustMethod = pAdjustMethod, universe = universe, minGSSize = minGSSize,
                           maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff, USER_DATA = DATA)
    if (is.null(res))
      return(res)
  } else {
    species <- clusterProfiler:::organismMapper(organism)
    if (use_internal_kegg) {
      DATA <- clusterProfiler:::get_data_from_KEGG_db(species)
    } else {
      DATA <- clusterProfiler:::prepare_KEGG(species, "KEGG", keyType)
    }
    res <- DOSE:::enricher_internal(gene, pvalueCutoff = pvalueCutoff,
                                    pAdjustMethod = pAdjustMethod, universe = universe, minGSSize = minGSSize,
                                    maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff, USER_DATA = DATA)
    if (is.null(res))
      return(res)
    res@ontology <- "KEGG"
    res@organism <- species
    res@keytype <- keyType
  }
  res <- dplyr::mutate(res, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
  return(res)
}





#' Enricher Custom
#'
#' This internal function performs pathway enrichment analysis using a custom set of pathways and genes/proteins.
#'
#' @param gene A character vector of gene or protein IDs found to be enriched.
#' @param pvalueCutoff The cutoff p-value for enrichment.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param universe background genes/proteins
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of each geneSet for analyzing.
#' @param qvalueCutoff The cutoff q-value (adjusted p value) for enrichment.
#' @param USER_DATA ontology information.
#'
#' @return An enriched results objects of class \code{enrichResult}
#'
#' @noRd
#' @keywords internal
enricher_custom <- function (gene, pvalueCutoff, pAdjustMethod = "BH", universe = NULL,
                             minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, USER_DATA)
{
  gene <- as.character(unique(gene))

  EXTID2NAME <- get("EXTID2NAME", envir = USER_DATA)
  qExtID2Name <- EXTID2NAME[gene]
  len <- sapply(qExtID2Name, length)
  notZero.idx <- len != 0
  qExtID2TermID <- qExtID2Name[notZero.idx]
  qTermID <- unlist(qExtID2TermID)
  if (is.null(qTermID)) {
    message("--> No gene can be mapped. Either no gene was found as part of any annotated custom pathway or gene IDs have the wrong format.")
    p2e <- get("NAME2EXTID", envir = USER_DATA)
    sg <- unlist(p2e[1:10])
    sg <- sample(sg, min(length(sg), 6))
    message("--> Expected input gene IDs (examples): ", paste0(sg,
                                                               collapse = ","))
    return(NULL)

  }
  qExtID2TermID.df <- data.frame(extID = rep(names(qExtID2TermID),
                                             times = lapply(qExtID2TermID, length)), termID = qTermID)
  qExtID2TermID.df <- unique(qExtID2TermID.df)
  extID <- get("NAME2EXTID", envir = USER_DATA) %>% unlist() %>% unique()

  qTermID2ExtID <- with(qExtID2TermID.df, split(as.character(extID),
                                                as.character(termID)))
  if (missing(universe)){
    universe <- NULL
  }
  if (!is.null(universe)) {
    if (is.character(universe)) {
      extID <- intersect(extID, universe)
    }
    else {
      message("`universe` is not in character and will be ignored...")
    }
  }
  qTermID2ExtID <- lapply(qTermID2ExtID, intersect, extID)
  qTermID <- unique(names(qTermID2ExtID))

  termID2ExtID <- get("PATHID2EXTID", envir = USER_DATA)
  termID2ExtID <- termID2ExtID[qTermID]
  termID2ExtID <- lapply(termID2ExtID, intersect, extID)
  geneSets <- lapply(termID2ExtID, intersect, extID)

  if (is.na(minGSSize) || is.null(minGSSize))
    minGSSize <- 1
  if (is.na(maxGSSize) || is.null(maxGSSize))
    maxGSSize <- Inf
  geneSet_size <- sapply(geneSets, length)
  idx <- minGSSize <= geneSet_size & geneSet_size <= maxGSSize


  if (sum(idx) == 0) {
    msg <- paste("No gene sets have size between",
                 minGSSize, "and", maxGSSize, "...")
    message(msg)
    message("--> return NULL...")
    return(NULL)
  }
  termID2ExtID <- termID2ExtID[idx]
  qTermID2ExtID <- qTermID2ExtID[idx]
  qTermID <- unique(names(qTermID2ExtID))
  k <- sapply(qTermID2ExtID, length)
  k <- k[qTermID]
  M <- sapply(termID2ExtID, length)
  M <- M[qTermID]
  N <- rep(length(extID), length(M))
  n <- rep(length(qExtID2TermID), length(M))
  args.df <- data.frame(numWdrawn = k - 1, numW = M, numB = N -
                          M, numDrawn = n)
  pvalues <- apply(args.df, 1, function(n) stats::phyper(n[1], n[2],
                                                         n[3], n[4], lower.tail = FALSE))
  GeneRatio <- apply(data.frame(a = k, b = n), 1, function(x) paste(x[1],
                                                                    "/", x[2], sep = "", collapse = ""))
  BgRatio <- apply(data.frame(a = M, b = N), 1, function(x) paste(x[1],
                                                                  "/", x[2], sep = "", collapse = ""))
  Over <- data.frame(ID = as.character(qTermID), GeneRatio = GeneRatio,
                     BgRatio = BgRatio, pvalue = pvalues, stringsAsFactors = FALSE)
  p.adj <- stats::p.adjust(Over$pvalue, method = pAdjustMethod)
  qobj <- tryCatch(qvalue::qvalue(p = Over$pvalue, lambda = 0.05, pi0.method = "bootstrap"),
                   error = function(e) NULL)
  if (class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  } else {
    qvalues <- NA
  }
  geneID <- sapply(qTermID2ExtID, function(i) paste(i, collapse = "/"))
  geneID <- geneID[qTermID]
  Over <- data.frame(Over, p.adjust = p.adj, qvalue = qvalues,
                     geneID = geneID, Count = k, stringsAsFactors = FALSE)
  Description <- qTermID
  if (length(qTermID) != length(Description)) {
    idx <- qTermID %in% names(Description)
    Over <- Over[idx, ]
  }
  Over$Description <- Description
  nc <- ncol(Over)
  Over <- Over[, c(1, nc, 2:(nc - 1))]
  Over <- Over[order(pvalues), ]
  Over$ID <- as.character(Over$ID)
  row.names(Over) <- as.character(Over$ID)
  x <- new("enrichResult", result = Over, pvalueCutoff = pvalueCutoff,
           pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,
           gene = as.character(gene), universe = extID, geneSets = geneSets,
           organism = "UNKNOWN", keytype = "UNKNOWN",
           ontology = "UNKNOWN", readable = FALSE)
  return(x)
}
