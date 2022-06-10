####____rna.read_data____####
rna.read_data <- function (data = "dat_prot.csv", # File or dataframe containing proteomics data
                            expdesign = NULL, # Experimental design as file path or data frame, if made previously
                            csvsep = ";", # optional: delimiter if reading CSV file
                            filter = c("Reverse", "Potential contaminant"),
                            name = 'Gene Symbol', # Header of column containing primary protein IDs
                            id = 'Ensembl Gene ID', # Header of column containing alternative protein IDs
                            pfx = "abundances.", # Prefix in headers of columns containing protein abundances
                            rsd_thresh = NULL,  # RSD filter in %!
                            filt_type = "condition",
                            filt_thr = 3, # keep genes that have a maximum of 'filt_thr' missing values in at least one condition.
                            filt_min = NULL # Sets the threshold for the minimum fraction of valid values allowed for any protein if type = "fraction".
) {
  assertthat::assert_that(is.character(name),
                          length(name) == 1,
                          is.character(id),
                          length(id) == 1)
  if (is.character(data)) {
    # Read table file
    if (file.exists(data)) {
      if (str_replace_all(data, ".{1,}\\.", "") == "csv") {
        prot <-
          utils::read.csv(
            data,
            sep = csvsep,
            header = T,
            stringsAsFactors = F,
            fill = T,
            na.strings = "",
            quote = "",
            comment.char = "",
            check.names = F
          )
      } else if (str_replace_all(data, ".{1,}\\.", "") == "xls" |
                 str_replace(data, ".{1,}\\.", "") == "xlsx") {
        prot <- readxl::read_excel(data)
      } else if (str_replace_all(data, ".{1,}\\.", "") == "tsv") {
        prot <-
          utils::read.csv(
            data,
            sep = "\t",
            header = T,
            stringsAsFactors = F,
            fill = T,
            na.strings = "",
            quote = "",
            comment.char = "",
            check.names = F
          )
      } else if (str_replace_all(data, ".{1,}\\.", "") == "txt") {
        prot <-
          utils::read.table(
            data,
            sep = "\t",
            header = T,
            stringsAsFactors = F,
            fill = T,
            na.strings = "",
            quote = "",
            comment.char = "",
            check.names = F
          )
      } else {
        stop(
          "No compatible file format provided.
             Supported formats are: \\.txt (tab delimited), \\.csv (delimiters can be specified with the argument \"csvsep = \", \\.tsv, \\.xls, and \\.xlsx"
        )
      }
    } else {
      stop(paste0("File \"", data, "\" does not exist."), call. = F)
    }
  } else {
    prot <- data
  }

  # Test for occurence of prefix for abundance columns.
  if (!(any(grepl(pfx, colnames(prot))))) {
    stop(paste0("The prefix '", pfx, "' does not exist in any column of '",
                data, "'. Please provide a valid prefix to identify columns with protein abundances."), call. = F)
  }

  # Filter out the positive genes (indicated by '+')
  # in the pre-defined "filter" columns
  cols_filt <- grep(paste("^", filter, "$", sep = "", collapse = "|"),  # The columns to filter on
                    colnames(prot))
  if (!is.null(cols_filt)) {
    if (length(cols_filt) == 1) {
      rows <- which(prot[, cols_filt] == "+")
      if(length(rows) > 0) prot <- prot[-rows,]
    } else {
      rows <- which(apply(prot[, cols_filt] == "+", 1, any))
      if(length(rows) > 0) prot <- prot[-rows,]
    }
  }

  # Make unique names using the annotations in the as name and id defined columns as primary and
  # secondary names, respectively.
  if (any(colnames(prot) %in% name) && any(colnames(prot) %in% id)) {
    # Remove rows that have 'NA' in both specified name and id columns
    row.rm <- which(is.na(prot[name]) & is.na(prot[id]), arr.ind = T)
    if (length(row.rm) > 1) {
      prot.rm <- prot[-c(row.rm),]
      message(
        paste0(
          "Removing genes with no defined name in column \"",
          name,
          "\" or \"",
          id,
          "\". (Removed: ",
          (nrow(prot) - nrow(prot.rm)),
          " genes)"
        )
      )
    }


    if (length(row.rm) > 1) {
      rna_unique <- prot.make_unique(proteins = prot.rm, names = name, ids = id, delim = ";")
    } else {
      rna_unique <- prot.make_unique(proteins = prot, names = name, ids = id, delim = ";")
    }

  } else if (!any(colnames(prot) %in% name) && !any(colnames(prot) %in% id)) {
    stop("\"", name, "\" and \"", id, "\" are not columns in \"", data, "\".", "Please provide valid column names with protein identifiers as \"name= \" and \"id= \".",
         call. = F)
  } else if (!any(colnames(prot) %in% name)) {
    stop("\"", name, "\" is not a column in \"", data, "\".", "Please provide a valid column name with protein names as \"name= \".",
         call. = F)
  } else if (!any(colnames(prot) %in% id)) {
    stop("\"", id, "\" is not a column in \"", data, "\".", "Please provide a valid column name with protein identifiers as \"id= \".",
         call. = F)
  }


  if (is.null(expdesign)) {
    # Create experimental design based on column names if neither file nor data frame is provided
    label <- rna_unique %>%
      select(., contains(pfx)) %>%
      colnames() %>%
      gsub(pfx, "", .)

    condition <- rna_unique %>%
      select(., contains(pfx)) %>%
      colnames() %>%
      gsub(pfx, "", .) %>%
      gsub(".[[:digit:]]+$", "", .)  # Remove prefix and replicate number from sample name

    replicate <- rna_unique %>%
      select(., contains(pfx)) %>%
      colnames() %>%
      str_extract(., "[:digit:]{1,}$")  # Remove string before replicate number

    experimental_design <- data.frame(label, condition, replicate, check.names = FALSE)
    message(paste0(
      "Writing experimental design file to: " ,
      getwd(),
      "/experimental_design.txt"
    ))
    utils::write.table(
      experimental_design,
      file = "experimental_design.txt",
      sep = "\t",
      row.names = F
      #col.names = NA
    )

  } else if (is.list(expdesign)) {
    experimental_design <- expdesign
  } else if (is.character(expdesign) && file.exists(expdesign)) {
    if (str_replace_all(expdesign, ".{1,}\\.", "") == "csv") {
      experimental_design <- utils::read.csv( expdesign, sep = csvsep, header = T,
                                              stringsAsFactors = F, fill = T, na.strings = "",
                                              quote = "", comment.char = "", check.names=FALSE )
    } else if (str_replace_all(expdesign, ".{1,}\\.", "") == "xls" |
               str_replace_all(expdesign, ".{1,}\\.", "") == "xlsx") {
      experimental_design <- readxl::read_excel(expdesign)
    } else if (str_replace_all(expdesign, ".{1,}\\.", "") == "tsv") {
      experimental_design <- utils::read.csv(expdesign, sep = "\t", header = T,
                                             stringsAsFactors = F, fill = T, na.strings = "",
                                             quote = "", comment.char = "", check.names=FALSE )
    } else if (str_replace_all(expdesign, ".{1,}\\.", "") == "txt") {
      experimental_design <- utils::read.table(expdesign, sep = "\t", header = T, stringsAsFactors = F,
                                               fill = T, na.strings = "", comment.char = "", check.names=FALSE )
    } else {
      stop(
        "No compatible file format for experimental design provided.
             Supported formats are: \\.txt (tab delimited), \\.csv (delimiters can be specified with the argument \"csvsep = \", \\.tsv, \\.xls, and \\.xlsx"
      )
    }
  } else {
    stop(paste0("File \"", expdesign, "\" does not exist."), call. = F)
  }

  # Generate a SummarizedExperiment object using an experimental design
  abundance_columns <- grep(pfx, colnames(rna_unique))  # get abundance column numbers
  message("Generating SummarizedExperiment.")
  rna_se <- rna.make_se(rna_unique, abundance_columns, experimental_design)

  # Apply RSD threshold
  if(!is.null(rsd_thresh)){
    rsd_thresh <- rsd_thresh/100

    int.mat <- assay(rna_se)
    # replace NA with 0
    int.mat[is.na(int.mat)] <- 0

    # Calculate the standard deviation for each compound in each group. Create a list of n vectors for n sample groups.
    ls.sd <- lapply(1:length(unique(se$condition)), function (x)
      apply(int.mat[,unique(se$condition)[x] == gsub("_[[:digit:]]+$", "", colnames(int.mat))], 1, stats::sd, na.rm = F))
    # Calculate the average for each compound in each group. Create a list of n vectors for n sample groups.
    ls.mns <-
      lapply(1:length(unique(se$condition)), function (x)
        apply(int.mat[,unique(se$condition)[x] == gsub("_[[:digit:]]+$", "", colnames(int.mat))], 1, mean, na.rm = F))
    # Calculate the RSD for each compound in each group. Create a list of n vectors for n sample groups.
    ls.rsd <- mapply("/", ls.sd, ls.mns, SIMPLIFY = FALSE)
    # Combine group vectors in RSD list into matrix
    mat.rsd <- do.call(rbind, ls.rsd)
    # Create vector with the maximum RSD value among all groups for each compound.
    filter.val <- apply(mat.rsd, 2, max, na.rm = T)
    # Check if each compound has an RSD value below the threshold.
    keep <- filter.val <= rsd_thresh

    # Apply RSD filter to data matrix
    filtered_rsd <- rna_se[keep, ]
  } else {
    filtered_rsd <- rna_se
  }
  if( (nrow(SummarizedExperiment::assay(rna_se)) - nrow(SummarizedExperiment::assay(filtered_rsd))) != 0 ){
    number_removed_rsd <- nrow(SummarizedExperiment::assay(rna_se)) - nrow(SummarizedExperiment::assay(filtered_rsd))
    cat(paste0(number_removed_rsd, " out of ",
               nrow(SummarizedExperiment::assay(se)), " genes were removed from the dataset due to too high RSD.\n\n"))
    filtered_rsd@metadata$n.filtered.rsd <- number_removed_rsd
  }

  ## Drop genes with missing values based on defined type filter "filt_thr"
  rna_se <- rna.filter_missing(filtered_rsd, type = filt_type, thr = filt_thr, min = filt_min)

  cat(paste0("Identified conditions:\n ", paste(str_c(unique(rna_se$condition), collapse = ", ")), "\n"))

  return(rna_se)
}

####____rna.filter_missing___####
rna.filter_missing <- function (se, type = c("complete", "condition", "fraction", NULL),
                                 thr = NULL, min = NULL)
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  type <- match.arg(type)
  if (any(!c("name", "ID") %in% colnames(SummarizedExperiment::rowData(se,
                                                                       use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)), "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in%
          colnames(SummarizedExperiment::colData(se)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)), "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  if (!is.null(type)){
    if (type == "complete") {
      keep <- !apply(SummarizedExperiment::assay(se), 1, function(x) any(is.na(x)))
      filtered <- se[keep, ]
    }
    if (type == "condition") {
      assertthat::assert_that(is.numeric(thr), length(thr) ==
                                1)
      max_repl <- max(SummarizedExperiment::colData(se)$replicate)
      if (thr < 0 | thr > max_repl) {
        stop("invalid filter threshold 'thr' applied",
             "\nRun filter() with a threshold ranging from 0 to ",
             max_repl)
      }
      filtered <- filter_missval(se, thr = thr)
    }
    if (type == "fraction") {
      assertthat::assert_that(is.numeric(min), length(min) ==
                                1)
      if (min < 0 | min > 1) {
        stop("invalid filter threshold 'min' applied",
             "\nRun filter() with a percent ranging from 0 to 1")
      }
      bin_data <- SummarizedExperiment::assay(se)
      idx <- is.na(SummarizedExperiment::assay(se))
      bin_data[!idx] <- 1
      bin_data[idx] <- 0
      keep <- bin_data %>% as.data.frame() %>% tibble::rownames_to_column() %>%
        gather(ID, value, -rowname) %>% group_by(rowname) %>%
        summarize(n = n(), valid = sum(value), frac = valid/n) %>%
        filter(frac >= min)
      filtered <- se[keep$rowname, ]
    }
  } else {
    filtered <- se
  }
  if( (nrow(SummarizedExperiment::assay(se)) - nrow(SummarizedExperiment::assay(filtered))) != 0 ){
    number_removed <- nrow(SummarizedExperiment::assay(se)) - nrow(SummarizedExperiment::assay(filtered))
    cat(paste0(number_removed, " out of ",
               nrow(SummarizedExperiment::assay(se)), " genes were removed from the dataset due to missing values.\n\n"))
    filtered@metadata$n.filtered <- number_removed
  }
  filtered@metadata$filt_type <- type
  filtered@metadata$n.pre_filt <- nrow(SummarizedExperiment::assay(se))
  return(filtered)
}

####____rna.workflow____####
rna.workflow <- function(se, # SummarizedExperiment, generated with read_prot().
                          imp_fun = c("zero", "man", "bpca", "knn", "QRILC", "MLE", "MinDet", # Method for imputing of missing values
                                      "MinProb", "min", "zero", "mixed", "nbavg", "SampMin"),
                          q = 0.01, # q value for imputing missing values with method "fun = 'MinProb'".
                          knn.rowmax = 0.5, # The maximum percent missing data allowed in any row (default 50%).
                          # For any rows with more than rowmax% missing are imputed using the overall mean per sample.
                          type = c("all", "control", "manual"), # Type of differential analysis to perform.
                          control = NULL, # Control condition; required if type = "control".
                          contrast = NULL, # Defined test for differential analysis "A_vs_B"; required if type = "manual".
                          controlGenes = NULL, # specifying those genes to use for size factor estimation (e.g. housekeeping or spike-in genes)
                          alpha = 0.05, # Significance threshold for adj. p values.
                          alpha_pathways = 0.1,
                          lfc = 1, # Relevance threshold for log2(fold change) values.
                          heatmap.show_all = TRUE, # Shall all samples be displayed in the heatmap or only the samples contained in the defined "contrast"?
                          # (only applicable for type = "manual")
                          heatmap.kmeans = F, # Shall the proteins be clustered in the heat map?
                          k = 6, # Number of protein clusters in heat map if kmeans = TRUE.
                          heatmap.col_limit = NA, # Define the breaks in the heat map legends.
                          heatmap.show_row_names = TRUE, # Show protein names in heat map?
                          heatmap.row_font_size = 6, # Font size of protein names if show_row_names = TRUE.
                          volcano.add_names = FALSE, # Display names next to symbols in volcano plot.
                          volcano.label_size = 2.5, # Size of labels in volcano plot if
                          volcano.adjusted = TRUE, # Display adjusted p-values on y axis of volcano plot?
                          plot = FALSE, # Shall plots be returned in the Plots pane?
                          export = FALSE, # Shall plots be exported as PDF and PNG files?
                          report = TRUE, # Shall a report (HTML and PDF) be created?
                          report.dir = NULL, # Folder name for created report (if report = TRUE)
                          pathway_enrichment = FALSE, # Perform pathway over-representation analysis for each tested contrast
                          pathway_kegg = FALSE, # Perform pathway over-representation analysis with gene sets in the KEGG database
                          kegg_organism = NULL, # Name of the organism in the KEGG database (if 'pathway_kegg = TRUE')
                          custom_pathways = NULL) # Dataframe providing custom pathway annotations
{
  # Show error if inputs are not the required classes
  if(is.integer(alpha)) alpha <- as.numeric(alpha)
  if(is.integer(lfc)) lfc <- as.numeric(lfc)
  assertthat::assert_that(is.character(imp_fun),
                          is.character(type),
                          is.numeric(alpha),
                          length(alpha) == 1,
                          is.numeric(lfc),
                          length(lfc) == 1)

  # Show error if inputs are not valid
  if (!type %in% c("all", "control", "manual")) {
    stop("run workflow_proteomics() with a valid type",
         "\nValid types are: 'all', 'control' and 'manual'.",
         call. = FALSE)
  }
  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
  }
  # Impute missing values
  if(length(imp_fun) == 1){
    if (imp_fun == "MinProb"){
        prot_imp <- rna.impute(se, fun = imp_fun, q = q)
    } else if (imp_fun == "knn"){
        prot_imp <- rna.impute(se, fun = imp_fun, rowmax = knn.rowmax)
    } else {
        prot_imp <- rna.impute(se, fun = imp_fun)
    }
  } else {
    prot_imp <- rna.impute(se, fun = "zero")
  }

  # Create DESeqDataSet
  ddsSE <- DESeq2::DESeqDataSet(prot_imp, design = ~ condition)

  # Pre-filtering:  keep only rows that have at least 10 reads total
  keep <- BiocGenerics::rowSums(BiocGenerics::counts(ddsSE)) >= 10
  ddsSE <- ddsSE[keep,]


  # Perform PCA Analysis
  rna_pca <- prot.pca(SummarizedExperiment::assay(ddsSE))

  # Differential expression analysis
  rna_diff <- DESeq2::DESeq(ddsSE, sfType = "ratio")














  function (object, test = c("Wald", "LRT"), fitType = c("parametric",
                                                         "local", "mean", "glmGamPoi"), sfType = c("ratio", "poscounts",
                                                                                                   "iterate"), betaPrior, full = design(object), reduced, quiet = FALSE,
            minReplicatesForReplace = 7, modelMatrixType, useT = FALSE,
            minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5, parallel = FALSE,
            BPPARAM = bpparam())
  {
    stopifnot(is(object, "DESeqDataSet"))
    test <- match.arg(test, choices = c("Wald", "LRT"))
    fitType <- match.arg(fitType, choices = c("parametric", "local",
                                              "mean", "glmGamPoi"))
    dispersionEstimator <- if (fitType == "glmGamPoi") {
      "glmGamPoi"
    }
    else {
      "DESeq2"
    }
    if (fitType == "glmGamPoi") {
      minReplicatesForReplace <- Inf
      if (parallel) {
        warning("parallelization of DESeq() is not implemented for fitType='glmGamPoi'")
      }
    }
    sfType <- match.arg(sfType, choices = c("ratio", "poscounts",
                                            "iterate"))
    stopifnot(is.logical(quiet))
    stopifnot(is.numeric(minReplicatesForReplace))
    stopifnot(is.logical(parallel))
    modelAsFormula <- !is.matrix(full) & is(design(object), "formula")
    if (missing(betaPrior)) {
      betaPrior <- FALSE
    }
    else {
      stopifnot(is.logical(betaPrior))
    }
    object <- sanitizeRowRanges(object)
    if (test == "LRT") {
      if (missing(reduced)) {
        stop("likelihood ratio test requires a 'reduced' design, see ?DESeq")
      }
      if (betaPrior) {
        stop("test='LRT' does not support use of LFC shrinkage, use betaPrior=FALSE")
      }
      if (!missing(modelMatrixType) && modelMatrixType == "expanded") {
        stop("test='LRT' does not support use of expanded model matrix")
      }
      if (is.matrix(full) | is.matrix(reduced)) {
        if (!(is.matrix(full) & is.matrix(reduced))) {
          stop("if one of 'full' and 'reduced' is a matrix, the other must be also a matrix")
        }
      }
      if (modelAsFormula) {
        checkLRT(full, reduced)
      }
      else {
        checkFullRank(full)
        checkFullRank(reduced)
        if (ncol(full) <= ncol(reduced)) {
          stop("the number of columns of 'full' should be more than the number of columns of 'reduced'")
        }
      }
    }
    if (test == "Wald" & !missing(reduced)) {
      stop("'reduced' ignored when test='Wald'")
    }
    if (dispersionEstimator == "glmGamPoi" && test == "Wald") {
      warning("glmGamPoi dispersion estimator should be used in combination with a LRT and not a Wald test.",
              call. = FALSE)
    }
    if (modelAsFormula) {
      designAndArgChecker(object, betaPrior)
      if (design(object) == formula(~1)) {
        warning("the design is ~ 1 (just an intercept). is this intended?")
      }
      if (full != design(object)) {
        stop("'full' specified as formula should equal design(object)")
      }
      modelMatrix <- NULL
    }
    else {
      if (!quiet)
        message("using supplied model matrix")
      if (betaPrior == TRUE) {
        stop("betaPrior=TRUE is not supported for user-provided model matrices")
      }
      checkFullRank(full)
      modelMatrix <- full
    }
    attr(object, "betaPrior") <- betaPrior
    stopifnot(length(parallel) == 1 & is.logical(parallel))
    if (!is.null(sizeFactors(object)) || !is.null(normalizationFactors(object))) {
      if (!quiet) {
        if (!is.null(normalizationFactors(object))) {
          message("using pre-existing normalization factors")
        }
        else {
          message("using pre-existing size factors")
        }
      }
    }
    else {
      if (!quiet)
        message("estimating size factors")
        object <- estimateSizeFactors(object, type = sfType,
                                      quiet = quiet, controlGenes = c(3, 5))
    }
    if (!parallel) {
      if (!quiet)
        message("estimating dispersions")
      object <- estimateDispersions(object, fitType = fitType,
                                    quiet = quiet, modelMatrix = modelMatrix, minmu = minmu)
      if (!quiet)
        message("fitting model and testing")
      if (test == "Wald") {
        object <- nbinomWaldTest(object, betaPrior = betaPrior,
                                 quiet = quiet, modelMatrix = modelMatrix, modelMatrixType = modelMatrixType,
                                 useT = useT, minmu = minmu)
      }
      else if (test == "LRT") {
        object <- nbinomLRT(object, full = full, reduced = reduced,
                            quiet = quiet, minmu = minmu, type = dispersionEstimator)
      }
    }
    else if (parallel) {
      if (!missing(modelMatrixType)) {
        if (betaPrior)
          stopifnot(modelMatrixType == "expanded")
      }
      object <- DESeqParallel(object, test = test, fitType = fitType,
                              betaPrior = betaPrior, full = full, reduced = reduced,
                              quiet = quiet, modelMatrix = modelMatrix, useT = useT,
                              minmu = minmu, BPPARAM = BPPARAM)
    }
    sufficientReps <- any(nOrMoreInCell(attr(object, "modelMatrix"),
                                        minReplicatesForReplace))
    if (sufficientReps) {
      object <- refitWithoutOutliers(object, test = test, betaPrior = betaPrior,
                                     full = full, reduced = reduced, quiet = quiet, minReplicatesForReplace = minReplicatesForReplace,
                                     modelMatrix = modelMatrix, modelMatrixType = modelMatrixType)
    }
    metadata(object)[["version"]] <- packageVersion("DESeq2")
    object
  }






























  # Denote significantly differentially expressed proteins
  prot_dep <- prot.add_rejections(prot_diff, alpha = alpha, lfc = lfc)
  contrasts <- SummarizedExperiment::rowData(prot_dep) %>%
    data.frame(check.names = FALSE) %>%
    select(ends_with("_diff")) %>%
    colnames() %>% str_replace_all("_diff", "")
  # Generate a results table
  results <- prot.get_results(prot_dep)
  n_significant <- results %>% filter(significant) %>% nrow()

  message(paste0(n_significant,
                 " proteins were found to be differentially expressed with ",
                 expression(alpha),
                 " = ", alpha,
                 " and |log2(fold change)| > ",
                 lfc, "."))
  if (pathway_enrichment) {
    ls.significant_df <- list()
    ls.significant_up <- list()
    ls.significant_dn <- list()

    for (i in 1:length(contrasts)) {
      ls.significant_df[[i]] <-
        SummarizedExperiment::rowData(prot_dep[SummarizedExperiment::rowData(prot_dep)[[paste0(contrasts[i], "_significant")]],]) %>% data.frame()
      ls.significant_up[[i]] <-
        ls.significant_df[[i]][ls.significant_df[[i]][paste0(contrasts[i], "_diff")] > 0,]
      ls.significant_dn[[i]] <-
        ls.significant_df[[i]][ls.significant_df[[i]][paste0(contrasts[i], "_diff")] < 0,]
    }
    if(!pathway_kegg && is.null(custom_pathways)) {
      stop(
        "Cannot perform custom pathway over-representation analysis without a table of pathways and corresponding genes.\nPlease provide a dataframe containing 'Pathway' and 'Accession' columns in the 'custom_pathways =' argument. Alternatively, choose 'pathway_kegg = TRUE' and a valid KEGG organism identifier in the 'kegg_organism = ' argument."
      )
    }
    if (pathway_kegg) {
      if (is.null(kegg_organism)) {
        stop(
          "Cannot perform KEGG pathway over-representation analysis without specifying a valid KEGG organism id in the 'kegg_organism' argument."
        )
      } else {

        ls.pora_kegg_up <- rep(list(0), length(contrasts))
        kegg_pathway_up <- function(x) {
          prot.pathway_enrich(
            gene = ls.significant_up[[x]]$ID,
            organism = kegg_organism,
            keyType = 'kegg',
            pvalueCutoff = alpha_pathways,
            pAdjustMethod = "BH",
            minGSSize = 2)
        }
        ls.pora_kegg_up <- suppressMessages(lapply(1:length(contrasts), kegg_pathway_up))

        ls.pora_kegg_dn <- rep(list(0), length(contrasts))
        kegg_pathway_dn <- function(x) {
          prot.pathway_enrich(
            gene = ls.significant_dn[[x]]$ID,
            organism = kegg_organism,
            keyType = 'kegg',
            pvalueCutoff = alpha_pathways,
            pAdjustMethod = "BH",
            minGSSize = 2)
        }
        ls.pora_kegg_dn <- suppressMessages(lapply(1:length(contrasts), kegg_pathway_dn))

        for (i in 1:length(contrasts)) {

          if(is.null(nrow(ls.pora_kegg_up[[i]]))){
            cat(paste0("No significantly upregulated KEGG pathways found for contrast:\n", contrasts[i], "\n"))
          }
          if(is.null(nrow(ls.pora_kegg_dn[[i]]))){
            cat(paste0("No significantly downregulated KEGG pathways found for contrast:\n", contrasts[i], "\n"))
          }

        }
        names(ls.pora_kegg_up) <- contrasts
        names(ls.pora_kegg_dn) <- contrasts
      }
    }
    if (!is.null(custom_pathways)) {

      ls.pora_custom_up <- rep(list(0), length(contrasts))
      custom_pathway_up <- function(x) {
        prot.pathway_enrich(
          gene = ls.significant_up[[x]]$ID,
          pvalueCutoff = alpha_pathways,
          pAdjustMethod = "BH",
          custom_gene_sets = T,
          custom_pathways = custom_pathways,
          minGSSize = 2)
      }
      ls.pora_custom_up <- suppressMessages(lapply(1:length(contrasts), custom_pathway_up))

      ls.pora_custom_dn <- rep(list(0), length(contrasts))
      custom_pathway_dn <- function(x) {
        prot.pathway_enrich(
          gene = ls.significant_dn[[x]]$ID,
          pvalueCutoff = alpha_pathways,
          pAdjustMethod = "BH",
          custom_gene_sets = T,
          custom_pathways = custom_pathways,
          minGSSize = 2)
      }
      ls.pora_custom_dn <- suppressMessages(lapply(1:length(contrasts), custom_pathway_dn))

      for (i in 1:length(contrasts)) {
        if(is.null(nrow(ls.pora_custom_up[[i]]))){
          cat(paste0("No significantly upregulated custom pathways found for contrast:\n", contrasts[i], "\n"))
        }
        if(is.null(nrow(ls.pora_custom_dn[[i]]))){
          cat(paste0("No significantly downregulated custom pathways found for contrast:\n", contrasts[i], "\n"))
        }
      }
      names(ls.pora_custom_up) <- contrasts
      names(ls.pora_custom_dn) <- contrasts
    }
  }

  suppress_ggrepel <- function(w) {
    if (any(grepl("ggrepel", w)))
      invokeRestart("muffleWarning")
  }

  if (export == TRUE |
      plot == TRUE) {
    if (export == TRUE){
      message(paste0("Rendering and exporting figures to:\n",
                     getwd(), "/Plots."))
    }
    suppressMessages(
      prot.boxplot_intensity(se, prot_norm, prot_imp, plot = plot, export = export)
    )
    try(suppressMessages(
      prot.plot_missval(prot_norm, plot = plot, export = export)
    ))
    try(suppressMessages(
      prot.plot_detect(prot_norm, basesize = 10, plot = plot, export = export)
    ))
    suppressMessages(
      prot.plot_imputation(se, prot_norm, prot_imp, plot = plot, export = export, basesize = 12)
    )
    suppressMessages(
      suppressWarnings(
        prot.plot_screeplot(prot_pca, axisLabSize = 18, titleLabSize = 22, plot = plot, export = export)
      ) )
    withCallingHandlers(suppressMessages(
      prot.plot_loadings(prot_pca,labSize = 3, plot = plot, export = export)
    ) , warning = suppress_ggrepel)
    suppressMessages(
      prot.plot_pca(prot_imp,x = 1, y = 2, point_size = 4, basesize = 14, title = "PC Scores - PC2 vs. PC1",
                    plot = plot,export = export)
    )
    suppressMessages(
      prot.plot_pca(prot_imp,x = 1, y = 3, point_size = 4, basesize = 14, title = "PC Scores - PC3 vs. PC1",
                    plot = plot,export = export)
    )
    suppressMessages(
      prot.plot_heatmap(prot_dep, type = "centered", kmeans = heatmap.kmeans, show_all = heatmap.show_all, contrast = contrast,
                        k = k, col_limit = heatmap.col_limit,show_row_names = heatmap.show_row_names,
                        row_font_size = heatmap.row_font_size, indicate = c("condition"),
                        plot = plot, export = export)
    )
    suppressMessages(
      prot.plot_heatmap(prot_dep, type = "contrast", contrast = contrast,  kmeans = heatmap.kmeans, k = k, col_limit = heatmap.col_limit,
                        show_row_names = heatmap.show_row_names, row_font_size = heatmap.row_font_size,
                        plot = plot, export = export)
    )

    for (i in 1:length(contrasts)){
      suppressMessages(
        suppressWarnings(
          prot.plot_volcano(prot_dep, contrast = contrasts[i],
                            add_names = volcano.add_names, label_size = volcano.label_size, adjusted =  volcano.adjusted,
                            plot = plot, export = export, lfc = lfc, alpha = alpha)
        ) )
    }
    if (pathway_kegg) {
      for (i in 1:length(contrasts)) {
        if(!(nrow(as.data.frame(ls.pora_kegg_up[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(ls.pora_kegg_up[[i]], title = paste0("Upregulated pathways", " - KEGG"),
                                   subtitle =  str_replace(contrasts[i], "_vs_", " vs. "), plot = plot, export = export, kegg = TRUE)
            ) )
        }
        if(!(nrow(as.data.frame(ls.pora_kegg_dn[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(ls.pora_kegg_dn[[i]], title = paste0("Downregulated pathways", " - KEGG"),
                                   subtitle = str_replace(contrasts[i], "_vs_", " vs. "), plot = plot, export = export, kegg = TRUE)
            ) )
        }
      }
    }
    if(!is.null(custom_pathways)){
      for (i in 1:length(contrasts)) {
        if(!(nrow(as.data.frame(ls.pora_custom_up[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(ls.pora_custom_up[[i]], title = "Upregulated pathways",
                                   subtitle =  str_replace(contrasts[i], "_vs_", " vs. "), plot = plot, export = export, kegg = FALSE)
            ) )
        }
        if(!(nrow(as.data.frame(ls.pora_custom_dn[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(ls.pora_custom_dn[[i]], title = "Downregulated pathways",
                                   subtitle = str_replace(contrasts[i], "_vs_", " vs. "), plot = plot, export = export, kegg = FALSE)
            ) )
        }
      }
    }
  }





  param <- data.frame(type, alpha, lfc, check.names = FALSE)
  results <- list(data = SummarizedExperiment::rowData(se), se = se, norm = prot_norm,
                  imputed = prot_imp, pca = prot_pca, diff = prot_diff, dep = prot_dep,
                  results = results, param = param)
  if (pathway_enrichment == T && pathway_kegg) {
    results <- c(results, pora_kegg_up = list(ls.pora_kegg_up), pora_kegg_dn = list(ls.pora_kegg_dn))
  }
  if (pathway_enrichment == T && !is.null(custom_pathways)) {
    results <- c(results, pora_custom_up = list(ls.pora_custom_up), pora_custom_dn = list(ls.pora_custom_dn))
  }

  if(report == TRUE){
    prot.report(results,
                volcano.adjusted = volcano.adjusted,
                pathway_enrichment = pathway_enrichment,
                heatmap.show_all = heatmap.show_all,
                heatmap.kmeans = heatmap.kmeans,
                k = k,
                report.dir = report.dir)
  }
  return(results)
}

####____rna.impute____####
rna.impute <- function (se, fun = c("zero", "bpca", "knn", "QRILC",
                                     "MLE", "MinDet", "MinProb", "man",
                                     "min", "zero", "mixed", "nbavg", "SampMin"), ...)
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(fun))
  fun <- match.arg(fun)
  if (any(!c("name", "ID") %in% colnames(SummarizedExperiment::rowData(se,
                                                                       use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)), "'\nRun prot.make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if (!any(is.na(SummarizedExperiment::assay(se)))) {
    warning("No missing values in '", deparse(substitute(se)),
            "'. ", "Returning the unchanged object.",
            call. = FALSE)
    return(se)
  }
  SummarizedExperiment::rowData(se)$imputed <- apply(is.na(SummarizedExperiment::assay(se)), 1, any)
  SummarizedExperiment::rowData(se)$num_NAs <- rowSums(is.na(SummarizedExperiment::assay(se)))
  if(fun != "zero" && length(fun) == 1){
      if (fun == "man") {
        se <- manual_impute(se, ...)
      } else if (fun == "SampMin"){
        SummarizedExperiment::assay(se) <- SummarizedExperiment::assay(se) %>% data.frame() %>%
          mutate_if(is.numeric, function(x) ifelse(is.na(x), min(x, na.rm = T), x)) %>%
          as.matrix()
        cat("Missing values were replaced with the minimum value in each sample.")
      } else {
        MSnSet_data <- methods::as(se, "MSnSet")
        MSnSet_imputed <- MSnbase::impute(MSnSet_data, method = fun,
                                          ...)
        SummarizedExperiment::assay(se) <- MSnbase::exprs(MSnSet_imputed)
      }
  } else {
    SummarizedExperiment::assay(se) <- SummarizedExperiment::assay(se) %>% data.frame() %>%
      mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x)) %>%
      as.matrix()
    cat("Missing values were replaced with 0.\n\n")
  }
  se@metadata["imp_fun"] <- fun
  return(se)
}

####____rna.make_ddsSE____####
rna.make_se <- function(genes_unique, columns, expdesign)
{
  assertthat::assert_that(is.data.frame(genes_unique), is.integer(columns),
                          is.data.frame(expdesign))
  if (any(!c("name", "ID") %in% colnames(genes_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)), "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns",
         "are not present in the experimental design", call. = FALSE)
  }
  if (any(!apply(genes_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric", "\nRun make_se_parse() with the appropriate columns as argument",
         call. = FALSE)
  }
  if (tibble::is_tibble(genes_unique))
    genes_unique <- as.data.frame(genes_unique)
  if (tibble::is_tibble(expdesign))
    expdesign <- as.data.frame(expdesign)
  rownames(genes_unique) <- genes_unique$name
  raw <- genes_unique[, columns]
  raw[is.na(raw)] <- 0
  expdesign <- mutate(expdesign, condition = make.names(condition)) %>%
    tidyr::unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID
  matched <- match(make.names(delete_prefix(expdesign$label)),
                   make.names(delete_prefix(colnames(raw))))
  if (any(is.na(matched))) {
    stop("None of the labels in the experimental design match ",
         "with column names in 'genes_unique'", "\nRun make_se() with the correct labels in the experimental design",
         "and/or correct columns specification")
  }
  colnames(raw)[matched] <- expdesign$ID
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
  row_data <- genes_unique[, -columns]
  rownames(row_data) <- row_data$name
  mat <- as.matrix(sapply(raw, as.integer))
  rownames(mat) <- rownames(raw)
  mat[mat==0] <- NA
  se <- SummarizedExperiment::SummarizedExperiment(assays = mat, colData = expdesign,
                                                   rowData = row_data)
  return(se)
}


rna.DESeq <- function (object, test = c("Wald", "LRT"), fitType = c("parametric", "local", "mean", "glmGamPoi"),
                       sfType = c("ratio", "poscounts", "iterate"), controlGenes = NULL, betaPrior, full = design(object), reduced, quiet = FALSE,
                       minReplicatesForReplace = 7, modelMatrixType, useT = FALSE,
                       minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5, parallel = FALSE,
                       BPPARAM = bpparam())
{
  stopifnot(is(object, "DESeqDataSet"))
  test <- match.arg(test, choices = c("Wald", "LRT"))
  fitType <- match.arg(fitType, choices = c("parametric", "local",
                                            "mean", "glmGamPoi"))
  dispersionEstimator <- if (fitType == "glmGamPoi") {
    "glmGamPoi"
  }
  else {
    "DESeq2"
  }
  if (fitType == "glmGamPoi") {
    minReplicatesForReplace <- Inf
    if (parallel) {
      warning("parallelization of DESeq() is not implemented for fitType='glmGamPoi'")
    }
  }
  sfType <- match.arg(sfType, choices = c("ratio", "poscounts",
                                          "iterate"))
  if(!is.null(controlGenes)){
    if(is.numeric(controlGenes)){
      controlGenes <- controlGenes
    } else if (is.character(controlGenes)){
      controlGenes <- match(controlGenes, rownames(assay(object)))
    }
  }
  stopifnot(is.logical(quiet))
  stopifnot(is.numeric(minReplicatesForReplace))
  stopifnot(is.logical(parallel))
  modelAsFormula <- !is.matrix(full) & is(design(object), "formula")
  if (missing(betaPrior)) {
    betaPrior <- FALSE
  }
  else {
    stopifnot(is.logical(betaPrior))
  }
  object <- sanitizeRowRanges(object)
  if (test == "LRT") {
    if (missing(reduced)) {
      stop("likelihood ratio test requires a 'reduced' design, see ?DESeq")
    }
    if (betaPrior) {
      stop("test='LRT' does not support use of LFC shrinkage, use betaPrior=FALSE")
    }
    if (!missing(modelMatrixType) && modelMatrixType == "expanded") {
      stop("test='LRT' does not support use of expanded model matrix")
    }
    if (is.matrix(full) | is.matrix(reduced)) {
      if (!(is.matrix(full) & is.matrix(reduced))) {
        stop("if one of 'full' and 'reduced' is a matrix, the other must be also a matrix")
      }
    }
    if (modelAsFormula) {
      checkLRT(full, reduced)
    }
    else {
      checkFullRank(full)
      checkFullRank(reduced)
      if (ncol(full) <= ncol(reduced)) {
        stop("the number of columns of 'full' should be more than the number of columns of 'reduced'")
      }
    }
  }
  if (test == "Wald" & !missing(reduced)) {
    stop("'reduced' ignored when test='Wald'")
  }
  if (dispersionEstimator == "glmGamPoi" && test == "Wald") {
    warning("glmGamPoi dispersion estimator should be used in combination with a LRT and not a Wald test.",
            call. = FALSE)
  }
  if (modelAsFormula) {
    designAndArgChecker(object, betaPrior)
    if (design(object) == formula(~1)) {
      warning("the design is ~ 1 (just an intercept). is this intended?")
    }
    if (full != design(object)) {
      stop("'full' specified as formula should equal design(object)")
    }
    modelMatrix <- NULL
  }
  else {
    if (!quiet)
      message("using supplied model matrix")
    if (betaPrior == TRUE) {
      stop("betaPrior=TRUE is not supported for user-provided model matrices")
    }
    checkFullRank(full)
    modelMatrix <- full
  }
  attr(object, "betaPrior") <- betaPrior
  stopifnot(length(parallel) == 1 & is.logical(parallel))
  if (!is.null(sizeFactors(object)) || !is.null(normalizationFactors(object))) {
    if (!quiet) {
      if (!is.null(normalizationFactors(object))) {
        message("using pre-existing normalization factors")
      }
      else {
        message("using pre-existing size factors")
      }
    }
  }
  else {
    if (!quiet)
      message("estimating size factors")
    object <- estimateSizeFactors(object, type = sfType, controlGenes = controlGenes, quiet = quiet)
  }
  if (!parallel) {
    if (!quiet)
      message("estimating dispersions")
    object <- estimateDispersions(object, fitType = fitType,
                                  quiet = quiet, modelMatrix = modelMatrix, minmu = minmu)
    if (!quiet)
      message("fitting model and testing")
    if (test == "Wald") {
      object <- nbinomWaldTest(object, betaPrior = betaPrior,
                               quiet = quiet, modelMatrix = modelMatrix, modelMatrixType = modelMatrixType,
                               useT = useT, minmu = minmu)
    }
    else if (test == "LRT") {
      object <- nbinomLRT(object, full = full, reduced = reduced,
                          quiet = quiet, minmu = minmu, type = dispersionEstimator)
    }
  }
  else if (parallel) {
    if (!missing(modelMatrixType)) {
      if (betaPrior)
        stopifnot(modelMatrixType == "expanded")
    }
    object <- DESeqParallel(object, test = test, fitType = fitType,
                            betaPrior = betaPrior, full = full, reduced = reduced,
                            quiet = quiet, modelMatrix = modelMatrix, useT = useT,
                            minmu = minmu, BPPARAM = BPPARAM)
  }
  sufficientReps <- any(nOrMoreInCell(attr(object, "modelMatrix"),
                                      minReplicatesForReplace))
  if (sufficientReps) {
    object <- refitWithoutOutliers(object, test = test, betaPrior = betaPrior,
                                   full = full, reduced = reduced, quiet = quiet, minReplicatesForReplace = minReplicatesForReplace,
                                   modelMatrix = modelMatrix, modelMatrixType = modelMatrixType)
  }
  metadata(object)[["version"]] <- packageVersion("DESeq2")
  object
}
