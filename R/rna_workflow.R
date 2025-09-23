#' @title RNA sequencing workflow
#'
#' @description This function performs a complete RNA sequencing workflow, including imputation of missing values, normalization,
#' principal component analysis, differential expression analysis, and pathway analysis. The function also provides
#' several options for plotting, exporting plots, and creating a report.
#'
#' @param se A SummarizedExperiment object, generated with read_prot().
#' @param imp_fun (Character string)  Function used for data imputation. "SampMin", "man", "bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb", "min", "zero", "mixed", or "nbavg". See (\code{\link{rna.impute}}) for details.
#' @param trans_method (Character string) Transformation method for the data. Options are "rlog" (for \code{\link[DESeq2]{rlog}}) or "vst" (for \code{\link[DESeq2]{vst}}).
#' @param q (Numeric) q value for imputing missing values with method \code{imp_fun = 'MinProb'}.
#' @param knn.rowmax (Numeric) The maximum percent missing data allowed in any row for \code{imp_fun = 'knn'}. Default: 0.5.
#' @param type (Character string) Type of differential analysis to perform. "all" (contrast each condition with every other condition), "control" (contrast each condition to a defined control condition), "manual" (manually define selected conditions).
#' @param design Formula for the design matrix.
#' @param size.factors Optional: Manually define size factors for normalization.
#' @param altHypothesis  Specify those genes you are interested in finding. The test provides p values for the null hypothesis, the complement of the set defined by altHypothesis. For further details, see \code{\link[DESeq2]{results}}.
#' @param control Control condition; required if type = "control".
#' @param contrast (String or vector of strings) Defined test(s) for differential analysis in the form "A_vs_B"; required if type = "manual".
#' @param controlGenes Specifying those genes to use for size factor estimation (e.g. housekeeping or spike-in genes).
#' @param pAdjustMethod Method for adjusting p values. Available options are "IHW" (Independent Hypothesis Weighting),"BH" (Benjamini-Hochberg).
#' @param alpha Significance threshold for adjusted p values.
#' @param alpha.independent Adjusted p value threshold for independent filtering or NULL. If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
#' @param alpha_pathways Significance threshold for pathway analysis.
#' @param lfcShrink Use shrinkage to calculate log2 fold change values.
#' @param shrink.method Method for shrinkage. Available options are "apeglm", "ashr", "normal". See \code{\link[DESeq2]{lfcShrink}} for details.
#' @param lfc Relevance threshold for absolute log2(fold change) values. Used to filter unshrunken lfc values
#' or in shrinkage method "apeglm" or "normal".
#' @param heatmap.show_all Shall all samples be displayed in the heatmap or only the samples contained in the defined
#' "contrast"? (only applicable for type = "manual")
#' @param heatmap.kmeans Shall the proteins be clustered in the heat map?
#' @param k Number of protein clusters in heat map if kmeans = TRUE.
#' @param heatmap.col_limit Define the breaks in the heat map legends.
#' @param heatmap.show_row_names Show protein names in heat map?
#' @param heatmap.row_font_size Font size of protein names if show_row_names = TRUE.
#' @param volcano.add_names Display names next to symbols in volcano plot.
#' @param volcano.label_size Size of labels in volcano plot.
#' @param volcano.adjusted Display adjusted p-values on y axis of volcano plot?
#' @param plot Shall plots be returned in the Plots pane?
#' @param plot_volcano (Logical) Generate volcano plots for each defined contrast (\code{TRUE}) or not \code{FALSE}).
#' @param export Shall plots be exported as PDF and PNG files?
#' @param report Shall a report (HTML and PDF) be created?
#' @param report.dir Folder name for created report (if report = TRUE)
#' @param pathway_enrichment Perform pathway over-representation analysis for each tested contrast
#' @param pathway_kegg Perform pathway over-representation analysis with gene sets in the KEGG database
#' @param kegg_organism Name of the organism in the KEGG database (if 'pathway_kegg = TRUE')
#' @param custom_pathways Dataframe providing custom pathway annotations.
#' @param sample_covariates Optional data.frame of sample-level covariates (rownames must match sample IDs / column names of the expression matrix). Numeric columns will be used for correlation analysis.
#' @param correlate_genes Logical; if TRUE compute per-gene correlations vs provided numeric sample covariates on transformed (rlog/vst) expression values.
#' @param correlate_pathways Logical; if TRUE compute per-pathway correlations (mean expression per pathway) vs covariates. Requires `custom_pathways` as gene set definition (data.frame with Pathway/Accession or GMT filepath).
#' @param pathway_min_genes Minimum number of genes required in a pathway to compute a pathway score (default 5).
#' @param correlation_method Correlation method for covariate analyses: "spearman" (default) or "pearson".
#' @param correlation_p_adjust Multiple testing adjustment method applied within each covariate (default "BH").
#' @param gsea Perform GSEA (Gene Set Enrichment Analysis) for each defined contrast.
#' @param gsea_gmt Path to a GMT file containing gene sets for GSEA.
#' @param gsea_pAdjustMethod Method for adjusting p-values in GSEA. Options are "BH" (Benjamini-Hochberg) or "none".
#' @param gsea_pvalueCutoff P-value cutoff for GSEA results.
#' @param gsea_minGSSize Minimum size of gene sets to be considered in GSEA.
#' @param gsea_maxGSSize Maximum size of gene sets to be considered in GSEA.
#' @param quiet Suppress messages and warnings.
#' @param allow_no_replicates Allow analysis of conditions with no replicates (default: TRUE).
#' @param de_if_no_reps Strategy for differential expression when no replicates per condition. One of "none" or "edgeR_fixed" (default: "edgeR_fixed").
#' @param assumed_dispersion Fixed dispersion value used by edgeR when `de_if_no_reps = "edgeR_fixed"` (default: 0.1).
#'
#' @return The function returns a SummarizedExperiment object with added columns for log2 fold change, p-values and adjusted p-values for each comparison.
#' It also includes a column for significant genes for each comparison and a column for significant genes overall.
#' Additionally, the function generates various plots and a report (if specified).
#' @export
rna.workflow <- function(se, # SummarizedExperiment, generated with read_prot().
                         imp_fun = c("zero", "man", "bpca", "knn", "QRILC", "MLE", "MinDet",
                                     "MinProb", "min", "zero", "mixed", "nbavg", "SampMin"),
                         trans_method = c("vst", "rlog"),
                         q = 0.01,
                         knn.rowmax = 0.5,
                         type = c("all", "control", "manual"),
                         design = "~ 0 + condition",
                         size.factors = NULL,
                         altHypothesis = c("greaterAbs", "lessAbs", "greater", "less"),
                         control = NULL,
                         contrast = NULL,
                         controlGenes = NULL,
                         pAdjustMethod = c("IHW","BH"),
                         alpha = 0.05,
                         alpha.independent = 0.1,
                         alpha_pathways = 0.1,
                         lfcShrink = TRUE,
                         shrink.method = c("apeglm", "ashr", "normal"),
                         lfc = 2,
                         heatmap.show_all = TRUE,
                         heatmap.kmeans = F,
                         k = 6,
                         heatmap.col_limit = NA,
                         heatmap.show_row_names = TRUE,
                         heatmap.row_font_size = 6,
                         volcano.add_names = FALSE,
                         volcano.label_size = 2.5,
                         volcano.adjusted = TRUE,
                         plot = FALSE,
                         plot_volcano = FALSE, # Shall volcano plots be generated?
                         export = FALSE,
                         report = TRUE,
                         report.dir = NULL,
                         pathway_enrichment = FALSE,
                         pathway_kegg = FALSE,
                         kegg_organism = NULL,
                         custom_pathways = NULL,
                         # Sample covariate correlation parameters
                         sample_covariates = NULL,              # data.frame with rownames = sample IDs, numeric columns are used
                         correlate_genes = FALSE,               # compute gene-level correlations vs covariates
                         correlate_pathways = FALSE,            # compute pathway-level correlations vs covariates (requires custom_pathways)
                         pathway_min_genes = 5,                 # min genes for pathway score
                         correlation_method = c("spearman","pearson"),
                         correlation_p_adjust = "BH",
                         # GSEA parameters
                         gsea = FALSE,
                         gsea_gmt = NULL,
                         gsea_pAdjustMethod = c("BH", "none"),
                         gsea_pvalueCutoff = 0.05,
                         gsea_minGSSize = 5,
                         gsea_maxGSSize = 500,
                         quiet = FALSE,
                         allow_no_replicates = TRUE,
                         de_if_no_reps = c("edgeR_fixed","none"),
                         assumed_dispersion = 0.1
)
{
  # Show error if inputs are not the required classes
  imp_fun <- match.arg(imp_fun)
  shrink.method <- match.arg(shrink.method)
  pAdjustMethod <- match.arg(pAdjustMethod)
  gsea_pAdjustMethod <- match.arg(gsea_pAdjustMethod)
  de_if_no_reps <- match.arg(de_if_no_reps)
  correlation_method <- match.arg(correlation_method)
  if(is.integer(alpha)) alpha <- as.numeric(alpha)
  if(is.integer(lfc)) lfc <- as.numeric(lfc)
  assertthat::assert_that(is.character(imp_fun),
                          is.character(type),
                          is.numeric(alpha),
                          length(alpha) == 1,
                          is.numeric(lfc),
                          length(lfc) == 1
  )
  # Assert that trans_method is one of "rlog" or "vst"
  if (!trans_method %in% c("rlog", "vst")) {
    stop("run rna.workflow() with a valid 'trans_method'",
         "\nValid trans_method values are: 'rlog' and 'vst'.",
         call. = FALSE)
  }

  # Show error if inputs are not valid
  if (!type %in% c("all", "control", "manual")) {
    stop("run rna.workflow() with a valid 'type'",
         "\nValid types are: 'all', 'control' and 'manual'.",
         call. = FALSE)
  }
  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
  }
  # Pre-filtering:  keep only rows that have at least 3 samples with a count of 10 or more
  message("Pre-filtering:  keep only rows that have at least 3 samples with a count of 10 or more")
  message("Number of genes before pre-filtering: ", nrow(assay(se)))
  keep <- MatrixGenerics::rowSums(assay(se) >= 10, na.rm = TRUE) >= 3
  se <- se[keep, ]
  message("Number of genes after pre-filtering for counts ≥ 10 in ≥ 3 samples: ", nrow(assay(se)))

  # Impute missing values
  if (imp_fun == "MinProb"){
    se_imp <- rna.impute(se, fun = imp_fun, q = q)
  } else if (imp_fun == "knn"){
    se_imp <- rna.impute(se, fun = imp_fun, rowmax = knn.rowmax)
  } else {
    se_imp <- rna.impute(se, fun = imp_fun)
  }

  # Second Pre-filtering:  keep only rows that have at least 3 samples with a CPM of 1 or more
  dge   <- edgeR::DGEList(counts = assay(se_imp))
  cpm_m <- edgeR::cpm(dge)
  keep2 <- rowSums(cpm_m >= 1) >= 3        # CPM ≥ 1 in ≥ 3 samples
  se_imp    <- se_imp[keep2, ]
  message("Number of genes after pre-filtering for CPM ≥ 1 in ≥ 3 samples: ", nrow(assay(se_imp)))
  # Create DESeqDataSet and define contrasts
  if(!quiet) message(">> Creating DESeqDataSet object <<")
  conditions <- as.character(unique(se$condition))
  se_imp$condition <- as.factor(se_imp$condition)
  if(!quiet){
    ddsSE <- DESeq2::DESeqDataSet(se_imp, design = as.formula(design))
  } else {
    ddsSE <- suppressMessages( DESeq2::DESeqDataSet(se_imp, design = as.formula(design)) )
  }
  if (type == "control") {
    if (is.null(control))
      stop("run test_diff(type = 'control') with a 'control' argument")
    cntrst <- paste(conditions[!conditions %in% control],
                    control, sep = "_vs_")
    ddsSE$condition <- relevel(ddsSE$condition, control)
  }
  if (type == "manual") {
    if (is.null(contrast)) {
      stop("run test_diff(type = 'manual') with a 'contrast' argument")
    }
    assertthat::assert_that(is.character(contrast))
    if (any(!unlist(strsplit(contrast, "_vs_")) %in% conditions)) {
      stop("run rna.workflow() with valid contrasts in 'contrast'",
           ".\nValid contrasts should contain combinations of: '",
           paste0(conditions, collapse = "', '"),
           "', for example '", paste0(conditions[1],
                                      "_vs_", conditions[2]), "'.", call. = FALSE)
    }
    #ddsSE$condition <- relevel(ddsSE$condition, unlist(strsplit(contrast, "_vs_"))[2])
    cntrst <- contrast
  }
  if (type == "all") {
    cntrst <- apply(utils::combn(conditions, 2), 2, paste,
                    collapse = "_vs_")
    if (!is.null(control)) {
      flip <- grep(paste("^", control, "_vs_", sep = ""),
                   cntrst)
      if (length(flip) >= 1) {
  # rewrite without magrittr placeholder to avoid visible binding NOTE
  tmp_flip <- gsub(paste(control, "_vs_", sep = ""), "", cntrst[flip])
  cntrst[flip] <- paste(tmp_flip, "_vs_", control, sep = "")
      }
    }
  }

  # --- detect replicate structure ---
  n_per_cond <- table(SummarizedExperiment::colData(se_imp)$condition)
  singleton_mode <- !all(n_per_cond >= 2)
  if (singleton_mode && !allow_no_replicates) {
    stop("No biological replicates per condition. Set allow_no_replicates=TRUE or add replicates.", call. = FALSE)
  }


  # Pre-filtering:  keep only rows that have at least 3 samples with a count of 5 or more
  keep <- MatrixGenerics::rowSums(BiocGenerics::counts(ddsSE) >= 5) >= 3
  ddsSE <- ddsSE[keep,]

  # Fit model if replicates exist; otherwise run normalization + transform only
  if (!singleton_mode) {
    # Differential expression analysis; Performs 'replaceOutliers'
    dds <- DESeq2::DESeq(ddsSE, sfType = "ratio", quiet = quiet)
    if(!quiet) message("performing PCA analysis")
    norm.counts <- BiocGenerics::counts(dds, normalized = TRUE)
    if ( trans_method == "rlog"){
      rlog.counts <- tryCatch(DESeq2::rlog(dds, fitType = 'parametric'),
                              error = function(e) { DESeq2::rlog(dds, fitType = 'parametric') })
    } else {
      rlog.counts <- tryCatch(DESeq2::vst(dds, fitType = 'parametric'),
                              error = function(e) { DESeq2::rlog(dds, fitType = 'parametric') })
    }
  } else {
    if(!quiet) message("No replicates detected. Skipping DESeq model fit. Performing size-factor normalization and blind transform.")
    # Estimate size factors only
    dds <- DESeq2::estimateSizeFactors(ddsSE)
    norm.counts <- BiocGenerics::counts(dds, normalized = TRUE)
    # Use blind transformation which does not require replicates
    if ( trans_method == "rlog"){
      rlog.counts <- tryCatch(DESeq2::rlog(dds, fitType = 'parametric', blind = TRUE),
                              error = function(e) { DESeq2::rlog(dds, fitType = 'parametric', blind = TRUE) })
    } else {
      rlog.counts <- tryCatch(DESeq2::vst(dds, fitType = 'parametric', blind = TRUE),
                              error = function(e) { DESeq2::rlog(dds, fitType = 'parametric', blind = TRUE) })
    }
  }
  rna.pca <- prot.pca(SummarizedExperiment::assay(rlog.counts), center=TRUE,scale=TRUE)

  # Create list with test results for defined contrasts
  if(!quiet) message("assembling results for defined contrasts")

  if (!singleton_mode) {
    results <- list()
  for(i in seq_along(cntrst)){
      trmt <- str_split(cntrst[i], "_vs_")[[1]][1]
      ctrl <- str_split(cntrst[i], "_vs_")[[1]][2]
      if(!is.null(altHypothesis)){
        results[[i]] <- DESeq2::results(dds, contrast = c("condition", trmt, ctrl),
                                        independentFiltering = ifelse(!is.null(alpha.independent), TRUE, FALSE),
                                        altHypothesis = altHypothesis, lfcThreshold = lfc,
                                        alpha = alpha.independent,
                                        filterFun = ifelse(pAdjustMethod == "IHW", IHW::ihw, NULL))
      } else {
        results[[i]] <- DESeq2::results(dds, contrast = c("condition", trmt, ctrl),
                                        independentFiltering = ifelse(!is.null(alpha.independent), TRUE, FALSE),
                                        alpha = alpha.independent,
                                        filterFun = ifelse(pAdjustMethod == "IHW", IHW::ihw, NULL))
      }
    }
    names(results) <- cntrst

    # Store condition levels in object
    levels <- levels(ddsSE$condition)

    # --- Revised lfcShrink block ---
    if(lfcShrink == TRUE){
      if(shrink.method == "apeglm"){
        if(!quiet) message("performing lfc shrinkage with method 'apeglm'")
        results_shrink <- list()
        for(i in seq_along(cntrst)) {
          parts <- str_split(cntrst[i], "_vs_")[[1]]
          trmt <- parts[1]
          ctrl <- parts[2]
          dds_temp <- dds
          dds_temp$condition <- relevel(dds_temp$condition, ref = ctrl)
          dds_temp <- DESeq2::nbinomWaldTest(dds_temp, quiet = TRUE)
          coef_name <- paste0("condition", trmt)
          coef_index <- match(coef_name, DESeq2::resultsNames(dds_temp))
          if(is.na(coef_index)) {
            stop("Coefficient ", coef_name, " not found in resultsNames(dds_temp)")
          }
          shrink_res <- DESeq2::lfcShrink(dds_temp, coef = coef_index, type = "apeglm", quiet = TRUE)
          results_shrink[[ cntrst[i] ]] <- shrink_res
        }
      } else {
        if(!quiet) message(paste0("performing lfc shrinkage with method '", shrink.method, "'"))
        results_shrink <- list()
        for(i in seq_along(cntrst)) {
          parts <- str_split(cntrst[i], "_vs_")[[1]]
          trmt <- parts[1]
          ctrl <- parts[2]
          results_shrink[[i]] <- DESeq2::lfcShrink(dds, contrast = c("condition", trmt, ctrl),
                                                   type = shrink.method,
                                                   lfcThreshold = ifelse(shrink.method == "normal", lfc, 0),
                                                   quiet = TRUE)
        }
        names(results_shrink) <- cntrst
      }
      for(i in seq_along(cntrst)){
        results[[ cntrst[i] ]]$lfc.shrink <- as.vector(results_shrink[[ cntrst[i] ]]$log2FoldChange)
        results[[ cntrst[i] ]]$lfcSE.shrink <- as.vector(results_shrink[[ cntrst[i] ]]$lfcSE)
      }
    } # --- End revised lfcShrink block ---

  } else {
    # singleton_mode: perform optional DE using edgeR with fixed dispersion
    if (de_if_no_reps == "edgeR_fixed") {
      if(!quiet) message(sprintf("No replicates: running edgeR GLM with fixed dispersion=%.3f", assumed_dispersion))
      grp <- droplevels(se_imp$condition)
      dge <- edgeR::DGEList(counts = SummarizedExperiment::assay(se_imp), group = grp)
      dge <- edgeR::calcNormFactors(dge, method = "TMM")
      design_edge <- stats::model.matrix(~0 + grp)
      colnames(design_edge) <- levels(grp)
      fit <- edgeR::glmFit(dge, design_edge, dispersion = assumed_dispersion)

      results <- setNames(vector("list", length(cntrst)), cntrst)
      for (i in seq_along(cntrst)) {
        parts <- str_split(cntrst[i], "_vs_")[[1]]
        trmt <- parts[1]
        ctrl <- parts[2]
        # L <- edgeR::makeContrasts(contrasts = paste0(trmt, "-", ctrl), levels = design_edge)
        #lrt <- edgeR::glmLRT(fit, contrast = L[,1])
        contr <- rep(0, ncol(design_edge))
        names(contr) <- colnames(design_edge)
        contr[trmt] <- 1
        contr[ctrl] <- -1
        lrt <- edgeR::glmLRT(fit, contrast = contr)
        tt  <- edgeR::topTags(lrt, n = Inf)$table
        tt  <- tt[match(rownames(se_imp), rownames(tt)), ]
        results[[i]] <- data.frame(
          log2FoldChange = tt$logFC,
          lfcSE = NA_real_,
          pvalue = tt$PValue,
          padj = stats::p.adjust(tt$PValue, method = pAdjustMethod),
          row.names = rownames(se_imp)
        )
      }
    } else {
      if(!quiet) message("No replicates: skipping DE as requested (de_if_no_reps='none').")
      results <- list()
    }
  }

  # add significant column for each contrast, if any
  if (length(results) > 0) {
    for(i in seq_along(results)){
      if ("lfc.shrink" %in% colnames(results[[i]])) {
        results[[i]]$significant <- if(!is.null(altHypothesis)) {
          results[[i]]$padj <= alpha
        } else {
          abs(results[[i]]$lfc.shrink) >= lfc & results[[i]]$padj <= alpha
        }
      } else {
        results[[i]]$significant <- if(!is.null(altHypothesis)) {
          results[[i]]$padj <= alpha
        } else {
          abs(results[[i]]$log2FoldChange) >= lfc & results[[i]]$padj <= alpha
        }
      }
    }
  }

  if(!quiet) message("gathering results into list")
  # store results values for each contrast in dds object
  for(i in seq_along(results)){
    SummarizedExperiment::rowData(dds)[,paste0("lfc.", names(results)[i])] <- results[[i]]$log2FoldChange
    SummarizedExperiment::rowData(dds)[,paste0("lfcSE.", names(results)[i])] <- results[[i]]$lfcSE
    if(lfcShrink == TRUE){
      SummarizedExperiment::rowData(dds)[,paste0("lfc_shrink.", names(results)[i])] <- results[[i]]$lfc.shrink
      SummarizedExperiment::rowData(dds)[,paste0("lfcSE_shrink.", names(results)[i])] <- results[[i]]$lfcSE.shrink
    }
    SummarizedExperiment::rowData(dds)[,paste0("pvalue.", names(results)[i])] <- results[[i]]$pvalue
    SummarizedExperiment::rowData(dds)[,paste0("padj.", names(results)[i])] <- results[[i]]$padj
    SummarizedExperiment::rowData(dds)[,paste0("significant.", names(results)[i])] <- results[[i]]$significant
  }
  # Add 'significant' column if gene was significant in any contrast
  if(length(grep("significant.", colnames(SummarizedExperiment::rowData(dds)))) > 1){
    SummarizedExperiment::rowData(dds)[,"significant"] <- apply(SummarizedExperiment::rowData(dds)[,grep("significant.", colnames(SummarizedExperiment::rowData(dds)))], 1, any)
  }
  else{
    SummarizedExperiment::rowData(dds)[,"significant"] <- SummarizedExperiment::rowData(dds)[,grep("significant.", colnames(SummarizedExperiment::rowData(dds)))]
  }
  if(!quiet) message("Significantly different genes identified:")
  if(!quiet){
  for (i in seq_along(cntrst)) {
      message(paste0("* ",
                     nrow(dds[SummarizedExperiment::rowData(dds)[[paste0("significant.", cntrst[i])]][!is.na(SummarizedExperiment::rowData(dds)[[paste0("significant.", cntrst[i])]])],]),
                     " genes for contrast ", cntrst[i], "\n"))
    }
  }
  # Perform pathway enrichment analysis
  if (pathway_enrichment) {
    if(!quiet) message("performing pathway enrichment analysis")
    res.pathway <- enrich_pathways(dds, contrasts = cntrst, alpha_pathways = alpha_pathways, pathway_kegg=pathway_kegg, kegg_organism, custom_pathways)
  }
  else {
    res.pathway <- NA
  }

  # Write output (robust to missing columns)
  df <- data.frame(SummarizedExperiment::rowData(dds), check.names = FALSE)
  n_df <- nrow(df)

  # Ensure expected metadata columns exist
  if (!"ID"   %in% names(df)) df$ID   <- rownames(df)
  if (!"name" %in% names(df)) df$name <- rownames(df)
  if (!"imputed"  %in% names(df)) df$imputed  <- NA
  if (!"baseMean" %in% names(df)) df$baseMean <- rep(NA_real_, n_df)
  if (!"baseVar"  %in% names(df)) df$baseVar  <- rep(NA_real_, n_df)

  # Collect DE result columns if present
  lfc_cols <- grep("^lfc\\.", names(df))
  end_sig  <- match("significant", names(df))
  if (length(lfc_cols) > 0 && !is.na(end_sig)) {
    res_cols <- df[, lfc_cols[1]:end_sig, drop = FALSE]
  } else {
    res_cols <- data.frame()
  }

  res.table <- cbind(
    data.frame(
      gene_symbol = df$ID,
      gene_name   = df$name,
      imputed     = df$imputed,
      baseMean    = df$baseMean,
      baseVar     = df$baseVar,
      check.names = FALSE
    ),
    res_cols
  )

  if(!quiet) message(paste0("Save results as tab-delimited table to: ", getwd(), "/results.txt"))
  utils::write.table(res.table, file.path(getwd(), "results.txt"), row.names = FALSE, sep = "\t")

  param <- data.frame(type, design, altHypothesis = ifelse(!is.null(altHypothesis), altHypothesis, NA), controlGenes = ifelse(is.null(controlGenes), "none", controlGenes), pAdjustMethod,
                      lfcShrink, shrink.method, alpha, lfc, check.names = FALSE)

  res <- list(data = SummarizedExperiment::rowData(se), se = se, imputed = se_imp,
              dds = dds, pca = rna.pca,
              results = results, param = param)

  #### ----- Sample covariate correlations (gene-level & pathway-level) ----- ####
  if (!is.null(sample_covariates) && (correlate_genes || correlate_pathways)) {
    if (!all(rownames(sample_covariates) %in% colnames(rlog.counts))) {
      warning("Some rownames in 'sample_covariates' do not match sample names; they will be ignored.")
    }
    # Keep only samples present in expression matrix
    common_samples <- intersect(rownames(sample_covariates), colnames(rlog.counts))
    if (length(common_samples) < 3) {
      warning("<3 overlapping samples between expression data and 'sample_covariates'; skipping correlation analysis.")
    } else {
      cov_df <- sample_covariates[common_samples, , drop = FALSE]
      expr_mat <- SummarizedExperiment::assay(rlog.counts)
      expr_mat <- expr_mat[, common_samples, drop = FALSE]

      # Gene-level correlations
      if (correlate_genes) {
        if (!quiet) message("Computing gene-level correlations with sample covariates")
        gene_cor <- tryCatch(
          correlate_matrix_with_covariates(
            mat = expr_mat,
            covariates = cov_df,
            method = correlation_method,
            adjust_method = correlation_p_adjust,
            feature_type = "gene"
          ),
          error = function(e){ warning("Gene correlation failed: ", e$message); NULL }
        )
        res$gene_covariate_correlations <- gene_cor
        if (!is.null(gene_cor) && export) {
          utils::write.table(gene_cor, file.path(getwd(), "gene_covariate_correlations.txt"), sep = "\t", row.names = FALSE)
        }
      }

      # Pathway-level correlations (requires custom_pathways gene sets)
      if (correlate_pathways) {
        if (is.null(custom_pathways)) {
          warning("'correlate_pathways=TRUE' requires 'custom_pathways'; skipping pathway correlation.")
        } else {
          if (!quiet) message("Computing pathway-level correlations with sample covariates")
          # Accept either dataframe with Pathway/Accession or GMT filepath
          pathway_list <- list()
          if (is.data.frame(custom_pathways)) {
            if (!all(c("Pathway","Accession") %in% colnames(custom_pathways))) {
              warning("custom_pathways lacks 'Pathway' and 'Accession' columns; skipping pathway correlation.")
            } else {
              pathway_list <- lapply(seq_len(nrow(custom_pathways)), function(i){
                acc <- custom_pathways$Accession[i]
                # split by common delimiters
                unlist(strsplit(acc, "[,;/] ?| // "))
              })
              names(pathway_list) <- custom_pathways$Pathway
            }
          } else if (is.character(custom_pathways) && length(custom_pathways) == 1 && grepl("\\.gmt$", custom_pathways) && file.exists(custom_pathways)) {
            gmt_df <- clusterProfiler::read.gmt(custom_pathways)
            pathway_list <- split(gmt_df$gene, gmt_df$term)
          } else {
            warning("Unsupported 'custom_pathways' format for pathway correlation; expected data.frame or GMT filepath.")
          }
          if (length(pathway_list) > 0) {
            pw_cor <- tryCatch(
              correlate_pathways(
                expr_mat = expr_mat,
                pathway_list = pathway_list,
                covariates = cov_df,
                min_genes = pathway_min_genes,
                method = correlation_method,
                adjust_method = correlation_p_adjust
              ),
              error = function(e){ warning("Pathway correlation failed: ", e$message); NULL }
            )
            res$pathway_covariate_correlations <- pw_cor
            if (!is.null(pw_cor) && export) {
              utils::write.table(pw_cor, file.path(getwd(), "pathway_covariate_correlations.txt"), sep = "\t", row.names = FALSE)
            }
          }
        }
      }
    }
  }

  if (pathway_enrichment == T && pathway_kegg) {
    res <- c(res, pora_kegg_up = list(res.pathway$ls.pora_kegg_up), pora_kegg_dn = list(res.pathway$ls.pora_kegg_dn))
    if(!quiet) message(paste0("Writing results of KEGG pathway enrichment analysis to: ", getwd(), "/pora_kegg_contrast...txt"))
  for(i in seq_along(res.pathway$ls.pora_kegg_up)){
      if(!is.null(res.pathway$ls.pora_kegg_up[[i]])){
        utils::write.table(res.pathway$ls.pora_kegg_up[[i]]@result, paste(getwd(), paste0("pora_kegg_", names(res.pathway$ls.pora_kegg_up)[i], "_up.txt"),
                                                                          sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
  for(i in seq_along(res.pathway$ls.pora_kegg_dn)){
      if(!is.null(res.pathway$ls.pora_kegg_dn[[i]])){
        utils::write.table(res.pathway$ls.pora_kegg_dn[[i]]@result, paste(getwd(), paste0("pora_kegg_", names(res.pathway$ls.pora_kegg_dn)[i], "_down.txt"),
                                                                          sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
  }
  if (pathway_enrichment == T && !is.null(custom_pathways)) {
    res <- c(res, pora_custom_up = list(res.pathway$ls.pora_custom_up), pora_custom_dn = list(res.pathway$ls.pora_custom_dn))
    if(!quiet) message(paste0("Writing results of custom pathway enrichment analysis to: ", getwd(), "/pora_custom_contrast...txt"))
  for(i in seq_along(res.pathway$ls.pora_custom_up)){
      if(!is.null(res.pathway$ls.pora_custom_up[[i]])){
        utils::write.table(res.pathway$ls.pora_custom_up[[i]]@result, paste(getwd(), paste0("pora_custom_", names(res.pathway$ls.pora_custom_up)[i], "_up.txt"),
                                                                            sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
  for(i in seq_along(res.pathway$ls.pora_custom_dn)){
      if(!is.null(res.pathway$ls.pora_custom_dn[[i]])){
        utils::write.table(res.pathway$ls.pora_custom_dn[[i]]@result, paste(getwd(), paste0("pora_custom_", names(res.pathway$ls.pora_custom_dn)[i], "_down.txt"),
                                                                            sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
  }

  # GSEA for each contrast (after ORA output)
  if (gsea) {
    # Decide which GMT to use: either gsea_gmt, or fallback to custom_pathways if it's a .gmt file
    gmt_file_to_use <- NULL
    if (!is.null(gsea_gmt)) {
      gmt_file_to_use <- gsea_gmt
    } else if (is.character(custom_pathways) && grepl("\\.gmt$", custom_pathways)) {
      gmt_file_to_use <- custom_pathways
    }
    if (is.null(gmt_file_to_use)) {
      stop("GSEA requested but no GMT file provided in 'gsea_gmt', and 'custom_pathways' is not a GMT file.")
    }

    if (!quiet) message("performing GSEA for each contrast")
    ls_gsea_up <- vector("list", length(cntrst))
    names(ls_gsea_up) <- cntrst
    #browser()

    for (i in seq_along(cntrst)) {
      ## 1) extract ranked vector from DESeq2 results
      res_df <- as.data.frame(results[[cntrst[i]]])
      gene_list <- res_df$log2FoldChange
      names <- rownames(res_df)
      IDs <- SummarizedExperiment::rowData(dds)[match(names, SummarizedExperiment::rowData(dds)$name), "ID"]
      names(gene_list) <- IDs
      gene_list <- sort(gene_list, decreasing = TRUE)

      set.seed(1234)
      ## 2) run gsea_custom()
      gsea_res <- gsea_custom(
        geneList      = gene_list,
        gmt_file      = gmt_file_to_use,
        pAdjustMethod = gsea_pAdjustMethod,
        pvalueCutoff  = gsea_pvalueCutoff,
        minGSSize     = gsea_minGSSize,
        maxGSSize     = gsea_maxGSSize
      )

      ls_gsea_up[[cntrst[i]]] <- gsea_res

      ## 3) write to disk if non‐empty
      if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
        utils::write.table(
          as.data.frame(gsea_res),
          paste0(getwd(), "/gsea_", cntrst[i], ".txt"),
          row.names = FALSE, sep = "\t"
        )
      }
    }

    ## 4) attach to return list
    res$gsea_results <- ls_gsea_up
  }

  # Create Plots
  if(export == TRUE || plot == TRUE){
    if (export == TRUE){
      dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
      if(!quiet) message(paste0("Rendering and exporting figures as PDF and PNG files to: ", getwd(), "/Plots"))
    }
    suppress_ggrepel <- function(w) {
      if (any(grepl("ggrepel", w)))
        invokeRestart("muffleWarning")
    }
    #### Box plot of Cook's distance for each sample (outlier detection)
    if(plot==TRUE){
      if (!singleton_mode) {
        par(mar=c(8,5,2,2))
  boxplot(log10(SummarizedExperiment::assays(dds)[["cooks"]]), range=0, las=2)
      }
    }
    if (export == TRUE){
      grDevices::pdf(paste0("Plots/BoxPlot_CooksDistance", ".pdf"),
                     width = 6*length(dds$condition)/20, height = 6)
      if (!singleton_mode) {
        par(mar=c(8,5,2,2))
  boxplot(log10(SummarizedExperiment::assays(dds)[["cooks"]]), range=0, las=2)
      }
      grDevices::dev.off()

      grDevices::png(paste0("Plots/BoxPlot_CooksDistance", ".png"),
                     width = 6*length(dds$condition)/20, height = 6, units = 'in', res = 300)
      if (!singleton_mode) {
        par(mar=c(8,5,2,2))
  boxplot(log10(SummarizedExperiment::assays(dds)[["cooks"]]), range=0, las=2)
      }
      grDevices::dev.off()
    }
    try(suppressMessages(
      prot.plot_missval(se, plot = plot, export = export)
    ))

    se_log2 <- se
  SummarizedExperiment::assay(se_log2) <- log2(SummarizedExperiment::assay(se_log2))
    try(suppressMessages(
      prot.plot_detect(se_log2, basesize = 10, plot = plot, export = export)
    ))
    suppressMessages(
      rna.plot_imputation(log2(SummarizedExperiment::assay(se)+1),
                          "Imputed" = log2(SummarizedExperiment::assay(se_imp)+1),
                          "Normalized" = log2(norm.counts+1), colData = SummarizedExperiment::colData(se),
                          plot = plot, export = export, basesize = 12)
    )
    suppressMessages(
      suppressWarnings(
        prot.plot_screeplot(rna.pca, axisLabSize = 18, titleLabSize = 22, plot = plot, export = export)
      ) )
    withCallingHandlers(suppressMessages(
      prot.plot_loadings(rna.pca,labSize = 3, plot = plot, export = export)
    ) , warning = suppress_ggrepel)
    suppressMessages(
      rna.plot_pca(dds,x = 1, y = 2, point_size = 4, basesize = 14, title = "PC Scores - PC2 vs. PC1",
                   plot = plot,export = export)
    )
    suppressMessages(
      rna.plot_pca(dds,x = 1, y = 3, point_size = 4, basesize = 14, title = "PC Scores - PC3 vs. PC1",
                   plot = plot,export = export)
    )
    suppressMessages(
      rna.plot_heatmap(dds, type = "centered", kmeans = heatmap.kmeans, show_all = heatmap.show_all, contrast = cntrst,
                       k = k, col_limit = heatmap.col_limit,show_row_names = heatmap.show_row_names,
                       row_font_size = heatmap.row_font_size, indicate = c("condition"),
                       plot = plot, export = export)
    )
    suppressMessages(
      rna.plot_heatmap(dds, type = "contrast", contrast = cntrst,  kmeans = heatmap.kmeans, k = k, col_limit = heatmap.col_limit,
                       show_row_names = heatmap.show_row_names, row_font_size = heatmap.row_font_size,
                       plot = plot, export = export)
    )

    # Define variable 'volcano_plotting' if plot is TRUE and plot_volcano is TRUE
    if (plot == TRUE && plot_volcano == TRUE) {
      volcano_plotting <- TRUE
    } else {
      volcano_plotting <- FALSE
    }
    if (export == TRUE && plot_volcano == TRUE) {
      volcano_exporting <- TRUE
    } else {
      volcano_exporting <- FALSE
    }
  if (volcano_plotting || volcano_exporting) {
      # for (i in 1:length(contrasts)){
      #   suppressMessages(
      #     suppressWarnings(
      #       prot.plot_volcano(prot_dep, contrast = contrasts[i],
      #                         add_names = volcano.add_names, label_size = volcano.label_size, adjusted =  volcano.adjusted,
      #                         plot = volcano_plotting, export = volcano_exporting, lfc = lfc, alpha = alpha)
      #     ) )
      # }
  for (i in seq_along(cntrst)){
        suppressMessages(
          suppressWarnings(
            volcano.tmp <- rna.plot_volcano(dds, contrast = cntrst[i],
                                            add_names = volcano.add_names, label_size = volcano.label_size, adjusted =  volcano.adjusted,
                                            plot = volcano_plotting, export = volcano_exporting, lfc = lfc, alpha = alpha)
          ) )
      }
    }
    # for (i in 1:length(cntrst)){
    #   suppressMessages(
    #     suppressWarnings(
    #       volcano.tmp <- rna.plot_volcano(dds, contrast = cntrst[i],
    #                                       add_names = volcano.add_names, label_size = volcano.label_size, adjusted =  volcano.adjusted,
    #                                       plot = plot, export = export, lfc = lfc, alpha = alpha)
    #     ) )
    # }

    if (pathway_kegg) {
  for (i in seq_along(cntrst)) {
        if(!(nrow(as.data.frame(res.pathway$ls.pora_kegg_up[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(res.pathway$ls.pora_kegg_up[[i]], title = paste0("Upregulated pathways", " - KEGG"),
                                   subtitle =  str_replace(cntrst[i], "_vs_", " vs. "), plot = plot, export = export, kegg = TRUE)
            ) )
        }
        if(!(nrow(as.data.frame(res.pathway$ls.pora_kegg_dn[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(res.pathway$ls.pora_kegg_dn[[i]], title = paste0("Downregulated pathways", " - KEGG"),
                                   subtitle = str_replace(cntrst[i], "_vs_", " vs. "), plot = plot, export = export, kegg = TRUE)
            ) )
        }
      }
    }
    if(!is.null(custom_pathways)){
  for (i in seq_along(cntrst)) {
        if(!(nrow(as.data.frame(res.pathway$ls.pora_custom_up[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(res.pathway$ls.pora_custom_up[[i]], title = "Upregulated pathways",
                                   subtitle =  str_replace(cntrst[i], "_vs_", " vs. "), plot = plot, export = export, kegg = FALSE)
            ) )
        }
        if(!(nrow(as.data.frame(res.pathway$ls.pora_custom_dn[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(res.pathway$ls.pora_custom_dn[[i]], title = "Downregulated pathways",
                                   subtitle = str_replace(cntrst[i], "_vs_", " vs. "), plot = plot, export = export, kegg = FALSE)
            ) )
        }
      }
    } # if(!is.null(custom_pathways))
  } #if(export == TRUE || plot == TRUE)

  if(report == TRUE){
    rna.report(res,
               volcano.adjusted = volcano.adjusted,
               volcano.add_names = volcano.add_names,
               pathway_enrichment = pathway_enrichment,
               heatmap.show_all = heatmap.show_all,
               heatmap.kmeans = heatmap.kmeans,
               heatmap.row_font_size = heatmap.row_font_size,
               export = export,
               k = k,
               report.dir = report.dir)
  }
  message("Done!")
  return(res)
}
