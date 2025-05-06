#' @title RNA sequencing workflow
#'
#' @description This function performs a complete RNA sequencing workflow, including imputation of missing values, normalization,
#' principal component analysis, differential expression analysis, and pathway analysis. The function also provides
#' several options for plotting, exporting plots, and creating a report.
#'
#' @param se A SummarizedExperiment object, generated with read_prot().
#' @param imp_fun (Character string)  Function used for data imputation. "SampMin", "man", "bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb", "min", "zero", "mixed", or "nbavg". See (\code{\link{rna.impute}}) for details.
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
#' @param export Shall plots be exported as PDF and PNG files?
#' @param report Shall a report (HTML and PDF) be created?
#' @param report.dir Folder name for created report (if report = TRUE)
#' @param pathway_enrichment Perform pathway over-representation analysis for each tested contrast
#' @param pathway_kegg Perform pathway over-representation analysis with gene sets in the KEGG database
#' @param kegg_organism Name of the organism in the KEGG database (if 'pathway_kegg = TRUE')
#' @param custom_pathways Dataframe providing custom pathway annotations.
#' @param quiet Suppress messages and warnings.
#'
#' @return The function returns a SummarizedExperiment object with added columns for log2 fold change, p-values and adjusted p-values for each comparison.
#' It also includes a column for significant genes for each comparison and a column for significant genes overall.
#' Additionally, the function generates various plots and a report (if specified).
#' @export
rna.workflow <- function(se, # SummarizedExperiment, generated with read_prot().
                         imp_fun = c("zero", "man", "bpca", "knn", "QRILC", "MLE", "MinDet",
                                     "MinProb", "min", "zero", "mixed", "nbavg", "SampMin"),
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
                         export = FALSE,
                         report = TRUE,
                         report.dir = NULL,
                         pathway_enrichment = FALSE,
                         pathway_kegg = FALSE,
                         kegg_organism = NULL,
                         custom_pathways = NULL,
                         quiet = FALSE
)
{
  # Show error if inputs are not the required classes
  imp_fun <- match.arg(imp_fun)
  shrink.method <- match.arg(shrink.method)
  pAdjustMethod <- match.arg(pAdjustMethod)
  if(is.integer(alpha)) alpha <- as.numeric(alpha)
  if(is.integer(lfc)) lfc <- as.numeric(lfc)
  assertthat::assert_that(is.character(imp_fun),
                          is.character(type),
                          is.numeric(alpha),
                          length(alpha) == 1,
                          is.numeric(lfc),
                          length(lfc) == 1
  )

  # Show error if inputs are not valid
  if (!type %in% c("all", "control", "manual")) {
    stop("run workflow_proteomics() with a valid type",
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
  message("Number of genes after pre-filtering: ", nrow(assay(se)))
  # Impute missing values
  if (imp_fun == "MinProb"){
    se_imp <- rna.impute(se, fun = imp_fun, q = q)
  } else if (imp_fun == "knn"){
    se_imp <- rna.impute(se, fun = imp_fun, rowmax = knn.rowmax)
  } else {
    se_imp <- rna.impute(se, fun = imp_fun)
  }
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
        cntrst[flip] <- cntrst[flip] %>% gsub(paste(control,
                                                    "_vs_", sep = ""), "", .) %>%
          paste("_vs_", control, sep = "")
      }
    }
  }


  # Pre-filtering:  keep only rows that have at least 3 samples with a count of 5 or more
  keep <- MatrixGenerics::rowSums(BiocGenerics::counts(ddsSE) >= 5) >= 3
  ddsSE <- ddsSE[keep,]

  # Differential expression analysis; Performs 'replaceOutliers' (replaces outlier counts flagged by extreme Cook's distances)
  dds <- DESeq2::DESeq(ddsSE, sfType = "ratio", quiet = quiet)

  # Perform PCA Analysis
  if(!quiet) message("performing PCA analysis")
  norm.counts <- BiocGenerics::counts(dds, normalized = TRUE)
  rlog.counts <- tryCatch(DESeq2::rlog(dds, fitType = 'mean'), error = function(e) { rlog(dds, fitType = 'mean') })
  rna.pca <- prot.pca(SummarizedExperiment::assay(rlog.counts))

  # Create list with test results for defined contrasts
  if(!quiet) message("assembling results for defined contrasts")
  results <- list()
  for(i in 1:length(cntrst)){
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

  # # Creation of shrunken results list after re-running nbinomWaldTest with re-ordered conditions
  # if(lfcShrink == TRUE){
  #   if(!quiet) message(paste0("performing lfc shrinkage with method '", shrink.method, "'"))
  #   results_shrink <- list()
  #   if(shrink.method == "apeglm"){
  #     cntrst_ndx <- lapply(2:length(DESeq2::resultsNames(dds)), function(x) match(gsub("condition_", "", DESeq2::resultsNames(dds)[x]), cntrst))
  #     cntrst_ndx <- cntrst_ndx[!sapply(cntrst_ndx,is.na)]
  #     names(cntrst_ndx) <- cntrst[unlist(cntrst_ndx)]
  #     if(length(cntrst_ndx) > 0){
  #       for(i in 1:length(cntrst_ndx)){
  #         results_shrink[[i]] <- DESeq2::lfcShrink(dds, coef = match(names(cntrst_ndx)[i], gsub("condition_", "", DESeq2::resultsNames(dds))),
  #                                                  type = "apeglm", returnList = F, quiet = T, )
  #       }
  #     }
  #     names(results_shrink) <- names(cntrst_ndx)
  #     if (type == "all") {
  #       for(i in 2:length(levels)){
  #         dds$condition <- relevel(dds$condition, levels[i])
  #         dds <- DESeq2::nbinomWaldTest(object = dds, quiet = T)
  #         cntrst_ndx <- lapply(2:length(DESeq2::resultsNames(dds)), function(x) match(gsub("condition_", "", DESeq2::resultsNames(dds)[x]), cntrst))
  #         cntrst_ndx <- cntrst_ndx[!sapply(cntrst_ndx,is.na)]
  #         names(cntrst_ndx) <- cntrst[unlist(cntrst_ndx)]
  #         if(length(cntrst_ndx) > 0){
  #           for(j in 1:length(cntrst_ndx)){
  #             results_shrink[[length(results_shrink)+1]] <- DESeq2::lfcShrink(dds, coef = match(names(cntrst_ndx)[j], gsub("condition_", "", DESeq2::resultsNames(dds))),
  #                                                                             type = "apeglm", returnList = F, quiet = T)
  #             names(results_shrink)[length(results_shrink)] <- names(cntrst_ndx)[j]
  #           }
  #         }
  #       }
  #     }
  #   }
  #   if(shrink.method == "ashr" || shrink.method == "normal"){
  #     for (i in 1:length(results)) {
  #       trmt <- str_split(cntrst[i], "_vs_")[[1]][1]
  #       ctrl <- str_split(cntrst[i], "_vs_")[[1]][2]
  #       results_shrink[[i]] <- DESeq2::lfcShrink(dds, contrast = c("condition", trmt, ctrl), type = shrink.method,
  #                                                lfcThreshold = ifelse(shrink.method == "normal", lfc, 0), quiet = T)
  #     }
  #     names(results_shrink) <- cntrst
  #   }
  #   # add shrunken lfc values to results object
  #   for(i in 1:length(cntrst)){
  #     results[[match(cntrst[i], names(results))]]$lfc.shrink <- as.vector(results_shrink[[match(cntrst[i], names(results_shrink))]]$log2FoldChange)
  #     results[[match(cntrst[i], names(results))]]$lfcSE.shrink <- as.vector(results_shrink[[match(cntrst[i], names(results_shrink))]]$lfcSE)
  #   }
  # } # if(lfcShrink == TRUE)

  # --- Revised lfcShrink block ---
  if(lfcShrink == TRUE){
    if(shrink.method == "apeglm"){
      if(!quiet) message("performing lfc shrinkage with method 'apeglm'")
      # For apeglm, we need to use a coefficient-based approach.
      # Loop over each contrast, relevel dds so that the control becomes the reference,
      # re-run the Wald test, and then call lfcShrink using the appropriate coefficient.
      results_shrink <- list()
      for(i in seq_along(cntrst)) {
        parts <- str_split(cntrst[i], "_vs_")[[1]]
        trmt <- parts[1]
        ctrl <- parts[2]
        # Make a temporary copy of dds
        dds_temp <- dds
        # Relevel so that ctrl is the reference level
        dds_temp$condition <- relevel(dds_temp$condition, ref = ctrl)
        # Re-run the Wald test after releveling
        dds_temp <- DESeq2::nbinomWaldTest(dds_temp, quiet = TRUE)
        # In an intercept model, the log fold change is given by the coefficient for the treatment level,
        # which will have the name "condition<trmt>"
        coef_name <- paste0("condition", trmt)
        coef_index <- match(coef_name, DESeq2::resultsNames(dds_temp))
        if(is.na(coef_index)) {
          stop("Coefficient ", coef_name, " not found in resultsNames(dds_temp)")
        }
        shrink_res <- DESeq2::lfcShrink(dds_temp, coef = coef_index, type = "apeglm", quiet = TRUE)
        results_shrink[[ cntrst[i] ]] <- shrink_res
      }
    } else {
      # For "ashr" or "normal", we can use the contrast argument directly.
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
    # Add shrunken lfc values to results object
    for(i in seq_along(cntrst)){
      results[[ cntrst[i] ]]$lfc.shrink <- as.vector(results_shrink[[ cntrst[i] ]]$log2FoldChange)
      results[[ cntrst[i] ]]$lfcSE.shrink <- as.vector(results_shrink[[ cntrst[i] ]]$lfcSE)
    }
  } # --- End revised lfcShrink block ---

  # add significant column for each contrast
  for(i in 1:length(cntrst)){
    if(lfcShrink == TRUE){
      if(!is.null(altHypothesis)){
        results[[match(cntrst[i], names(results))]]$significant <-
          results[[match(cntrst[i], names(results))]]$padj <= alpha
      } else {
        results[[match(cntrst[i], names(results))]]$significant <-
          abs(results[[match(cntrst[i], names(results))]]$lfc.shrink) >= lfc &
          results[[match(cntrst[i], names(results))]]$padj <= alpha
      }
    }
    else{
      if(!is.null(altHypothesis)){
        results[[match(cntrst[i], names(results))]]$significant <-
          results[[match(cntrst[i], names(results))]]$padj <= alpha
      } else {
        results[[match(cntrst[i], names(results))]]$significant <-
          abs(results[[match(cntrst[i], names(results))]]$log2FoldChange) >= lfc &
          results[[match(cntrst[i], names(results))]]$padj <= alpha
      }
    }
  }

  if(!quiet) message("gathering results into list")
  # store results values for each contrast in dds object
  for(i in 1:length(results)){
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
    for (i in 1:length(cntrst)) {
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

  # Write output
  df <- data.frame(SummarizedExperiment::rowData(dds))
  res.table <- data.frame("gene_symbol" = df$ID,
                          "gene_name" = df$name,
                          "imputed" = df$imputed,
                          "baseMean" = df$baseMean,
                          "baseVar" = df$baseVar,
                          df[,grep("lfc\\.", colnames(df))[1]:match("significant", colnames(df))]
  )
  if(!quiet) message(paste0("Save results as tab-delimited table to: ", getwd(), "/results.txt"))
  utils::write.table(res.table, paste(getwd(), "results.txt",
                                      sep = "/"), row.names = FALSE, sep = "\t")

  param <- data.frame(type, design, altHypothesis = ifelse(!is.null(altHypothesis), altHypothesis, NA), controlGenes = ifelse(is.null(controlGenes), "none", controlGenes), pAdjustMethod,
                      lfcShrink, shrink.method, alpha, lfc, check.names = FALSE)

  res <- list(data = SummarizedExperiment::rowData(se), se = se, imputed = se_imp,
              dds = dds, pca = rna.pca,
              results = results, param = param)

  if (pathway_enrichment == T && pathway_kegg) {
    res <- c(res, pora_kegg_up = list(res.pathway$ls.pora_kegg_up), pora_kegg_dn = list(res.pathway$ls.pora_kegg_dn))
    if(!quiet) message(paste0("Writing results of KEGG pathway enrichment analysis to: ", getwd(), "/pora_kegg_contrast...txt"))
    for(i in 1:length(res.pathway$ls.pora_kegg_up)){
      if(!is.null(res.pathway$ls.pora_kegg_up[[i]])){
        utils::write.table(res.pathway$ls.pora_kegg_up[[i]]@result, paste(getwd(), paste0("pora_kegg_", names(res.pathway$ls.pora_kegg_up)[i], "_up.txt"),
                                                                          sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
    for(i in 1:length(res.pathway$ls.pora_kegg_dn)){
      if(!is.null(res.pathway$ls.pora_kegg_dn[[i]])){
        utils::write.table(res.pathway$ls.pora_kegg_dn[[i]]@result, paste(getwd(), paste0("pora_kegg_", names(res.pathway$ls.pora_kegg_dn)[i], "_down.txt"),
                                                                          sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
  }
  if (pathway_enrichment == T && !is.null(custom_pathways)) {
    res <- c(res, pora_custom_up = list(res.pathway$ls.pora_custom_up), pora_custom_dn = list(res.pathway$ls.pora_custom_dn))
    if(!quiet) message(paste0("Writing results of custom pathway enrichment analysis to: ", getwd(), "/pora_custom_contrast...txt"))
    for(i in 1:length(res.pathway$ls.pora_custom_up)){
      if(!is.null(res.pathway$ls.pora_custom_up[[i]])){
        utils::write.table(res.pathway$ls.pora_custom_up[[i]]@result, paste(getwd(), paste0("pora_custom_", names(res.pathway$ls.pora_custom_up)[i], "_up.txt"),
                                                                            sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
    for(i in 1:length(res.pathway$ls.pora_custom_dn)){
      if(!is.null(res.pathway$ls.pora_custom_dn[[i]])){
        utils::write.table(res.pathway$ls.pora_custom_dn[[i]]@result, paste(getwd(), paste0("pora_custom_", names(res.pathway$ls.pora_custom_dn)[i], "_down.txt"),
                                                                            sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
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
      par(mar=c(8,5,2,2))
      boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
    }
    if (export == TRUE){
      grDevices::pdf(paste0("Plots/BoxPlot_CooksDistance", ".pdf"),
                     width = 6*length(dds$condition)/20, height = 6)
      par(mar=c(8,5,2,2))
      boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
      grDevices::dev.off()

      grDevices::png(paste0("Plots/BoxPlot_CooksDistance", ".png"),
                     width = 6*length(dds$condition)/20, height = 6, units = 'in', res = 300)
      par(mar=c(8,5,2,2))
      boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
      grDevices::dev.off()
    }
    try(suppressMessages(
      prot.plot_missval(se, plot = plot, export = export)
    ))

    se_log2 <- se
    assay(se_log2) <- log2(assay(se_log2))
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
    for (i in 1:length(cntrst)){
      suppressMessages(
        suppressWarnings(
          volcano.tmp <- rna.plot_volcano(dds, contrast = cntrst[i],
                                          add_names = volcano.add_names, label_size = volcano.label_size, adjusted =  volcano.adjusted,
                                          plot = plot, export = export, lfc = lfc, alpha = alpha)
        ) )
    }

    if (pathway_kegg) {
      for (i in 1:length(cntrst)) {
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
      for (i in 1:length(cntrst)) {
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
