#' @title Run a complete proteomics analysis workflow.
#'
#' @description \code{prot.workflow} performs variance stabilization normalization (\code{\link{prot.normalize_vsn}}), missing value imputation (\code{\link{prot.impute}}), principal component analysis (\code{\link{prot.pca}}), differential enrichment test (\code{\link{prot.test_diff}}), and pathway enrichment analysis. If desired, standardized plots and a report are generated and exported as separate files.
#'
#' @param se \code{SummarizedExperiment} object, proteomics data parsed with \code{\link{prot.read_data}}.
#' @param normalize (Logical) Should the data be normalized via variance stabilization normalization?
#' @param imp_fun (Character string)  Function used for data imputation. "SampMin", "man", "bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb", "min", "zero", "mixed", or "nbavg". See (\code{\link{prot.impute}}) for details.
#' @param q (Numeric) q value for imputing missing values with method \code{imp_fun = 'MinProb'}.
#' @param knn.rowmax (Numeric) The maximum percent missing data allowed in any row for \code{imp_fun = 'knn'}. Default: 0.5.
#' @param type (Character string) Type of differential analysis to perform. "all" (contrast each condition with every other condition), "control" (contrast each condition to a defined control condition), "manual" (manually define selected conditions).
#' @param control (Character string) The name of the control condition if \code{type = control}.
#' @param contrast (Character string or vector of strings) Define the contrasts to be tested if \code{type = manual} in the form: \code{"ConditionA_vs_ConditionB"}, or \code{c("ConditionA_vs_ConditionC", "ConditionB_vs_ConditionC")}.
#' @param alpha (Numeric) Significance threshold for adjusted p values.
#' @param alpha_pathways (Numeric) Significance threshold for adjusted p values in pathway enrichment analysis.
#' @param lfc (Numeric) Relevance threshold for log2(fold change) values. Only proteins with a |log2(fold change)| value above \code{lfc} for a given contrast are considered "significant" (if they additionally fullfil the \code{alpha} criterion).
#' @param heatmap.show_all (Logical) Shall all samples be displayed in the heat map (\code{TRUE}) or only the samples contained in the defined \code{contrast} (\code{FALSE})?
#' @param heatmap.kmeans (Logical) Shall the proteins be clustered in the heat map with the k-nearest neighbour method (\code{TRUE}) or not \code{FALSE})?
#' @param k (Integer) Number of protein clusters in heat map if \code{heatmap.kmeans = TRUE}.
#' @param heatmap.col_limit (Integer) Define the outer breaks in the heat map legend. Example: if \code{heatmap.col_limit = 3}, the color scale will span from -3 to 3. Alls values below -3 will have the same color as -3, and all values above 3 will have the same color as 3.
#' @param heatmap.show_row_names (Logical) Show protein names in heat map (\code{TRUE}) or not \code{FALSE}).
#' @param heatmap.row_font_size (Numeric) Font size of protein names in heat maps if \code{heatmap.show_row_names = TRUE}.
#' @param volcano.add_names (Logical) Show protein names in volcano plots (\code{TRUE}) or not \code{FALSE}).
#' @param volcano.label_size (Numeric) Font size of protein names in volcano plots if \code{volcano.add_names = TRUE}.
#' @param volcano.adjusted (Logical) Shall adjusted p values be shown on the y axis of volcano plots (\code{TRUE}) or raw p values (\code{FALSE})?.
#' @param plot (Logical) Show the generated plots in the \code{Plots} pane of RStudio (\code{TRUE}) or not \code{FALSE}).
#' @param plot_volcano (Logical) Show the volcano plot in the \code{Plots} pane of RStudio (\code{TRUE}) or not \code{FALSE}).
#' @param export (Logical) Exported the generated plots as PNG and PDF files (\code{TRUE}) or not \code{FALSE}).
#' @param report (Logical) Render and export a report in PDF and HTML format that summarizes the results (\code{TRUE}) or not \code{FALSE}).
#' @param report.dir (Character string) Provide the name of or path to a folder into which the report will be saved.
#' @param pathway_enrichment (Logical) Perform pathway over-representation analysis for each tested contrast (\code{TRUE}) or not \code{FALSE}).
#' @param pathway_kegg (Logical) Perform pathway over-representation analysis with gene sets in the KEGG database (\code{TRUE}) or not \code{FALSE}).
#' @param kegg_organism (Character string) Identifier of the organism in the KEGG database (if \code{pathway_kegg = TRUE})
#' @param custom_pathways (a R dataframe object) Data frame providing custom pathway annotations. The table must contain a "Pathway" column listing identified pathway in the studies organism, and an "Accession" column listing the proteins (or genes) each pathway is composed of. The **Accession** entries must match with protein **IDs**.
#' @param out.dir (Character string) absolute path to the location where result TXT files should be exported to.
#' @return A list containing `SummarizedExperiment` object for every computation step of the workflow, a \code{pca} object, and lists of up- or down regulated pathways for each tested contrast and method (KEGG and/or custom).
#' @export
#' @importFrom assertthat assert_that
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom dplyr select filter
prot.workflow <- function(se, # SummarizedExperiment, generated with read_prot().
                          normalize = TRUE,
                          imp_fun = c("SampMin", "man", "bpca", "knn", "QRILC", "MLE", "MinDet", # Method for imputing of missing values
                                      "MinProb", "min", "zero", "mixed", "nbavg"),
                          q = 0.01, # q value for imputing missing values with method "fun = 'MinProb'".
                          knn.rowmax = 0.5, # The maximum percent missing data allowed in any row (default 50%).
                          # For any rows with more than rowmax% missing are imputed using the overall mean per sample.
                          type = c("all", "control", "manual"), # Type of differential analysis to perform.
                          control = NULL, # Control condition; required if type = "control".
                          contrast = NULL, # Defined test for differential analysis "A_vs_B"; required if type = "manual".
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
                          plot_volcano,
                          export = FALSE, # Shall plots be exported as PDF and PNG files?
                          report = TRUE, # Shall a report (HTML and PDF) be created?
                          report.dir = NULL, # Folder name for created report (if report = TRUE)
                          pathway_enrichment = FALSE, # Perform pathway over-representation analysis for each tested contrast
                          pathway_kegg = FALSE, # Perform pathway over-representation analysis with gene sets in the KEGG database
                          kegg_organism = NULL, # Name of the organism in the KEGG database (if 'pathway_kegg = TRUE')
                          custom_pathways = NULL, # Dataframe providing custom pathway annotations
                          out.dir = NULL)
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

  imp_fun <- match.arg(imp_fun)
  out_dir <- out.dir

  # If out.dir is defined, create the directory and set it as the working directory
  old_wd = getwd()
  if (!is.null(out.dir)) {
    dir.create(out.dir, showWarnings = F)
    setwd(out.dir)
    message("Running proteomics workflow in the directory:", as.character(getwd()))
  }

  # if contrast is defined, set type to "manual"
  if (!is.null(contrast)) {
    type <- "manual"
  }

  # Show error if inputs are not valid
  if (!type %in% c("all", "control", "manual")) {
    stop("run workflow_proteomics() with a valid type",
         "\nValid types are: 'all', 'control' and 'manual'.",
         call. = FALSE)
  }
  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
  }
  # Variance stabilization
  if(normalize == TRUE){
    prot_norm <- suppressMessages(prot.normalize_vsn(se, plot = plot, export = export))
  } else {
    prot_norm <- se
  }
  # Impute missing values
  if (imp_fun == "MinProb"){
    prot_imp <- prot.impute(prot_norm, fun = imp_fun, q = q)
  } else if (imp_fun == "knn"){
    prot_imp <- prot.impute(prot_norm, fun = imp_fun, rowmax = knn.rowmax)
  } else {
    prot_imp <- prot.impute(prot_norm, fun = imp_fun)
  }
  # Perform PCA Analysis
  prot_pca <- prot.pca(SummarizedExperiment::assay(prot_imp), in_workflow = TRUE)
  # Test for differential expression by empirical Bayes moderation
  # of a linear model and defined contrasts
  prot_diff <- prot.test_diff(prot_imp, type = type, control = control, test = contrast)
  # Denote significantly differentially expressed proteins
  prot_dep <- prot.add_rejections(prot_diff, alpha = alpha, lfc = lfc)
  contrasts <- SummarizedExperiment::rowData(prot_dep) %>%
    data.frame(check.names = FALSE) %>%
    select(ends_with("_diff")) %>%
    colnames() %>% str_replace_all("_diff", "")
  # Generate a results table
  results <- prot.get_results(prot_dep)
  n_significant <- results %>% dplyr::filter(significant) %>% nrow()

  message(paste0(n_significant,
                 " proteins were found to be differentially expressed with ",
                 expression(alpha),
                 " = ", alpha,
                 " and |log2(fold change)| > ",
                 lfc, "."))
  # Perform pathway enrichment analysis
  if (pathway_enrichment) {
    res.pathway <- enrich_pathways(prot_dep, contrasts, alpha_pathways = alpha_pathways,
                                   pathway_kegg = pathway_kegg, kegg_organism = kegg_organism,
                                   custom_pathways = custom_pathways)
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
    # Define variable 'volcano_plotting' if plot is TRUE and plot_volcano is TRUE
    if (plot == TRUE && plot_volcano == TRUE) {
      volcano_plotting <- TRUE
    } else {
      volcano_plotting <- FALSE
    }
    for (i in 1:length(contrasts)){
      suppressMessages(
        suppressWarnings(
          prot.plot_volcano(prot_dep, contrast = contrasts[i],
                            add_names = volcano.add_names, label_size = volcano.label_size, adjusted =  volcano.adjusted,
                            plot = volcano_plotting, export = export, lfc = lfc, alpha = alpha)
        ) )
    }
    if (pathway_kegg) {
      for (i in 1:length(contrasts)) {
        if(!(nrow(as.data.frame(res.pathway$ls.pora_kegg_up[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(res.pathway$ls.pora_kegg_up[[i]], title = paste0("Upregulated pathways", " - KEGG"),
                                   subtitle =  str_replace(contrasts[i], "_vs_", " vs. "), plot = plot, export = export, kegg = TRUE)
            ) )
        }
        if(!(nrow(as.data.frame(res.pathway$ls.pora_kegg_dn[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(res.pathway$ls.pora_kegg_dn[[i]], title = paste0("Downregulated pathways", " - KEGG"),
                                   subtitle = str_replace(contrasts[i], "_vs_", " vs. "), plot = plot, export = export, kegg = TRUE)
            ) )
        }
      }
    }
    if(!is.null(custom_pathways)){
      for (i in 1:length(contrasts)) {
        if(!(nrow(as.data.frame(res.pathway$ls.pora_custom_up[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(res.pathway$ls.pora_custom_up[[i]], title = "Upregulated pathways",
                                   subtitle =  str_replace(contrasts[i], "_vs_", " vs. "), plot = plot, export = export, kegg = FALSE)
            ) )
        }
        if(!(nrow(as.data.frame(res.pathway$ls.pora_custom_dn[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(res.pathway$ls.pora_custom_dn[[i]], title = "Downregulated pathways",
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

  message(paste0("Save results as tab-delimited table to: ", getwd(), "/results.txt"))
  utils::write.table(results$results, paste(getwd(), "results.txt",
                                            sep = "/"), row.names = FALSE, sep = "\t")

  message("Save RData object")
  save(results, file = paste(out_dir, "results.RData", sep = "/"))

  if (pathway_enrichment == T && pathway_kegg) {
    results <- c(results, pora_kegg_up = list(res.pathway$ls.pora_kegg_up), pora_kegg_dn = list(res.pathway$ls.pora_kegg_dn))
    message(paste0("Writing results of KEGG pathway enrichment analysis to: ", out_dir, "/pora_kegg_contrast...txt'"))
    for(i in 1:length(res.pathway$ls.pora_kegg_up)){
      if(!is.null(res.pathway$ls.pora_kegg_up[[i]])){
        utils::write.table(res.pathway$ls.pora_kegg_up[[i]]@result, paste(out_dir, paste0("pora_kegg_", names(res.pathway$ls.pora_kegg_up)[i], "_up.txt"),
                                                                          sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
    for(i in 1:length(res.pathway$ls.pora_kegg_dn)){
      if(!is.null(res.pathway$ls.pora_kegg_dn[[i]])){
        utils::write.table(res.pathway$ls.pora_kegg_dn[[i]]@result, paste(out_dir, paste0("pora_kegg_", names(res.pathway$ls.pora_kegg_dn)[i], "_down.txt"),
                                                                          sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
  }
  if (pathway_enrichment == T && !is.null(custom_pathways)) {
    results <- c(results, pora_custom_up = list(res.pathway$ls.pora_custom_up), pora_custom_dn = list(res.pathway$ls.pora_custom_dn))
    message(paste0("Writing results of custom pathway enrichment analysis to: ", out_dir, "/pora_custom_contrast...txt"))
    for(i in 1:length(res.pathway$ls.pora_custom_up)){
      if(!is.null(res.pathway$ls.pora_custom_up[[i]])){
        utils::write.table(res.pathway$ls.pora_custom_up[[i]]@result, paste(out_dir, paste0("pora_custom_", names(res.pathway$ls.pora_custom_up)[i], "_up.txt"),
                                                                            sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
    for(i in 1:length(res.pathway$ls.pora_custom_dn)){
      if(!is.null(res.pathway$ls.pora_custom_dn[[i]])){
        utils::write.table(res.pathway$ls.pora_custom_dn[[i]]@result, paste(out_dir, paste0("pora_custom_", names(res.pathway$ls.pora_custom_dn)[i], "_down.txt"),
                                                                            sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
  }

  if(report == TRUE){
    prot.report(results,
                volcano.adjusted = volcano.adjusted,
                pathway_enrichment = pathway_enrichment,
                heatmap.show_all = heatmap.show_all,
                heatmap.kmeans = heatmap.kmeans,
                volcano.add_names = volcano.add_names,
                k = k,
                report.dir = report.dir)
  }
  setwd(old_wd)
  return(results)
}
