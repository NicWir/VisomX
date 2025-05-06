## ----library, echo = T, eval = F, results = "hide"----------------------------
# # Load the VisomX package
# library(VisomX)

## ----libraryCode, echo = F, eval = T, include = F, results = "hide"-----------
# Load the VisomX package
suppressPackageStartupMessages(library(VisomX))
suppressPackageStartupMessages(library(dplyr))

## ----SamplePath, echo = T, eval = T, include = T, results = "asis"------------
system.file("prot_KT2440_glc_ac_LB.txt", package = "VisomX")


## ----read, include = T, eval = T----------------------------------------------
# Read data, incl. filtering
data <- prot.read_data(
          data = system.file("prot_KT2440_glc_ac_LB.txt", package = "VisomX"),
          pfx = "abundance.",
          id = "Ensembl Gene ID",
          name = "Gene Symbol",
          filt_type = "fraction",
          filt_min = 0.66)

## ----workflow, include = T, eval = T, results = "asis"------------------------
# Read TXT file with BioCyc pathways for P. putida KT2440
custom_df <-
  read.table(
    system.file("BioCyc_pathways_KT2440.txt", package = "VisomX"),
    sep = "\t",
    header = T,
    stringsAsFactors = F,
    fill = T,
    na.strings = "",
    quote = "",
    comment.char = "",
    check.names = F
  )

# Run workflow and export results
results <- prot.workflow(se = data,
        imp_fun = "SampMin",
        type = "all",
        alpha = 0.05,
        lfc = 1,
        pathway_enrichment = TRUE,
        pathway_kegg = TRUE,
        kegg_organism = "ppu",
        custom_pathways = custom_df,
        volcano.add_names = FALSE, # show protein labels in volcano
        report = FALSE,
        out.dir = tempdir()
        )

## ----CustomWorkflow, include = T, eval = F------------------------------------
# # 1.  Perform data normalization:
# prot_norm <- prot.normalize_vsn(se)
# 
# # 2.  Missing value imputation
# 
# #    a.   with imp_fun == "MinProb":
#         prot_imp <- prot.impute(prot_norm, fun = imp_fun, q = q)
# 
# #    b.  with imp_fun == "knn":
#         prot_imp <- prot.impute(prot_norm, fun = imp_fun, rowmax = knn.rowmax)
# 
# #    c.  for all other methods:
#         prot_imp <- prot.impute(prot_norm, fun = imp_fun)
# 
# # 3.  Principal component analysis
# prot_pca <- prot.pca(SummarizedExperiment::assay(prot_imp))
# 
# # 4.  Test for differential expression by empirical Bayes moderation of a linear model and defined contrasts:
# prot_diff <- prot.test_diff(prot_imp, type = type, control = control, test = contrast)
# 
# # 5.  Denote significantly differential proteins:
# prot_dep <- prot.add_rejections(prot_diff, alpha = alpha, lfc = lfc)
# 
# # 6.  Generate a results table:
# prot_res <- prot.get_results(prot_dep)
# 
# # 7.  Pathway enrichment analysis
# 
# #    a. Get vector of tested contrasts:
#         contrasts <- SummarizedExperiment::rowData(prot_dep) %>%
#           data.frame(check.names = FALSE) %>%
#           select(ends_with("_diff")) %>%
#           colnames() %>% str_replace_all("_diff", "")
# 
# #    b. Perform pathway enrichment analysis
#         res.pathway <- enrich_pathways(prot_dep, contrasts, alpha_pathways = alpha_pathways,
#                                    pathway_kegg = pathway_kegg, kegg_organism = kegg_organism,
#                                    custom_pathways = custom_pathways)
# 
# # 8.  Assemble results into list
#         results <- list(data = SummarizedExperiment::rowData(se), se = se, norm = prot_norm,
#                   imputed = prot_imp, pca = prot_pca, diff = prot_diff, dep = prot_dep,
#                   results = results, param = param)

## ----echo=FALSE---------------------------------------------------------------
pg_width = ncol(data) / 2
if (pg_width > 10) { pg_width = 10 }
suppress_ggrepel <- function(w) {
    if (any(grepl("ggrepel", w)))
      invokeRestart("muffleWarning")
  }

## ----numbers, echo=FALSE, warning=FALSE, fig.height = 4, fig.width = pg_width, fig.align='center'----
prot.plot_numbers(data, 
                  plot = T, 
                  export = F)

## ----coverage, echo=FALSE, warning=FALSE, fig.height = 4, fig.width = 6, fig.align='center'----
plot_coverage(data)

## ----meanSDplot1,echo=FALSE, warning=FALSE, fig.width=5, fig.height=3, fig.align='center', warning=F----
  meanSdPlot(data, 
             ranks = TRUE, 
             plot = T, 
             xlab = "Rank(mean)", 
             ylab = "SD")

## ----meanSDplot2,echo=FALSE, warning=FALSE, fig.width=5, fig.height=3, fig.align='center', warning=F----
  meanSdPlot(results$norm, 
             plot = T, 
             xlab = "Rank(mean)", 
             ylab = "SD")

## ----echo=FALSE---------------------------------------------------------------
norm_height = ncol(data) / 2.2
if (norm_height > 9.5) {norm_height = 9.5}

## ----norm, echo=FALSE, warning=FALSE, fig.height = norm_height, fig.align='center'----
prot.boxplot_intensity(data, 
                       results$norm, 
                       plot = T, 
                       export = F)

## ----missval, echo=FALSE, fig.width = 7, fig.height = 6, fig.align='center'----
prot.plot_missval(results$norm, 
                  fontsize = 12, 
                  export = F)

## ----missval2, echo=FALSE, fig.width = 7, fig.height = 7, fig.align='center'----
prot.plot_detect(results$norm,
                 basesize = 10,
                 plot = T,
                 export = F)

## ----density, echo=FALSE, fig.width = 7, fig.height = 9, fig.align='center'----
prot.plot_imputation(
  data,
  results$norm,
  results$imp,
  basesize = 13,
  plot = T,
  export = F
)

## ----screeplot, echo=FALSE, fig.width = 7, fig.height= 6, warning=F, message=F, fig.align='center'----
prot.plot_screeplot(
  results$pca,
  axisLabSize = 16,
  titleLabSize = 17,
  plot = T,
  export = F,
  components = PCAtools::getComponents(pca)[1:ncomp]
)

## ----plot_pca_1_2, echo=FALSE, fig.width = 7, fig.height= 6, warning=F, message=F, fig.align='center'----
prot.plot_pca(
  results$imp,
  x = 1,
  y = 2,
  point_size = 4,
  basesize = 12,
  plot = T,
  export = F,
  title = "PC Scores - PC2 vs. PC1"
)

## ----plot_pca_1_3, echo=FALSE, fig.width = 7, fig.height= 6, warning=F, message=F, fig.align='center'----
prot.plot_pca(
  results$imp,
  x = 1,
  y = 3,
  point_size = 4,
  basesize = 12,
  plot = T,
  export = F,
  title = "PC Scores - PC3 vs. PC1"
)

## ----loadingsPlot, echo=FALSE, fig.width = 7, fig.height= 8, warning=F, message=F, fig.align='center'----
prot.plot_loadings(results$pca,labSize = 3, plot = T, export = F)

