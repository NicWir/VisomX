## ----library, echo = T, eval = F, results = "hide"----------------------------
#  # Load the VisomX package
#  library(VisomX)

## ----libraryCode, echo = F, eval = T, include = F, results = "hide"-----------
# Load the VisomX package
suppressPackageStartupMessages(library(VisomX))

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
#  # 1.  Perform data normalization:
#  prot_norm <- prot.normalize_vsn(se)
#  
#  # 2.  Missing value imputation
#  
#  #    a.   with imp_fun == "MinProb":
#          prot_imp <- prot.impute(prot_norm, fun = imp_fun, q = q)
#  
#  #    b.  with imp_fun == "knn":
#          prot_imp <- prot.impute(prot_norm, fun = imp_fun, rowmax = knn.rowmax)
#  
#  #    c.  for all other methods:
#          prot_imp <- prot.impute(prot_norm, fun = imp_fun)
#  
#  # 3.  Principal component analysis
#  prot_pca <- prot.pca(SummarizedExperiment::assay(prot_imp))
#  
#  # 4.  Test for differential expression by empirical Bayes moderation of a linear model and defined contrasts:
#  prot_diff <- prot.test_diff(prot_imp, type = type, control = control, test = contrast)
#  
#  # 5.  Denote significantly differential proteins:
#  prot_dep <- prot.add_rejections(prot_diff, alpha = alpha, lfc = lfc)
#  
#  # 6.  Generate a results table:
#  prot_res <- prot.get_results(prot_dep)
#  
#  # 7.  Pathway enrichment analysis
#  
#  #    a. Get vector of tested contrasts:
#          contrasts <- SummarizedExperiment::rowData(prot_dep) %>%
#            data.frame(check.names = FALSE) %>%
#            select(ends_with("_diff")) %>%
#            colnames() %>% str_replace_all("_diff", "")
#  
#  #    b. Perform pathway enrichment analysis
#          res.pathway <- enrich_pathways(prot_dep, contrasts, alpha_pathways = alpha_pathways,
#                                     pathway_kegg = pathway_kegg, kegg_organism = kegg_organism,
#                                     custom_pathways = custom_pathways)
#  
#  # 8.  Assemble results into list
#          results <- list(data = SummarizedExperiment::rowData(se), se = se, norm = prot_norm,
#                    imputed = prot_imp, pca = prot_pca, diff = prot_diff, dep = prot_dep,
#                    results = results, param = param)

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

## ---- echo=FALSE--------------------------------------------------------------
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

## ----correlation_heatmap_all, echo=FALSE, fig.width = 7, fig.height= 7, warning=F, message=F, fig.align='center'----
prot.plot_corrheatmap(
  results$dep,
  lower = 0.8,
  upper = 1,
  font_size = 12,
  significant = FALSE
)

## ----correlation_heatmap, echo=FALSE, fig.width = 7, fig.height= 7, warning=F, message=F, fig.align='center'----
prot.plot_corrheatmap(
  results$dep,
  lower = 0.8,
  upper = 1,
  font_size = 12,
  significant = T
)

## ----echo=FALSE, warning=F, message=F-----------------------------------------
width = ncol(data) / 2
if (width < 7) { width = 7 }
if (width > 12) { width = 12 }
len = nrow(results$dep[SummarizedExperiment::rowData(results$dep)$significant, ]) / 13
if (len < 5) { len = 5 }

## ----heatmap_1, echo=FALSE, fig.width = width, fig.height = len, warning=F, message=F, fig.align='center'----
prot.plot_heatmap(
  results$dep,
  type = "contrast",
  kmeans = TRUE,
  k = 4,
  show_row_names = T,
  row_font_size = 5,
  plot = T,
  export = F
)

## ----heatmap_2, echo=FALSE, fig.width = width*1.2, fig.height = len, fig.align='center'----
prot.plot_heatmap(
  results$dep,
  type = "centered",
  kmeans = TRUE,
  k = 4,
  show_all = T,
  show_row_names = T,
  row_font_size = 5,
  indicate = "condition",
  col_limit = 2,
  plot = T,
  export = F
)

## ----volcano, echo=FALSE, fig.width = 8, fig.height = 8.8, warning=F, message=F, fig.align='center', warning = FALSE----
prot.plot_volcano(
  results$dep,
  contrast = "Glc_vs_Ac",
  add_names = TRUE,
  label_size = 2.5,
  adjusted =  FALSE,
  plot = TRUE,
  export = FALSE,
  lfc = 1,
  alpha = 0.05
)

## ----echo=FALSE---------------------------------------------------------------
w_pora <- 11
h_pora <- 7

## ----enrichment_kegg, echo=FALSE, fig.width = w_pora, fig.height = h_pora, warning=F, message=F, results='asis'----
  # Lollipop chart
  prot.plot_enrichment(
    results$pora_kegg_up[["Glc_vs_Ac"]],
    title = "Upregulated pathways - KEGG",
    subtitle =  "Glc_vs_Ac",
    plot = T,
    export = F
  )

  # UpSet plot
  prot.plot_upset(results$pora_kegg_up[["Glc_vs_Ac"]],
                  line.size = 0.9,
                  mb.ratio = c(0.6, 0.4)
                  )


## ----echo=T-------------------------------------------------------------------
df <- as.data.frame(results$pora_kegg_up[["Glc_vs_Ac"]])
# Inspect the structure of the dataframe
glimpse(df)

## ----echo=T-------------------------------------------------------------------
df <- as.data.frame(results$pora_kegg_dn[["Glc_vs_Ac"]])
# Inspect the structure of the dataframe
glimpse(df)

## ----echo=FALSE---------------------------------------------------------------
w_pora <- 11
h_pora <- 9

## ----enrichment_custom_up, echo=FALSE, fig.width = w_pora, fig.height = h_pora, warning=F, message=F, results='asis'----
# Plot overrepresented KEGG pathways
  ## Lollipop chart
  prot.plot_enrichment(
    results$pora_custom_up[["Glc_vs_Ac"]],
    title = "Upregulated pathways - KEGG",
    subtitle =  "Glc_vs_Ac",
    plot = T,
    export = F
  )
  ## UpSet plot
  prot.plot_upset(results$pora_custom_up[["Glc_vs_Ac"]],
                  line.size = 0.9,
                  mb.ratio = c(0.6, 0.4)
                  )

## ----echo=FALSE---------------------------------------------------------------
w_pora <- 11
h_pora <- 6

## ----enrichment_custom_dn, echo=FALSE, fig.width = w_pora, fig.height = h_pora, warning=F, message=F, results='asis'----
# Plot underrepresented Custom pathways
  ## Lollipop chart
  prot.plot_enrichment(
    results$pora_custom_dn[["Glc_vs_Ac"]],
    title = "Upregulated pathways - KEGG",
    subtitle =  "Glc_vs_Ac",
    plot = T,
    export = F
  )
  ## UpSet plot
  prot.plot_upset(results$pora_custom_dn[["Glc_vs_Ac"]],
                  line.size = 0.9,
                  mb.ratio = c(0.6, 0.4)
                  )

## ----PlotBar, echo=FALSE, fig.width = 13, fig.height = 8, warning=F, message=F, results='asis'----
# Load the SummarizedExperiment package
library(SummarizedExperiment)

# Get all KEGG pathways and their respective genes
kegg_pathways <- prot.get_kegg_pathways("ppu")

# Extract genes constituting the TCA cycle
TCA_genes <- kegg_pathways[["Citrate cycle (TCA cycle)"]]

# Get the names instead of IDs for the genes
TCA_names <- rowData(results$dep)[rowData(results$dep)$ID %in% TCA_genes, "name"]

# Plot TCA cycle protein abundances
prot.plot_bar(
  dep = results$dep,
  type = "centered",
  combine = T,
  proteins = TCA_names,
  col.id = "name",
  export = F,
  plot = T
)

## ----PlotBarFC, echo=FALSE, fig.width = 13, fig.height = 8, warning=F, message=F, results='asis'----

# Plot TCA cycle protein abundances
prot.plot_bar(
  dep = results$dep,
  type = "contrast",
  contrast = "Glc_vs_Ac",
  combine = T,
  convert_name = T,
  proteins = TCA_names,
  col.id = "name",
  export = F,
  plot = T
)

