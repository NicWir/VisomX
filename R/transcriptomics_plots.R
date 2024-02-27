#' Plot a correlation heatmap
#'
#' This function plots a Pearson correlation heatmap for a given SummarizedExperiment object.
#'
#' @param dep A DESeqDataSet object.
#' @param significant Logical. If TRUE, only significant correlations are plotted. Default is TRUE.
#' @param lower Numeric. Lower limit of the color scale. Default is 0.
#' @param upper Numeric. Upper limit of the color scale. Default is 1.
#' @param pal Character. Brewer color palette to be used. Default is "PRGn". See \code{\link[RColorBrewer]{brewer.pal.info}} for options
#' @param pal_rev Logical. If TRUE, color palette is reversed. Default is FALSE.
#' @param indicate Character. Column name from colData to use for annotations. Default is NULL.
#' @param font_size Numeric. Font size for labels. Default is 12.
#' @param plot Logical. If FALSE, no plot is generated, only data frame is returned. Default is TRUE.
#' @param ... Other arguments to be passed to \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @return (Invisibly) a data frame with Pearson correlation values is returned.
#'
#' @export
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom assertthat assert_that
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
#' @importFrom tibble rownames_to_column
#' @importFrom SummarizedExperiment assay colData rowData
#'
rna.plot_corrheatmap <- function (dds, lower = 0, upper = 1, pal = "PRGn",
                                   pal_rev = FALSE, indicate = NULL, font_size = 12, plot = TRUE,
                                   ...)
{
  assertthat::assert_that(inherits(dds, "DESeqDataSet"),
                          is.logical(significant), length(significant) == 1, is.numeric(lower),
                          length(lower) == 1, is.numeric(upper), length(upper) ==
                            1, is.character(pal), length(pal) == 1, is.logical(pal_rev),
                          length(pal_rev) == 1, is.numeric(font_size), length(font_size) ==
                            1, is.logical(plot), length(plot) == 1)
  if (!(lower >= -1 && upper >= -1 && lower <= 1 && upper <= 1)) {
    stop("'lower' and/or 'upper' arguments are not valid\n         Run rna.plot_corrheatmap() with 'lower' and 'upper' between -1 and 1",
         call. = FALSE)
  }
  pals <- RColorBrewer::brewer.pal.info %>% tibble::rownames_to_column() %>%
    filter(category != "qual")
  if (!pal %in% pals$rowname) {
    stop("'", pal, "' is not a valid color panel",
         " (qualitative panels also not allowed)\n",
         "Run plot_pca() with one of the following 'pal' options: ",
         paste(pals$rowname, collapse = "', '"), "'",
         call. = FALSE)
  }
  if (!is.null(indicate)) {
    assertthat::assert_that(is.character(indicate))
    col_data <- SummarizedExperiment::colData(dds) %>% as.data.frame()
    columns <- colnames(col_data)
    if (any(!indicate %in% columns)) {
      stop("'", paste0(indicate, collapse = "' and/or '"),
           "' column(s) is/are not present in ", ddsarse(substitute(dds)),
           ".\nValid columns are: '", paste(columns,
                                            collapse = "', '"), "'.", call. = FALSE)
    }
    anno <- SummarizedExperiment::colData(dds) %>% data.frame() %>% select(indicate)
    names <- colnames(anno)
    anno_col <- vector(mode = "list", length = length(names))
    names(anno_col) <- names
    for (i in names) {
      var = anno[[i]] %>% unique() %>% sort()
      if (length(var) == 1)
        cols <- c("black")
      if (length(var) == 2)
        cols <- c("orangered", "cornflowerblue")
      if (length(var) < 7 & length(var) > 2)
        cols <- RColorBrewer::brewer.pal(length(var),
                                         "Pastel1")
      if (length(var) >= 7)
        cols <- RColorBrewer::brewer.pal(length(var),
                                         "Set3")
      names(cols) <- var
      anno_col[[i]] <- cols
    }
    ha1 = HeatmapAnnotation(df = anno, col = anno_col, show_annotation_name = TRUE)
  } else {
    ha1 <- NULL
  }
  cor_mat <- cor(SummarizedExperiment::assay(dds))
  ht1 = ComplexHeatmap::Heatmap(cor_mat, col = circlize::colorRamp2(seq(lower,
                                                                        upper, ((upper - lower)/7)), if (pal_rev) {
                                                                          rev(RColorBrewer::brewer.pal(8, pal))
                                                                        }
                                                                    else {
                                                                      RColorBrewer::brewer.pal(8, pal)
                                                                    }), heatmap_legend_param = list(color_bar = "continuous",
                                                                                                    legend_direction = "horizontal", legend_width = unit(5,
                                                                                                                                                         "cm"), title_position = "topcenter"),
                                name = "Pearson correlation", column_names_gp = grid:::gpar(fontsize = font_size),
                                row_names_gp = grid:::gpar(fontsize = font_size), top_annotation = ha1,
                                ...)
  if (plot) {
    ComplexHeatmap::draw(ht1, heatmap_legend_side = "top")
  }
  else {
    df <- as.data.frame(cor_mat)
    return(df)
  }
}

#' Plots a heatmap of the differentially expressed genes
#'
#' @param dds A DESeqDataSet object containing the normalized gene abundances and the associated metadata
#' @param type Type of heatmap to plot. Either "centered" or "contrast".
#' @param contrast (String or vector of strings) Analyze only significant genes contained in the specified contrast(s). Default is NULL. Contrasts must be given in the form "ConditionA_vs_ConditionB"
#' @param show_all Show gene abundances of all conditions or only of those in the specified constrast(s). Default is TRUE
#' @param pal Color palette to use. Default is "RdBu"
#' @param kmeans Perform k-means clustering on the rows. Default is FALSE
#' @param k Number of clusters to use when kmeans is TRUE. Default is 6
#' @param col_limit Specify a custom range for the color scale of the heatmap. Default is NA, in which case the minimum and maximum of the color scale are calculated as the 5% and 95% percentiles.
#' @param indicate Annotate the heatmap with the specified column from the metadata.
#' @param clustering_distance Distance metric used for clustering. Default is "euclidean".
#' @param show_row_names Show the gene names on the left side of the heatmap. Default is FALSE.
#' @param row_font_size Font size for the gene names. Default is 6
#' @param col_font_size Font size for the condition names. Default is 10
#' @param plot Plot the heatmap. Default is TRUE
#' @param export Export the heatmap as pdf and png. Default is TRUE
#' @param ... Other parameters passed to the ComplexHeatmap::Heatmap function
#' @return (Invisibly) returns a data frame with the genes and the normalized abundances
#' @export
rna.plot_heatmap <- function (dds,
                              type = c("centered", "contrast"),
                              contrast = NULL, # Analyze only genes contained in the specified contrast(s)
                              show_all = TRUE, # Show gene abundances of all conditions or only of those in the specified constrast(s)
                              pal = "RdBu",
                              kmeans = FALSE,
                              k = 6,
                              col_limit = NA,
                              indicate = NULL,
                              clustering_distance = c("euclidean",
                                                      "maximum", "manhattan", "canberra",
                                                      "binary", "minkowski", "pearson", "spearman",
                                                      "kendall", "gower"),
                              show_row_names = FALSE,
                              row_font_size = 6,
                              col_font_size = 10,
                              plot = TRUE,
                              export = TRUE, ...)
{
  if(is.null(contrast) & show_all == FALSE){
    contrast <- SummarizedExperiment::rowData(dds, use.names = FALSE) %>% data.frame(check.names = FALSE) %>% select(starts_with("lfc."))  %>%
      colnames(.) %>% gsub("lfc.", "", .)
  }
  if (length(grep(paste("lfc.", paste(contrast, collapse = "|"), sep = ""),
                  colnames(SummarizedExperiment::rowData(dds, use.names = FALSE)))) == 0) {
    valid_cntrsts <-
      SummarizedExperiment::rowData(dds, use.names = FALSE) %>% data.frame(check.names = FALSE) %>% select(starts_with("lfc.")) %>%
      colnames(.) %>% gsub("diff.", "", .)
    valid_cntrsts_msg <- paste0("Valid contrasts are: '",
                                paste0(valid_cntrsts, collapse = "', '"),
                                "'")
    stop(
      "Not a valid contrast, please run `rna.plot_heatmap()` with a valid contrast as argument\n",
      valid_cntrsts_msg,
      call. = FALSE
    )
  }
  if (!is.null(contrast) && show_all == FALSE){
    contrast_samples <- contrast %>% str_split("_vs_") %>% unlist() %>% unique()
  }
  if (is.integer(k)) k <- as.numeric(k)
  if (is.integer(row_font_size)) row_font_size <- as.numeric(row_font_size)
  if (is.integer(col_font_size)) col_font_size <- as.numeric(col_font_size)
  assertthat::assert_that(inherits(dds, "SummarizedExperiment"),
                          is.character(type), is.logical(kmeans), is.numeric(k),
                          length(k) == 1, is.numeric(row_font_size), length(row_font_size) ==
                            1, is.numeric(col_font_size), length(col_font_size) ==
                            1, is.logical(plot), length(plot) == 1)
  type <- match.arg(type)
  clustering_distance <- match.arg(clustering_distance)
  row_data <- SummarizedExperiment::rowData(dds, use.names = FALSE)
  col_data <- SummarizedExperiment::colData(dds) %>% data.frame(check.names = FALSE)
  if (any(!c("label", "condition", "replicate") %in%
          colnames(col_data))) {
    stop(paste0("'label', 'condition' and/or 'replicate' columns are not present in '",
                deparse(substitute(dds)), "'"), call. = FALSE)
  }
  if (length(grep("^lfc.", colnames(row_data))) < 1) {
    stop(paste0("'lfc.[...]'", "columns are not present in '",
                deparse(substitute(dds)), "'.\nRun rna.workflow() to obtain the required columns."),
         call. = FALSE)
  }
  if (!"significant" %in% colnames(row_data)) {
    stop(paste0("'significant' column is not present in '",
                deparse(substitute(dds)), "'.\nRun rna.workflow() to obtain the required column."),
         call. = FALSE)
  }
  if (!is.null(indicate) && type == "contrast") {
    warning("Heatmap annotation only applicable for type = 'centered'",
            call. = FALSE)
  }
  if (!is.null(indicate) && type == "centered" && show_all == TRUE) {
    ha1 <- get_annotation(dds, indicate)
  } else if (!is.null(indicate) && type == "centered" && show_all == FALSE) {
    ha1 <- get_annotation_contrast(dds, indicate, contrast = contrast_samples)
  } else {
    ha1 <- NULL
  }
  filtered <- dds[row_data$significant[!is.na(row_data$significant)], ]
  if (nrow(SummarizedExperiment::assay(filtered)) == 0){
    stop("No genes with significantly different abundance were found.")
  }
  if (any(is.na(SummarizedExperiment::assay(filtered)))) {
    warning("Missing values in '", deparse(substitute(dds)),
            "'. ", "Using clustering_distance = 'gower'",
            call. = FALSE)
    clustering_distance <- "gower"
    obs_NA <- TRUE
  } else {
    obs_NA <- FALSE
  }
  if (type == "centered") {
    SummarizedExperiment::rowData(filtered)$mean <- log2(rowMeans(SummarizedExperiment::assay(filtered)+1, na.rm = TRUE))
    df <- log2(SummarizedExperiment::assay(filtered)+1) - SummarizedExperiment::rowData(filtered, use.names = FALSE)$mean
    if (!is.null(contrast) && show_all == FALSE){
      df <- data.frame(df, check.names = FALSE)[,str_detect(colnames(df),paste(contrast_samples, collapse = "|"))]
    }
  }
  if (type == "contrast") {
    lfc.pfx <- ifelse(length(grep("lfc_shrink", colnames(SummarizedExperiment::rowData(dds))))>0,
                      "lfc_shrink\\.", "lfc\\.")
    df <- SummarizedExperiment::rowData(filtered, use.names = FALSE) %>% data.frame(check.names = FALSE) %>%
      tibble::column_to_rownames(var = "name") %>% select(starts_with(gsub("\\\\.", ".", lfc.pfx)))
    if(length(grep(lfc.pfx, colnames(SummarizedExperiment::rowData(filtered)))) == 1){
      df[[1]] <- df[[1]] %>% sort()
    }
    colnames(df) <- gsub(ifelse(length(grep("lfc_shrink", colnames(SummarizedExperiment::rowData(dds))))>0,
                                "lfc_shrink\\.", "lfc."), "", colnames(df)) %>%
      gsub("_vs_", " vs \n", .)

  }
  df <- as.matrix(df)
  if (kmeans && obs_NA) {
    warning("Cannot perform kmeans clustering with missing values",
            call. = FALSE)
    kmeans <- FALSE
  }
  if (kmeans && !obs_NA) {
    set.seed(1)
    df_kmeans <- kmeans(df, k)
    if (type == "centered") {
      order <- data.frame(df, check.names = FALSE) %>% cbind(., cluster = df_kmeans$cluster) %>%
        mutate(row = apply(.[, seq_len(ncol(.) - 1)],
                           1, function(x) max(x))) %>% group_by(cluster) %>%
        dplyr::summarize(index = sum(row)/n()) %>% arrange(desc(index)) %>%
        pull(cluster) %>% match(seq_len(k), .)
      df_kmeans$cluster <- order[df_kmeans$cluster]
    }
    if (type == "contrast") {
      order <- data.frame(df, check.names = FALSE) %>% cbind(., cluster = df_kmeans$cluster) %>%
        gather(condition, diff, -cluster) %>% group_by(cluster) %>%
        dplyr::summarize(row = mean(diff)) %>% arrange(desc(row)) %>%
        pull(cluster) %>% match(seq_len(k), .)
      df_kmeans$cluster <- order[df_kmeans$cluster]
    }
  }
  if (ncol(df) == 1) {
    col_clust = FALSE
  } else {
    col_clust = TRUE
  }
  if (nrow(df) == 1) {
    row_clust = FALSE
  } else {
    row_clust = TRUE
  }
  if (clustering_distance == "gower") {
    clustering_distance <- function(x) {
      dist <- cluster::daisy(x, metric = "gower")
      dist[is.na(dist)] <- max(dist, na.rm = TRUE)
      return(dist)
    }
  }
  col_limit = if(!is.na(col_limit)){
    col_limit
  } else {
    df_cent <- log2(SummarizedExperiment::assay(filtered)+1) - log2(rowMeans(SummarizedExperiment::assay(filtered)+1))
    col_lim <- ceiling(stats::quantile(df_cent, probs= 0.95, na.rm = TRUE)) - ceiling(stats::quantile(df_cent, probs= 0.05, na.rm = TRUE))
  }
  breaks =  if (col_limit == 1){
    seq(-1, 1, 0.5)
  } else if (1 < col_limit && col_limit <= 4){
    seq(-col_limit, col_limit, 1)
  } else if (4 < col_limit && col_limit <= 6){
    seq(-6, 6, 2)
  }else if (6 < col_limit && col_limit <= 8){
    seq(-8,8,2)
  } else if (8 < col_limit && col_limit <= 10){
    seq(-10,10,5)
  } else {
    seq(-plyr::round_any(col_limit, 5, f=ceiling),plyr::round_any(col_limit, 5, f=ceiling),5)
  }
  legend <- ifelse(type == "contrast", "log2 Fold change",
                   "log2 Centered intensity")
  col_fun = circlize::colorRamp2(seq(-col_limit,
                                     col_limit, (col_limit /
                                                   5)),
                                 rev(RColorBrewer::brewer.pal(11, pal)), space = "LAB")
  ht1 = ComplexHeatmap::Heatmap(
    df,
    col = col_fun,
    split = if (kmeans) {
      df_kmeans$cluster
    } else {
      NULL
    },
    cluster_rows = col_clust,
    cluster_columns = row_clust,
    row_names_side = "left",
    column_names_side = "top",
    clustering_distance_rows = clustering_distance,
    clustering_distance_columns = clustering_distance,
    heatmap_legend_param = list(
      at = breaks,
      color_bar = "continuous",
      legend_direction = "horizontal",
      legend_width = unit(5,
                          "cm"),
      title_position = "lefttop",
      border = "black"
    ),
    name = legend,
    row_names_gp = grid::gpar(fontsize = row_font_size),
    column_names_gp = grid::gpar(fontsize = col_font_size),
    border_gp = (grid::gpar(col = "black", lty = 1)),
    top_annotation = ha1,
    show_row_names = show_row_names,
    ...
  )
  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
    w <- 6+length(dds$condition)/10
    if (show_row_names == TRUE){
      len = nrow(filtered) / 12
    } else {
      len = nrow(filtered) / 45
    }
    if (len < 6) { len = 6 }
    message(paste0("Exporting heat map to:\n\"",
                   getwd(),
                   "/Plots/HeatMap_",
                   type,
                   ".pdf",
                   " and \"HeatMap_",
                   type,
                   ".png\""))

    grDevices::pdf(paste0("Plots/HeatMap_",
                          type,
                          ".pdf"),
                   width = w, height = len)
    ComplexHeatmap::draw(ht1, heatmap_legend_side = "top")
    grDevices::dev.off()

    grDevices::png(paste0("Plots/HeatMap_",
                          type,
                          ".png"),
                   width = w, height = len, units = 'in', res = 300)
    ComplexHeatmap::draw(ht1, heatmap_legend_side = "top")
    grDevices::dev.off()
  }
  if (plot) {
    ComplexHeatmap::draw(ht1, heatmap_legend_side = "top")
  } else {
    colnames(df) <- gsub(" ", "_", colnames(df))
    suppressWarnings(
      df <- df[, unlist(ComplexHeatmap::column_order(ht1))]
    )
    if (kmeans) {
      suppressWarnings(
        df <- cbind(df, k = df_kmeans$cluster)
      )
    }
    try(suppressWarnings(
      return <- df[unlist(ComplexHeatmap::row_order(ht1)), ]))
    data.frame(gene = row.names(return), return, check.names = FALSE) %>% mutate(order = row_number())
  }
}

#' Plot Results of Principal Component Analysis
#'
#' Performs 'regularized log' transformation followed by PCA on a given DESeqDataSet object and plots the results.
#'
#' @param dep DESeqDataSet object
#' @param x x-axis PC (default: 1)
#' @param y y-axis PC (default: 2)
#' @param indicate features to indicate in the plot ("condition" or "replicate")
#' @param title title of the plot (default: "PCA plot - top \code{n} variable genes")
#' @param label whether to label points in the plot (default: FALSE)
#' @param n number of variables to take into account, sorted according to their variance in descending order. Only the \code{n} most variable genes are used to perform PCA. default: number of columns in \code{dep})
#' @param point_size size of points in the plot (default: 4)
#' @param label_size size of labels in the plot (default: 3)
#' @param hline position of horizontal line (default: 0)
#' @param hlineType type of horizontal line (default: 'longdash')
#' @param hlineCol color of horizontal line (default: 'black')
#' @param hlineWidth width of horizontal line (default: 0.4)
#' @param vline position of vertical line (default: 0)
#' @param vlineType type of vertical line (default: 'longdash')
#' @param vlineCol color of vertical line (default: 'black')
#' @param vlineWidth width of vertical line (default: 0.4)
#' @param basesize base size of the plot (default: 15)
#' @param plot whether to return the PCA plot (default: TRUE)
#' @param export whether to export the PCA plot to the Plots directory as PNG and PDF file (default: TRUE)
#'
#' @return (invisibly) a data frame containing the PCA coordinates
#'
#' @importFrom DESeq2 rlog
#' @importFrom SummarizedExperiment assay colData
#'
#' @export
rna.plot_pca <- function (dds,
                           x = 1,
                           y = 2,
                           indicate = c("condition", "replicate"),
                           title = NULL,
                           label = FALSE,
                           n = ncol(t(SummarizedExperiment::assay(dds))),
                           point_size = 4,
                           label_size = 3,
                           hline = 0,
                           hlineType = 'longdash',
                           hlineCol = 'black',
                           hlineWidth = 0.4,
                           vline = 0,
                           vlineType = 'longdash',
                           vlineCol = 'black',
                           vlineWidth = 0.4,
                           basesize = 15,
                           plot = TRUE,
                           export = TRUE)
{
  if (is.integer(x))
    x <- as.numeric(x)
  if (is.integer(y))
    y <- as.numeric(y)
  if (is.integer(n))
    n <- as.numeric(n)
  if (is.integer(point_size))
    point_size <- as.numeric(point_size)
  if (is.integer(label_size))
    label_size <- as.numeric(label_size)
  assertthat::assert_that(
    inherits(dds, "SummarizedExperiment"),
    is.numeric(x),
    length(x) == 1,
    is.numeric(y),
    length(y) ==
      1,
    is.numeric(n),
    length(n) == 1,
    is.character(indicate),
    is.logical(label),
    is.numeric(point_size),
    length(point_size) ==
      1,
    is.numeric(label_size),
    length(label_size) ==
      1,
    is.logical(plot),
    length(plot) == 1
  )
  if (x > ncol(dds) || y > ncol(dds)) {
    stop(
      paste0(
        "'x' and/or 'y' arguments are not valid\n",
        "Run rna.plot_pca() with 'x' and 'y' <= ",
        ncol(dds),
        "."
      ),
      call. = FALSE
    )
  }
  if (n > nrow(dds)) {
    stop(
      paste0(
        "'n' argument is not valid.\n",
        "Run rna.plot_pca() with 'n' <= ",
        nrow(dds),
        "."
      ),
      call. = FALSE
    )
  }
  columns <- colnames(SummarizedExperiment::colData(dds))
  if (!is.null(indicate)) {
    if (length(indicate) > 3) {
      stop(
        "Too many features in 'indicate'\n        Run rna.plot_pca() with a maximum of 3 indicate features"
      )
    }
    if (any(!indicate %in% columns)) {
      stop(
        paste0(
          "'",
          paste0(indicate, collapse = "' and/or '"),
          "' column(s) is/are not present in ",
          deparse(substitute(dds)),
          ".\nValid columns are: '",
          paste(columns, collapse = "', '"),
          "'."
        ),
        call. = FALSE
      )
    }
  }
  rlog.counts <- tryCatch(DESeq2::rlog(dds, fitType = 'mean'), error = function(e) { rlog(dds, fitType = 'mean') })
  var <- apply(SummarizedExperiment::assay(rlog.counts), 1, sd)
  df <- SummarizedExperiment::assay(rlog.counts)[order(var, decreasing = TRUE)[seq_len(n)],]
  pca <- stats::prcomp(t(df), scale = FALSE)
  pca_df <- pca$x %>% data.frame(check.names = FALSE) %>% tibble::rownames_to_column() %>%
    left_join(., data.frame(SummarizedExperiment::colData(rlog.counts), check.names = FALSE), by = c(rowname = "ID"))
  percent <- round(100 * pca$sdev ^ 2 / sum(pca$sdev ^ 2), 1)
  if (!is.null(title)){
    title <- as.character(title)
  } else {
    title <- paste0("PCA plot - top ", n, " variable genes")
  }
  for (feat in indicate) {
    pca_df[[feat]] <- as.factor(pca_df[[feat]])
  }
  p <- ggplot(pca_df, aes(get(paste0("PC", x)), get(paste0("PC",
                                                           y)))) + labs(
                                                             title = title,
                                                             x = paste0("PC", x, ": ", percent[x], "%"),
                                                             y = paste0("PC", y, ": ", percent[y], "%")
                                                           ) + coord_fixed() + theme_DEP1(basesize = basesize)

  if (length(indicate) == 0) {
    p <- p + geom_point(size = point_size)
  }
  if (length(indicate) == 1) {
    p <- p + geom_point(aes(fill=.data[[indicate[1]]]), col = "black", size = point_size, shape = 21) +
      labs(fill = indicate[1])
  }
  if (length(indicate) == 2) {
    p <- p + geom_point(aes(col = .data[[indicate[1]]],
                            shape = .data[[indicate[2]]], fill=.data[[indicate[1]]]), size = point_size) +
      labs(col = indicate[1], shape = indicate[2]) +
      scale_shape_manual(values=c(21, 22, 23, 24, 25, 0, 1, 2, 5, 6))

    if (length(unique(rlog.counts$condition)) <= 8) {
      p <- p + scale_fill_brewer(palette = "Dark2", guide = guide_legend(override.aes = list(shape = 21)))+
        scale_color_manual(values=c(rep("black", length(unique(rlog.counts$condition)))))
    } else if (length(unique(rlog.counts$condition)) <= 12) {
      p <- p + scale_fill_brewer(palette = "Set3", guide = guide_legend(override.aes = list(shape = 21))) +
        scale_color_manual(values=c(rep("black", length(unique(rlog.counts$condition)))))
    } else {
      pal <- c(
        "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
        "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
        "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
        "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
        "green1", "yellow4", "yellow3", "darkorange4", "brown")
      p <- p + scale_fill_manual(values = pal[1:length(unique(rlog.counts$condition))],
                                 guide = guide_legend(override.aes = list(shape = 21))) +
        scale_color_manual(values=c(rep("black", length(unique(rlog.counts$condition)))))
    }


  }
  if (length(indicate) == 3) {
    p <- p + geom_point(aes(col = .data[[indicate[1]]],
                            shape = .data[[indicate[2]]], fill=.data[[indicate[1]]]), size = point_size) +
      facet_wrap( ~ .data[[indicate[3]]])
    labs(col = indicate[1], shape = indicate[2])+
      scale_shape_manual(values=c(15, 16, 17, 18, 19, 0, 1, 2, 5, 6))
  }
  if (label) {
    p <- p + geom_text(aes(label = rowname), size = label_size)
  }
  # add elements to the plot for vlines and hlines
  if (!is.null(vline)) {
    p <- p + geom_vline(
      xintercept = vline,
      linetype = vlineType,
      colour = vlineCol,
      size = vlineWidth
    )
  }
  if (!is.null(hline)) {
    p <- p + geom_hline(
      yintercept = hline,
      linetype = hlineType,
      colour = hlineCol,
      size = hlineWidth
    )
  }
  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
    message(paste0("Exporting PCA plot to:\n\"", getwd(),
                   "/Plots/PCAPlot_", "PC", x, "-PC", y, ".pdf\" and \".../PCAPlot_",
                   "PC", x, "-PC", y, ".png\""))

    grDevices::pdf(paste0("Plots/PCAPlot_", "PC", x, "-PC", y, ".pdf"))
    print(p)
    grDevices::dev.off()

    grDevices::png(paste0("Plots/PCAPlot_", "PC", x, "-PC", y, ".png"),
                   width = 8, height = 8, units = 'in', res = 300)
    print(p)
    grDevices::dev.off()
  }
  if (plot) {
    print(p)
  }
  else {
    df <- pca_df %>% select(rowname, paste0("PC", c(x, y)),
                            match(indicate, colnames(pca_df)))
    colnames(df)[1] <- "sample"
    return(df)
  }
}

#' Plot imputation results
#'
#' This function takes a matrix or array of imputation results and plots the distribution of intensities. It is recommended to log-transform the data before using this function.
#'
#' @param assay A matrix or array of imputation results. If working with a DESeqDataSet or SummarizedExperiment object, the array can be retrieved via \code{SummarizedExperiment::assay(object)}.
#' @param ... Additional matrices or arrays. Can be input as named argument, for example: \code{"matrix2" = log(matrix2)}
#' @param colData A data frame of sample metadata.  If working with a DESeqDataSet or SummarizedExperiment object, the metadata can be retrieved via \code{SummarizedExperiment::colData(object)}.
#' @param plot A logical indicating whether to plot the results.
#' @param basesize Base font size for plotting.
#' @param export A logical indicating whether to export the plot as PDF and PNG file into the /Plots folder.
#' @return (Invisibly) A ggplot object.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom assertthat assert_that
#' @importFrom purrr map_df
#' @importFrom tidyr gather
#'
rna.plot_imputation <- function(assay, ..., colData, plot = TRUE, basesize = 12, export = TRUE)
{
  call <- match.call()
  call$colData <- NULL
  call$export <- NULL
  call$plot <- NULL
  call$basesize <- NULL
  arglist <- lapply(call[-1], function(x) x)
  var.names <- vapply(arglist, deparse, character(1))
  arglist <- lapply(arglist, eval.parent, n = 2)
  names(arglist) <- names(var.names)
  lapply(arglist, function(x) {
    assertthat::assert_that(inherits(x, "matrix"),
                            msg = "input objects need to be of class 'matrix' or 'array'")
  })
  gather_join <- function(counts) {
    counts %>% data.frame(check.names = FALSE) %>% gather(ID, val) %>% left_join(., data.frame(colData, check.names = FALSE), by = "ID")
  }
  df <- purrr::map_df(arglist, gather_join, .id = "var") %>%
    mutate(var = factor(var, levels = names(arglist)))
  p <- ggplot(df, aes(val, col = condition)) +
    geom_density(na.rm = TRUE, show.legend = F) +
    stat_density(geom="line", position ="identity", na.rm = TRUE) +
    facet_wrap(~var, ncol = 1) +
    labs(x = expression(log[2] ~ "Intensity"), y = "Density") +
    theme_DEP1(basesize = basesize)

  if (length(unique(df$condition)) <= 8) {
    p <- p + scale_color_brewer(palette = "Dark2")
  } else if (length(unique(se$condition)) <= 12) {
    p <- p + scale_color_brewer(palette = "Set3")
  } else {
    pal <- c(
      "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
      "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
      "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
      "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
      "green1", "yellow4", "yellow3", "darkorange4", "brown")
    p <- p + scale_color_manual(values = pal[1:length(unique(se$condition))])
  }

  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
    message(paste0("Exporting imputation plot to:\n\"", getwd(), "/Plots/ImputationPlot.pdf\" and \".../ImputationPlot.png\""))
    grDevices::png("Plots/ImputationPlot.png",
                   width = 7, height = 8, units = 'in', res = 300)
    print(p)
    grDevices::dev.off()
    grDevices::pdf("Plots/ImputationPlot.pdf")
    print(p)
    grDevices::dev.off()
  }

  if (plot == TRUE){
    print(p)
  }
    invisible(p)

}


#' Plots a volcano plot for a given contrast
#'
#' @param dds A DESeqDataSet object
#' @param contrast The contrast of interest in the format "conditionA_vs_conditionB"
#' @param label_size The size of labels for the gene names
#' @param alpha The alpha value used to determine significance.
#' @param lfc The log2 fold change threshold.
#' @param add_names Logical, whether to add gene names to the plot
#' @param adjusted Logical, whether to plot adjusted p-values. The \code{alpha} threshold is applied to adjusted p values even if \code{adjusted = FALSE}.
#' @param plot Logical, whether to return the volcano plot
#' @param export Logical, whether to export the volcano plot as PND and PDF file
#'
#' @return (invisibly) A data frame with columns for gene, log2 fold change, and -log10 p-value
#'
#' @export
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom assertthat assert_that
#' @importFrom ggplot2 geom_point geom_vline geom_hline labs xlab ylab scale_colour_manual theme_bw theme coord_cartesian
#' @importFrom ggh4x force_panelsizes
#' @importFrom scales pretty_breaks
#' @importFrom stringr str_replace
#' @importFrom grDevices pdf png dev.off
rna.plot_volcano <-
  function (dds,
            contrast,
            label_size = 3,
            alpha = 0.05,
            lfc = 2,
            x.lim = NULL,
            add_names = TRUE,
            adjusted = FALSE,
            plot = TRUE,
            export = TRUE)
  {
    if (is.integer(label_size))
      label_size <- as.numeric(label_size)
    assertthat::assert_that(
      inherits(dds, "SummarizedExperiment"),
      is.character(contrast),
      length(contrast) == 1,
      is.numeric(label_size),
      length(label_size) == 1,
      is.logical(add_names),
      length(add_names) ==
        1,
      is.logical(adjusted),
      length(adjusted) == 1,
      is.logical(plot),
      length(plot) == 1
    )
    lfc.pfx <- ifelse(length(grep("lfc_shrink", colnames(SummarizedExperiment::rowData(dds))))>0,
                      "lfc_shrink\\.", "lfc\\.")
    row_data <- SummarizedExperiment::rowData(dds, use.names = FALSE)
    if (any(!c("name", "ID") %in% colnames(row_data))) {
      stop(
        paste0(
          "'name' and/or 'ID' columns are not present in '",
          deparse(substitute(dds)),
          "'.\nRun make_unique() to obtain required columns."
        ),
        call. = FALSE
      )
    }
    if (length(grep(paste0("^padj\\.|",lfc.pfx), colnames(row_data))) < 1) {
      stop(
        paste0(
          "'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
          deparse(substitute(dds)),
          "'.\nRun test_diff() to obtain the required columns."
        ),
        call. = FALSE
      )
    }
    if (length(grep("^significant\\.", colnames(row_data))) <
        1) {
      stop(
        paste0(
          "'[contrast]_significant' columns are not present in '",
          deparse(substitute(dds)),
          "'.\nRun add_rejections() to obtain the required columns."
        ),
        call. = FALSE
      )
    }
    if (length(grep(paste(lfc.pfx, contrast, sep = ""),
                    colnames(row_data))) == 0) {
      valid_cntrsts <-
        row_data %>% data.frame(check.names = FALSE) %>% select(starts_with(gsub("\\\\.",".", lfc.pfx))) %>%
        colnames(.) %>% gsub(lfc.pfx, "", .)
      valid_cntrsts_msg <- paste0("Valid contrasts are: '",
                                  paste0(valid_cntrsts, collapse = "', '"),
                                  "'")
      stop(
        "Not a valid contrast, please run `rna.plot_volcano()` with a valid contrast as argument\n",
        valid_cntrsts_msg,
        call. = FALSE
      )
    }
    diff <- match(paste(gsub("\\\\.",".", lfc.pfx), contrast, sep = ""),
                  colnames(row_data))
    if (adjusted) {
      p_values <- match(paste("padj.", contrast, sep = ""),
                        colnames(row_data))
    }
    else {
      p_values <- match(paste("pvalue.", contrast, sep = ""),
                        colnames(row_data))
    }
    signif <- match(paste("significant.", contrast, sep = ""),
                    colnames(row_data))

    # replace NA p values with the lowest p value
    na.ndx <- which(is.na(row_data[, p_values]))
    row_data[is.na(row_data[, p_values]), p_values] <- max(row_data[, p_values][!is.na(row_data[, p_values])])

    # replace 0 p values with the machine-specific lowest value
    row_data[,signif][is.na(row_data[,signif])] <- FALSE

    df <- data.frame(
      x = row_data[, diff],
      y = -log10(row_data[, p_values]),
      significant = row_data[, signif],
      name = row_data$name,
      na = seq(1:length(row_data[, signif])) %in% na.ndx,
      check.names = FALSE
    ) %>%
      filter(!is.na(significant)) %>% arrange(significant)

    quant99.5 <-  quantile(df$y[(df$y!=0)&!is.infinite(df$y)], probs = 0.995, na.rm = TRUE)

    # create new columns for shapes "circle" and "triangle" if genes lie within or outside (for genes with very small adj.p value < 3.162278e-05) of the y bounds, respectively.
    if (adjusted) {
      if (quant99.5 > 1.5*quantile(df$y[(df$y!=0)&!is.infinite(df$y)], probs = 0.985, na.rm = TRUE)) {
        df$shape <- ifelse(df$y > 0.98 * quant99.5, "triangle", "circle")
      }
      else{
        df$shape <- "circle"
      }
    }
    else{
      df$shape <- "circle"
    }
    # change the -Log10(q value) of genes exceeding the y plot bound so that they are displayed at the axis border
    if (adjusted) {
      if (quant99.5 > 1.5*quantile(df$y[(df$y!=0)&!is.infinite(df$y)], probs = 0.985, na.rm = TRUE)) {
        df$y[df$y > 0.98 * quant99.5] <- 0.98 * quant99.5
      }
    }
    df_volcano <- df %>%
      # Categorize genes based on the log2(fold change) and q value thresholds
      mutate(Genes = factor(
        case_when(
          x >= lfc & significant == "TRUE" ~ "Upregulated",
          x <= -lfc & significant == "TRUE" ~ "Downregulated",
          na == "TRUE" ~ "p-value NA",
          TRUE ~ "Not significant"
        )
      ))
    df_volcano$Genes <-relevel(df_volcano$Genes, "Upregulated")
    # levels(df_volcano$Genes) <- rev(levels(df_volcano$Genes))
    # Plot with points colored according to the thresholds
    plot_volcano <- df_volcano %>%
      ggplot(aes(x, y, group = name, colour = Genes)) +
      geom_point(aes(shape = shape), size = 2.5, alpha = 0.6) +   # Define circle shape, size and transparency
      scale_shape(guide = "none") +   # Hide shapes from legend
      geom_vline(xintercept = lfc,
                 linetype = "dashed",
                 alpha = 0.5) +       # Add dotted lines to indicate the log2(fc) threshold, semi-transparent
      geom_vline(
        xintercept = -lfc,
        linetype = "dashed",
        alpha = 0.5 ) +      # Add dotted lines to indicate the log2(fc) threshold, semi-transparent
      scale_colour_manual(
        # Adjust the circle colors
        labels = c("Upregulated", "Downregulated", "Not significant", "p-value NA"),
        values = c(
          "Upregulated" = "#e31f26",
          "Not significant" = "#6f716f",
          "Downregulated" = "#387fb9",
          "p-value NA" = "#353535")) +
      labs(title =  str_replace(contrast, "_vs_", " vs. ")) +
      xlab(expression(log[2]("Fold change"))) + ylab(expression(-Log[10]("q-value"))) +       # Label the axes
      theme_bw(base_size = 20) + #Set the theme, define label font size
      theme(legend.position = "bottom", title = element_text(size = exp(-nchar(str_replace(contrast, "_vs_", " vs. "))/40)*28),
            axis.title = element_text(size = 22),
            legend.title = element_text(size = 18)) +
      ggh4x::force_panelsizes(rows = unit(6.5, "in"),
                              cols = unit(6.5, "in")) +
      guides(color = guide_legend(reverse = FALSE))

    if(!is.null(x.lim))
      plot_volcano <- plot_volcano + scale_x_continuous(breaks = scales::pretty_breaks(n = max(abs(df$x)) + 1), limits = x.lim)
    else
      plot_volcano <- plot_volcano + scale_x_continuous(breaks = scales::pretty_breaks(n = max(abs(df$x)) + 1))


    if (adjusted) {
      if(is.null(x.lim)){
        plot_volcano <- plot_volcano + coord_cartesian(
          xlim = c(-max(abs(df$x)) - 0.3,  # Lower x bound [defined dynamically based on the largest absolute log2(fc) value]
                   max(abs(df$x)) + 0.3),
          # Upper x bound [defined dynamically based on the largest absolute log2(fc) value]
          ylim = c((min(df$y) - 0.1), ifelse(quant99.5 > 1.5*quantile(df$y[(df$y!=0)&!is.infinite(df$y)], probs = 0.985, na.rm = TRUE), quant99.5, max(df$y[(df$y!=0)&!is.infinite(df$y)])*1.02)),
          expand = FALSE
        )     # Plot y bounds (the bottom axis is defined dynamically based on the smallest q value)
      } else {
        plot_volcano <- plot_volcano + coord_cartesian(
          # Upper x bound [defined dynamically based on the largest absolute log2(fc) value]
          ylim = c((min(df$y) - 0.1), ifelse(quant99.5 > 1.5*quantile(df$y[(df$y!=0)&!is.infinite(df$y)], probs = 0.985, na.rm = TRUE), quant99.5, max(df$y[(df$y!=0)&!is.infinite(df$y)])*1.02)),
          expand = FALSE
        )     # Plot y bounds (the bottom axis is defined dynamically based on the smallest q value)
      }

    }
    else{
      if(is.null(x.lim)){
        plot_volcano <- plot_volcano + coord_cartesian(
          xlim = c(-max(abs(df$x)) - 0.3,  # Lower x bound [defined dynamically based on the largest absolute log2(fc) value]
                   max(abs(df$x)) + 0.3),
          # Upper x bound [defined dynamically based on the largest absolute log2(fc) value]
          ylim = c((min(df$y) - 0.1), max(df$y[!is.infinite(df$y)])),
          expand = TRUE
        )     # Plot y bounds (the bottom axis is defined dynamically based on the smallest q value)
      } else {
        plot_volcano <- plot_volcano + coord_cartesian(
          ylim = c((min(df$y) - 0.1), max(df$y[!is.infinite(df$y)])),
          expand = TRUE
        )     # Plot y bounds (the bottom axis is defined dynamically based on the smallest q value)
      }

    }
    if (adjusted) {
      plot_volcano <- plot_volcano +
        geom_hline(yintercept = -log10(alpha),
                   linetype = "dashed",
                   alpha = 0.5)      # Add dotted line to indicate the adj.p value threshold, semi-transparent
    }
    if (add_names) {
      plot_volcano <- plot_volcano + ggrepel::geom_text_repel(
        data = filter(df_volcano, significant|na),
        aes(label = name),
        size = label_size,
        box.padding = unit(0.25,
                           "lines"),
        point.padding = unit(0.1, "lines"),
        segment.size = 0.5,
        show.legend = FALSE,
        max.overlaps = 15
      )
    }
    if (adjusted) {
      plot_volcano <-
        plot_volcano + labs(y = expression(-log[10] ~ "(adj. p-value)"))
    }
    else {
      plot_volcano <-
        plot_volcano + labs(y = expression(-log[10] ~ "(p-value)"))
    }
    if (export == TRUE) {
      dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
      message(
        paste0(
          "Exporting volcano plot to:\n\"",
          getwd(),
          "/Plots/VolcanoPlot_", contrast,
          ".pdf\" and \".../VolcanoPlot_", contrast, ".png\""
        )
      )

      grDevices::pdf(paste0("Plots/VolcanoPlot_", contrast, ".pdf"),
                     width = 8.1,
                     height = 8.6)
      print(plot_volcano)
      grDevices::dev.off()

      grDevices::png(
        paste0("Plots/VolcanoPlot_", contrast, ".png"),
        width = 8.1,
        height = 8.6,
        units = 'in',
        res = 300
      )
      print(plot_volcano)
      grDevices::dev.off()
    }
    if (plot) {
      print(plot_volcano)

    }
    else {
      df <- df %>% select(name, x, y, significant) %>% arrange(desc(x))
      colnames(df)[c(1, 2, 3)] <- c("gene", "log2_fold_change",
                                    "p_value_-log10")
      if (adjusted) {
        colnames(df)[3] <- "adjusted_p_value_-log10"
      }
      return(df)
    }
  }

#' Plot the number of genes identified in a DESeqDataSet object
#'
#' @param se A DESeqDataSet object
#' @param plot Logical, if TRUE a plot is returned
#' @param export Logical, if TRUE the plot is exported to the /Plots directory
#'
#' @return If plot = TRUE, returns a ggplot object. If plot = FALSE, returns a data frame with plotted values.
#'
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom SummarizedExperiment colData assay
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 geom_col geom_hline labs scale_fill_brewer scale_fill_manual
#' @importFrom grDevices png pdf
#'
rna.plot_numbers <- function (se, plot = TRUE, export = FALSE, basesize = 12)
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.logical(plot), length(plot) == 1)
  df <- SummarizedExperiment::assay(se) %>% data.frame() %>% tibble::rownames_to_column() %>%
    gather(ID, bin, -rowname) %>% mutate(bin = ifelse(is.na(bin),
                                                      0, 1))
  stat <- df %>% group_by(ID) %>% dplyr::summarize(n = n(), sum = sum(bin)) %>%
    left_join(., data.frame(SummarizedExperiment::colData(se)), by = "ID")
  p <- ggplot(stat, aes(x = ID, y = sum, fill = condition)) +
    geom_col() + geom_hline(yintercept = unique(stat$n)) +
    labs(title = "Genes per sample", x = "",
         y = "Number of genes") + theme_DEP2()

  if (length(unique(se$condition)) <= 8) {
    p <- p + scale_fill_brewer(palette = "Dark2")
  } else if (length(unique(se$condition)) <= 12) {
    p <- p + scale_fill_brewer(palette = "Set3")
  } else {
    p <- p + scale_fill_manual(
      values = c(
        "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
        "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
        "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
        "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
        "green1", "yellow4", "yellow3", "darkorange4", "brown"
      )
    )
  }
  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
    fig_width = ncol(se)/1.5
    message(paste0("Exporting plot with identified genes per sample to:\n\"", getwd(), "/Plots/GenesIdentified.pdf\" and \".../GenesIdentified.png\""))
    grDevices::png("Plots/GenesIdentified.png",
                   width = fig_width, height = 8, units = 'in', res = 300)
    print(p)
    grDevices::dev.off()
    grDevices::pdf("Plots/GenesIdentified.pdf",
                   width = fig_width, height = 8)
    print(p)
    grDevices::dev.off()
  }
  if (plot) {
    return(p)
  }
  else {
    df <- as.data.frame(stat)
    colnames(df)[seq_len(3)] <- c("sample", "total_genes",
                                  "genes_in_sample")
    return(df)
  }
}

#' Plot bar plots for gene abundance or fold change
#'
#' This function generates bar plots for gene abundance or fold change, normalized by either a reference gene or centered to the mean. The user can also plot the contrast of two conditions.
#'
#' @param dds A DESeqDataSet object containing the data.
#' @param genes A vector of gene names or IDs.
#' @param type Type of plot: "centered" or "contrast".
#' @param shrunken.lfc A logical indicating if shrunken log2 fold changes should be used if \code{type = "contrast"}. Default is TRUE.
#' @param contrast Name of the contrast for plotting when type is "contrast".
#' @param col.id Name of the column with gene IDs in \code{dep}.
#' @param name_table A data frame containing the mapping of gene names to reference names.
#' @param match.id The name of a column in \code{name_table} to find matches with \code{col.id}.
#' @param match.name The name of a column in \code{name_table} to assign new names based on matches of \code{col.id} and \code{match.id}.
#' @param convert_name Logical value indicating whether to convert the gene names to reference names based on \code{name_table} (e.g., gene names).
#' @param shape.size Numeric value defining the (replicate) symbol size in the plot
#' @param collage A logical indicating if multiple plots for different genes should be combined into collages of up to 8 plots. Default is TRUE.
#' @param plot wether to return the plot or not.
#' @param export Logical; whether to export the plot as PDF an PNG files
#' @param export.nm Name of the output PDF and PNG files if \code{export = TRUE}.
#'
#' @return If `plot` is set to FALSE and the length of `genes` is less than or equal to 8, a data frame with the log2 fold changes or log2 intensities, along with their lower and upper bounds of 95% confidence intervals, is returned.
#'
#' @export
#'
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom tidyr gather spread
#' @importFrom stats qnorm
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stringr str_split
#'
rna.plot_bar <- function (dds,
                          genes,
                          type = c("contrast", "centered"),
                          shrunken.lfc = T,
                          contrast = NULL,
                          col.id = "ID",
                          name_table = NULL,
                          match.id = "Accession",
                          match.name = "Name",
                          convert_name = FALSE,
                          shape.size = 2.5,
                          collage = TRUE,
                          plot = TRUE,
                          export = TRUE,
                          export.nm = NULL)
{
  assertthat::assert_that(inherits(dds, "SummarizedExperiment"),
                          is.character(genes), is.character(type), is.logical(plot),
                          length(plot) == 1)
  type <- match.arg(type)
  if(shrunken.lfc == TRUE && length(grep("lfc_shrink", colnames(SummarizedExperiment::rowData(dds))))>0){
    lfc.pfx <- "lfc_shrink."
  } else {
    lfc.pfx <- "lfc."
  }
  # Replace rownames of dds with column of original data frame defined in "col.id" argument
  row_data <- SummarizedExperiment::rowData(dds, use.names = F)
  GeneID <- match(col.id,
                  colnames(row_data))
  rownames(dds) <- row_data[, GeneID]
  rownames(dds)[is.na(rownames(dds))] <- 0
  row_data <- SummarizedExperiment::rowData(dds, use.names = FALSE)

  if (any(!c("label", "condition", "replicate") %in%
          colnames(SummarizedExperiment::colData(dds)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(dds)), "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  if (length(grep("padj\\.|diff\\.", colnames(row_data))) < 1 && type == "contrast") {
    stop("'diff.[contrast]' and 'padj.[contrast]' columns are not present in '",
         deparse(substitute(dds)), "'\nRun rna.workflow() to obtain the required columns",
         call. = FALSE)
  }
  if (!"name" %in% colnames(row_data)) {
    stop("'name' column not present in '", deparse(substitute(dds)),
         "'\nRun rna.read_data() to obtain the required columns",
         call. = FALSE)
  }
  if (all(!genes %in% row_data[,col.id])) {
    if (length(genes) == 1) {
      rows <- grep(substr(genes, 1, nchar(genes) -
                            1), row_data[,col.id])
      possibilities <- row_data[,col.id][rows]
    }
    else {
      rows <- lapply(genes, function(x) grep(substr(x,
                                                       1, nchar(x) - 1), row_data[,col.id]))
      possibilities <- row_data[,col.id][unlist(rows)]
    }
    if (length(possibilities) > 0) {
      possibilities_msg <- paste0("Do you mean: '",
                                  paste0(possibilities, collapse = "', '"),
                                  "'")
    }
    else {
      possibilities_msg <- NULL
    }
    stop(paste0("The genes ", paste(genes, collapse=", "), " were not found in the column '", col.id,
                "'. Please run `prot.plot_bar()` with a valid gene names in the 'genes' argument or provide the valid column ID in the 'col.id' argument.\n"),
         possibilities_msg, call. = FALSE)
  }
  if (any(!genes %in% row_data[,col.id])) {
    genes <- genes[genes %in% row_data[,col.id]]
    warning("Only used the following gene(s): '",
            paste0(genes, collapse = "', '"), "'", call. = F)
  }
  subset_list <- list()
  if(collage){
    if(length(genes)<=8){
      subset_list[[1]] <- dds[genes]
    } else if(length(genes)<=16){
      subset_list[[1]] <- dds[genes[1:8]]
      subset_list[[2]] <- dds[genes[9:length(genes)]]
    } else if(length(genes)<=24){
      subset_list[[1]] <- dds[genes[1:8]]
      subset_list[[2]] <- dds[genes[9:16]]
      subset_list[[3]] <- dds[genes[17:length(genes)]]
    }else if(length(genes)<=32){
      subset_list[[1]] <- dds[genes[1:8]]
      subset_list[[2]] <- dds[genes[9:16]]
      subset_list[[3]] <- dds[genes[17:24]]
      subset_list[[4]] <- dds[genes[25:length(genes)]]
    } else if(length(genes)<=40){
      subset_list[[1]] <- dds[genes[1:8]]
      subset_list[[2]] <- dds[genes[9:16]]
      subset_list[[3]] <- dds[genes[17:24]]
      subset_list[[4]] <- dds[genes[25:32]]
      subset_list[[5]] <- dds[genes[33:length(genes)]]
    } else if(length(genes)<=48){
      subset_list[[1]] <- dds[genes[1:8]]
      subset_list[[2]] <- dds[genes[9:16]]
      subset_list[[3]] <- dds[genes[17:24]]
      subset_list[[4]] <- dds[genes[25:32]]
      subset_list[[5]] <- dds[genes[33:40]]
      subset_list[[6]] <- dds[genes[41:length(genes)]]
    } else if(length(genes)<=56){
      subset_list[[1]] <- dds[genes[1:8]]
      subset_list[[2]] <- dds[genes[9:16]]
      subset_list[[3]] <- dds[genes[17:24]]
      subset_list[[4]] <- dds[genes[25:32]]
      subset_list[[5]] <- dds[genes[33:40]]
      subset_list[[6]] <- dds[genes[41:48]]
      subset_list[[7]] <- dds[genes[49:length(genes)]]
    } else {
      subset_list[[1]] <- dds[genes[1:8]]
      subset_list[[2]] <- dds[genes[9:16]]
      subset_list[[3]] <- dds[genes[17:24]]
      subset_list[[4]] <- dds[genes[25:32]]
      subset_list[[5]] <- dds[genes[33:40]]
      subset_list[[6]] <- dds[genes[41:48]]
      subset_list[[7]] <- dds[genes[49:56]]
      subset_list[[8]] <- dds[genes[57:length(genes)]]
    }
  } else {
    subset_list <- lapply(1:length(genes), function(x) dds[genes[x]])
  }

  for(i in 1:length(subset_list)) {
    if (type == "centered") {
      means <- rowMeans(log2(SummarizedExperiment::assay(subset_list[[i]])), na.rm = TRUE)
      df_reps <-
        data.frame(log2(SummarizedExperiment::assay(subset_list[[i]])) - means, check.names = FALSE) %>% tibble::rownames_to_column() %>%
        gather(ID, val, -rowname) %>% left_join(., data.frame(SummarizedExperiment::colData(subset_list[[i]]), check.names = FALSE),
                                                by = "ID")
      if (convert_name == TRUE) {
        df_reps$rowname <- transform(df_reps,
                                     rowname = name_table[match(paste(df_reps$rowname),
                                                                paste(unlist(str_split(
                                                                  Reduce(c, name_table[match.id]), ", "
                                                                )))),
                                                          match.name]) %>%
          select(rowname) %>% unlist(., use.names = F)
      }
      df_reps$replicate <- as.factor(df_reps$replicate)
      df <-
        df_reps %>% group_by(condition, rowname) %>% dplyr::summarize(
          mean = mean(val,
                      na.rm = TRUE),
          sd = sd(val, na.rm = TRUE),
          n = n()
        ) %>%
        mutate(
          error = stats::qnorm(0.975) * sd / sqrt(n),
          CI.L = mean -
            error,
          CI.R = mean + error
        ) %>% data.frame(check.names = FALSE)



      p <-
        ggplot(df, aes(condition, mean)) + geom_hline(yintercept = 0) +
        geom_col(aes(y = mean, fill = condition), colour = "black")
      if (length(unique(dds$condition)) <= 8) {
        p <- p + scale_fill_manual(values = RColorBrewer::brewer.pal(n = length(unique(dds$condition)), name = "Dark2"))
      } else if (length(unique(dds$condition)) <= 12) {
        p <- p + scale_fill_manual(values = RColorBrewer::brewer.pal(n = length(unique(dds$condition)), name = "Set3"))
      } else {
        p <- p + scale_fill_manual(values = c(
          "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
          "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
          "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
          "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
          "green1", "yellow4", "yellow3", "darkorange4", "brown"
        ))
      }
      p <- p +
        ggnewscale::new_scale_fill() +
        geom_point(
          data = df_reps,
          aes(condition, val, fill = replicate),
          shape = 23,
          size = shape.size,
          color = "black",
          position = position_dodge(width = 0.3)
        ) +
        scale_fill_manual(values = case_when(
          as.numeric(max(dds$replicate))<=8       ~ RColorBrewer::brewer.pal(n=as.numeric(max(dds$replicate)), name="Greys"),
          as.numeric(max(dds$replicate))>8        ~ grDevices::colorRampPalette(RColorBrewer::brewer.pal(n=8, name="Greys"))(as.numeric(max(dds$replicate))))
        ) +
        geom_errorbar(aes(ymin = CI.L, ymax = CI.R), width = 0.3) +
        labs(
          x = "Baits",
          y = expression(log[2] ~ "Centered intensity" ~
                           "(\u00B195% CI)"),
          col = "Rep"
        ) + facet_wrap( ~ rowname) +
        theme(basesize = 12) + theme_DEP2()
      w <- 8+log(max(str_count(df[,"condition"]))-3, base = 1.6)
      h <- 10
    }
    if (type == "contrast") {
      if(is.null(contrast)){
        contrast <- SummarizedExperiment::rowData(subset_list[[i]], use.names = FALSE) %>% data.frame(check.names = FALSE) %>% select(starts_with("lfc."))  %>%
          colnames(.) %>% gsub("lfc.", "", .)
      }
      if (convert_name == TRUE) {
          df <-SummarizedExperiment::rowData(subset_list[[i]], use.names = FALSE) %>% data.frame(check.names = FALSE)
          df <- df %>%
            select(
              name,
              paste0(lfc.pfx, contrast),
              paste0("lfcSE", gsub("lfc", "", lfc.pfx), contrast)) %>%
            gather(var, val, -name) %>%
            mutate(
              contrast = gsub(paste(lfc.pfx, paste0(paste0("lfcSE", gsub("lfc", "", lfc.pfx))), sep = "|"), "", var),
              var = gsub("\\..+", "", var))  %>%
            spread(var, val)
          error = stats::qnorm(0.975) * df[[gsub("\\.", "", paste0("lfcSE", gsub("lfc", "", lfc.pfx)))]]
          df$CI_L <- df[[gsub("\\.", "", lfc.pfx)]] - error
          df$CI_R <- df[[gsub("\\.", "", lfc.pfx)]] + error

        p <-
          ggplot(df, aes(contrast, .data[[gsub("\\.", "", lfc.pfx)]])) + geom_hline(yintercept = 0) +
          geom_col(aes(y = .data[[gsub("\\.", "", lfc.pfx)]], fill = contrast), colour = "black")
        if (length(unique(dds$condition)) <= 8) {
          p <- p + scale_fill_manual(values = RColorBrewer::brewer.pal(n = length(unique(dds$condition)), name = "Dark2"))
        } else if (length(unique(dds$condition)) <= 12) {
          p <- p + scale_fill_manual(values = RColorBrewer::brewer.pal(n = length(unique(dds$condition)), name = "Set3"))
        } else {
          p <- p + scale_fill_manual(values = c(
            "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
            "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
            "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
            "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
            "green1", "yellow4", "yellow3", "darkorange4", "brown"
          ))
        }
        p <- p +
          geom_errorbar(aes(ymin = CI_L, ymax = CI_R), width = 0.3) +
          labs(
            x = "Contrasts",
            y = expression(log[2] ~ "Fold change" ~
                             "(\u00B195% CI)")) + facet_wrap( ~ name) +
          theme(basesize = 12) + theme_DEP2()


      } else {
          df <-SummarizedExperiment::rowData(subset_list[[i]], use.names = FALSE) %>% data.frame(check.names = FALSE)
          df <- df %>%
            select(
              name,
              paste0(lfc.pfx, contrast),
              paste0("lfcSE", gsub("lfc", "", lfc.pfx), contrast)) %>%
            gather(var, val, -name) %>%
            mutate(
              contrast = gsub(paste(lfc.pfx, paste0(paste0("lfcSE", gsub("lfc", "", lfc.pfx))), sep = "|"), "", var),
              var = gsub("\\..+", "", var))  %>%
            spread(var, val)
          error = stats::qnorm(0.975) * df[[gsub("\\.", "", paste0("lfcSE", gsub("lfc", "", lfc.pfx)))]]
          df$CI_L <- df[[gsub("\\.", "", lfc.pfx)]] - error
          df$CI_R <- df[[gsub("\\.", "", lfc.pfx)]] + error

        p <-
          ggplot(df, aes(contrast, .data[[gsub("\\.", "", lfc.pfx)]])) + geom_hline(yintercept = 0) +
          geom_col(aes(y = .data[[gsub("\\.", "", lfc.pfx)]], fill = contrast), colour = "black")
        if (length(unique(dds$condition)) <= 8) {
          p <- p + scale_fill_manual(values = RColorBrewer::brewer.pal(n = length(unique(dds$condition)), name = "Dark2"))
        } else if (length(unique(dds$condition)) <= 12) {
          p <- p + scale_fill_manual(values = RColorBrewer::brewer.pal(n = length(unique(dds$condition)), name = "Set3"))
        } else {
          p <- p + scale_fill_manual(values = c(
            "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
            "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
            "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
            "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
            "green1", "yellow4", "yellow3", "darkorange4", "brown"
          ))
        }
        p <- p +
          geom_errorbar(aes(ymin = CI_L, ymax = CI_R), width = 0.3) +
          labs(
            x = "Contrasts",
            y = expression(log[2] ~ "Fold change" ~
                             "(\u00B195% CI)")) + facet_wrap( ~ name) +
          theme(basesize = 12) + theme_DEP2()
      }
      w <- 8+log(max(str_count(df[,"contrast"]))-3, base = 1.6)
      h <- 10
    }

    if (export == TRUE) {
      if(collage){
        if (!is.null(export.nm)) {
          if (length(genes) <= 8) {
            export_name <- paste0("Plots/", export.nm)
          } else{
            export_name <- paste0("Plots/", export.nm, "_", i)
          }
        } else {
          if (length(genes) <= 8) {
            export_name <- paste0("Plots/BarPlot_",
                                  paste(genes, collapse = "_"))
          } else{
            export_name <- paste0("Plots/BarPlot_",
                                  paste(genes, collapse = "_"), "_", i)
          }
        }
      } else {
        export_name <- paste0("Plots/BarPlot_",
                              paste(genes[i]))
      }
      dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
      message(
        paste0(
          "Exporting bar plot(s) to:\n\"",
          getwd(),
          export_name,
          ".pdf",
          " and \"...",
          export_name,
          ".png\""
        )
      )

      grDevices::pdf(paste0(export_name, ".pdf"),
                     width = w,
                     height = h)
      print(p)
      grDevices::dev.off()

      grDevices::png(
        paste0(export_name, ".png"),
        width = w,
        height = h,
        units = 'in',
        res = 300
      )
      print(p)
      grDevices::dev.off()
    }
    if(length(genes)>8){
      plot = FALSE
    }
    if(plot){
      print(p)
    } else {
      if (type == "centered") {
        df <- df %>% select(rowname, condition, mean, CI.L,
                            CI.R)
        colnames(df) <- c("gene", "condition",
                          "log2_intensity", "CI.L", "CI.R")
      }
      if (type == "contrast") {
        df <- df %>% select(name, contrast, diff, CI.L, CI.R) %>%
          mutate(contrast = paste0(contrast))
        colnames(df) <- c("gene", "contrast",
                          "log2_fold_change", "CI.L", "CI.R")
      }
      if (length(genes) <= 8) {
        return(df)
      }
    }
  }
}


