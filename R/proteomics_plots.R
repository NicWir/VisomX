#' @title  theme_DEP1
#' @description A function to create a theme for ggplot2
#' @param basesize numeric. Base font size for the theme.
#' @return A ggplot2 theme object.
#' @export
#' @importFrom ggplot2 theme_bw
theme_DEP1 <- function (basesize = 15)
{
  theme <- ggplot2::theme_bw(base_size = basesize)
  theme$plot.title$face <- "bold"
  theme$plot.title$size <- basesize + 2
  theme$plot.title$hjust <- 0.5
  theme$axis.title.x$size <- basesize + 2
  theme$axis.title.y$size <- basesize + 2
  theme$axis.text$size <- basesize
  theme$axis.text$colour <- "black"
  theme$legend.title$size <- basesize + 2
  theme$legend.text$size <- basesize
  theme$strip.text$face <- "bold"
  theme$strip.text$size <- basesize + 2
  theme$strip.text$colour <- "black"
  return(theme)
}



#' Plot the number of proteins identified in a SummarizedExperiment object
#'
#' @param se A SummarizedExperiment object
#' @param basesize Numeric. Base font size for the plot
#' @param plot Logical, if TRUE a plot is returned
#' @param export Logical, if TRUE the plot is exported to the /Plots directory
#'
#' @return If plot = TRUE, returns a ggplot object. If plot = FALSE, returns a data frame.
#'
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom SummarizedExperiment colData assay
#' @importFrom ggplot2 ggplot aes geom_col geom_hline labs scale_fill_brewer scale_fill_manual
#' @importFrom grDevices png pdf
#'
prot.plot_numbers <- function (se, basesize = 12, plot = TRUE, export = FALSE)
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.logical(plot), length(plot) == 1)
  df <- SummarizedExperiment::assay(se) %>% data.frame() %>% tibble::rownames_to_column() %>%
    gather(ID, bin, -rowname) %>% mutate(bin = ifelse(is.na(bin),
                                                      0, 1))
  stat <- df %>% group_by(ID) %>% summarize(n = n(), sum = sum(bin)) %>%
    left_join(., data.frame(SummarizedExperiment::colData(se)), by = "ID")
  p <- ggplot(stat, aes(x = ID, y = sum, fill = condition)) +
    geom_col() +
    geom_hline(yintercept = unique(stat$n)) +
    labs(title = "Proteins per sample", x = "",
         y = "Number of proteins") +
    theme_DEP2(basesize=basesize)

  if (length(unique(se$condition)) <= 8) {
    p <- p + scale_fill_brewer(palette = "Dark2")
  } else if (length(unique(se$condition)) <= 12) {
    p <- p + scale_fill_brewer(palette = "Set3")
  } else if (length(unique(se$condition)) <= 25) {
    p <- p + scale_fill_manual(
      values = c(
        "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
        "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
        "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
        "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
        "green1", "yellow4", "yellow3", "darkorange4", "brown"
      )
    )
  } else {
    p <- p + scale_fill_manual(values = rainbow(length(unique(se$condition))))
  }

  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
    fig_width = ncol(se)/1.5
    message(paste0("Exporting plot with identified proteins per sample to:\n\"", getwd(), "/Plots/ProteinsIdentified.pdf\" and \".../ProteinsIdentified.png\""))
    grDevices::png("Plots/ProteinsIdentified.png",
                   width = fig_width, height = 8, units = 'in', res = 300)
    print(p)
    grDevices::dev.off()
    grDevices::pdf("Plots/ProteinsIdentified.pdf",
                   width = fig_width, height = 8)
    print(p)
    grDevices::dev.off()
  }
  if (plot) {
    return(p)
  }
  else {
    df <- as.data.frame(stat)
    colnames(df)[seq_len(3)] <- c("sample", "total_proteins",
                                  "proteins_in_sample")
    return(df)
  }
}



#' @title Plot imputation results
#' @description Plot the results of imputation as a density plot. It is recommended to log-transform the data before using this function.
#' @param se SummarizedExperiment object
#' @param ... Additional SummarizedExperiment objects. Can be input as named argument, for example: \code{"se2" = se2}
#' @param plot logical, whether to plot the results or not (default: TRUE)
#' @param basesize numeric, base font size (default: 12)
#' @param export logical, whether to export the results as PNG and PDF or not (default: TRUE)
#' @return (invisible) a ggplot object
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom purrr map_df
#' @importFrom ggplot2 ggplot aes geom_density stat_density facet_wrap labs scale_color_brewer scale_color_manual
#' @importFrom grDevices png pdf dev.off
#' @importFrom SummarizedExperiment assay colData
#' @importFrom magrittr %>%
prot.plot_imputation <- function(se, ..., plot = TRUE, basesize = 12, export = FALSE)
{
  call <- match.call()
  call$export <- NULL
  call$plot <- NULL
  call$basesize <- NULL
  arglist <- lapply(call[-1], function(x) x)
  var.names <- vapply(arglist, deparse, character(1))
  arglist <- lapply(arglist, eval.parent, n = 2)
  names(arglist) <- var.names
  lapply(arglist, function(x) {
    assertthat::assert_that(inherits(x, "SummarizedExperiment"),
                            msg = "input objects need to be of class 'SummarizedExperiment'")
    if (any(!c("label", "ID", "condition",
               "replicate") %in% colnames(SummarizedExperiment::colData(x)))) {
      stop("'label', 'ID', 'condition' and/or 'replicate' ",
           "columns are not present in (one of) the input object(s)",
           "\nRun make_se() or make_se_parse() to obtain the required columns",
           call. = FALSE)
    }
  })
  gather_join <- function(se) {
    SummarizedExperiment::assay(se) %>% data.frame(check.names = FALSE) %>% gather(ID, val) %>% left_join(.,
                                                                                                          data.frame(SummarizedExperiment::colData(se)), by = "ID")
  }
  df <- purrr::map_df(arglist, gather_join, .id = "var") %>%
    mutate(var = factor(var, levels = names(arglist)))
  p <- ggplot(df, aes(val, col = condition)) +
    geom_density(na.rm = TRUE, show.legend = F) +
    stat_density(geom="line", position ="identity", na.rm = TRUE) +
    facet_wrap(~var, ncol = 1) +
    labs(x = expression(log[2] ~ "Intensity"), y = "Density") +
    theme_DEP1(basesize = basesize)

  if (length(unique(se$condition)) <= 8) {
    p <- p + scale_color_brewer(palette = "Dark2")
  } else if (length(unique(se$condition)) <= 12) {
    p <- p + scale_color_brewer(palette = "Set3")
  } else if (length(unique(se$condition)) <= 25) {
    p <- p + scale_fill_manual(
      values = c(
        "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
        "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
        "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
        "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
        "green1", "yellow4", "yellow3", "darkorange4", "brown"
      )
    )
  } else {
    p <- p + scale_fill_manual(values = rainbow(length(unique(se$condition))))
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



#' Generates density plots for a SummarizedExperiment object
#'
#' @param se1 A SummarizedExperiment object with normalized data.
#' @param se2 (optional) A second SummarizedExperiment object with normalized data.
#' @param se3 (optional) A third SummarizedExperiment object with normalized data.
#' @param title1 A title for the first density plot.
#' @param title2 (optional) A title for the second density plot.
#' @param title3 (optional) A title for the third density plot.
#' @param group A logical value indicating whether to group the plots by condition (TRUE) or not (FALSE). If TRUE, \code{\link{prot.plot_imputation}} is called.
#' @param basesize The base font size for the plots.
#' @param plot A logical value indicating whether to plot the plots (TRUE) or not (FALSE).
#' @param export A logical value indicating whether to export the plots (TRUE) or not (FALSE).
#'
#' @export
#'
#' @return (Invisibly) returns a ggplot object.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices pdf
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @importFrom gridExtra marrangeGrob
#' @importFrom ggplot2 ggplot aes stat_density geom_line geom_density ggtitle scale_color_manual scale_fill_manual xlab ylab scale_linetype_manual theme
#'
prot.plot_density <- function(se1,
                              se2 = NULL,
                              se3 = NULL,
                              title1 = "",
                              title2 = "",
                              title3 = "",
                              group = FALSE,
                              basesize = 15,
                              plot = TRUE,
                              export = FALSE) {
  ls <- list()
  if (group == F) {
    # Create data frame from SummarizedExperiment with normalized data
    se.df1 <- data.frame(SummarizedExperiment::assay(se1), check.names = FALSE) %>% tibble::rownames_to_column()
    if (!is.null(se2)){
      se.df2 <- data.frame(SummarizedExperiment::assay(se2), check.names = FALSE) %>% tibble::rownames_to_column()
    }
    if (!is.null(se3)){
      se.df3 <- data.frame(SummarizedExperiment::assay(se3), check.names = FALSE) %>% tibble::rownames_to_column()
    }
    # Convert into long format
    se.tall1 <-
      se.df1 %>% pivot_longer(c(2:ncol(se.df1)), names_to = "colname", values_to = "val")
    if (!is.null(se2)){
      se.tall2 <-
        se.df2 %>% pivot_longer(c(2:ncol(se.df2)), names_to = "colname", values_to = "val")
    }
    if (!is.null(se3)){
      se.tall3 <-
        se.df3 %>% pivot_longer(c(2:ncol(se.df3)), names_to = "colname", values_to = "val")
    }
    if (length(unique(se1$condition)) <= 8) {
      pal <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
    } else if (length(unique(se1$condition)) <= 12) {
      pal <- RColorBrewer::brewer.pal(n = 12, name = "Set3")
    } else if (length(unique(se$condition)) <= 25) {
      p <- p + scale_fill_manual(
        values = c(
          "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
          "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
          "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
          "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
          "green1", "yellow4", "yellow3", "darkorange4", "brown"
        )
      )
    } else {
      p <- p + scale_fill_manual(values = rainbow(length(unique(se$condition))))
    }
    pal <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
    colors <- c()
    for (i in 1:length(unique(se1$condition))) {
      colors <-
        append(colors, (replicate(sum(
          se1$condition == unique(se1$condition[i])
        ),
        pal[i])))
    }

    linetype.ls <-
      c("solid", "dashed", "dotted", "dotdash", "twodash", "longdash")
    linetypes <- c()
    for (i in 1:length(se1$condition)) {
      linetypes[i] <- linetype.ls[as.double(se1$replicate[i])]
    }

    # Create a density plot of all samples
    densplot1 <- ggplot(se.tall1) +
      stat_density(aes(x = val, col = colname, lty = colname), lwd = 0.7, geom="line", position ="identity") +
      theme_DEP1(basesize = basesize) +
      theme(legend.key.width = unit(1.2, 'cm')) +
      ggtitle(title1) +
      scale_color_manual(values = colors) +
      xlab(expression(log[2]("Abundance"))) + ylab("Density") +
      scale_linetype_manual(values = linetypes)
    densplot1$labels$colour <- "Sample"
    densplot1$labels$linetype <- "Sample"

    ls[[1]] <- densplot1

    if (!is.null(se2)){
      densplot2 <- ggplot(se.tall2) +
        stat_density(aes(x = val, col = colname, lty = colname), lwd = 0.7, geom="line", position ="identity") +
        theme_DEP1(basesize = basesize) +
        theme(legend.key.width = unit(1.2, 'cm')) +
        ggtitle(title2) +
        scale_color_manual(values = colors) +
        xlab(expression(log[2]("Abundance"))) + ylab("Density") +
        scale_linetype_manual(values = linetypes)
      densplot2$labels$colour <- "Sample"
      densplot2$labels$linetype <- "Sample"
    }
    if (!is.null(se3)){
      densplot3 <- ggplot(se.tall3) +
        stat_density(aes(x = val, col = colname, lty = colname), lwd = 0.7, geom="line", position ="identity") +
        theme_DEP1(basesize = basesize) +
        theme(legend.key.width = unit(1.2, 'cm')) +
        ggtitle(title3) +
        scale_color_manual(values = colors) +
        xlab(expression(log[2]("Abundance"))) + ylab("Density") +
        scale_linetype_manual(values = linetypes)
      densplot3$labels$colour <- "Sample"
      densplot3$labels$linetype <- "Sample"
    }
    if (is.null(se2) && is.null(se3)){
      ls[[1]] <- densplot1
      ls.dens <- 1
      lay <- rbind( c(1,1,1), c(1,1,1) )
      width = 12
      height = 6
    } else if (!is.null(se2) && is.null(se3)){
      ls[[1]] <- densplot1
      ls[[2]] <- densplot2
      ls.dens <- c(1:2)
      lay <- rbind( c(1,1), c(1,1), c(2,2), c(2,2) )
      width = 12
      height = 12
    } else if (!is.null(se2) && !is.null(se3)){
      ls[[1]] <- densplot1
      ls[[2]] <- densplot2
      ls[[3]] <- densplot3
      ls.dens <- c(1:3)
      lay <- rbind( c(1,1,1), c(1,1,1), c(2,2,2), c(2,2,2), c(3,3,3), c(3,3,3) )
      width = 12
      height = 18
    }

    if (export == TRUE){
      dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
      message(paste0("Exporting density plots to:\n\"", getwd(), "/Plots/DensityPlot.pdf\" and \".../DensityPlot.png\""))
      grDevices::pdf("Plots/DensityPlot.pdf", height = height, width = width)
      suppressWarnings(print(gridExtra::marrangeGrob(ls[ls.dens], layout_matrix = lay, top = NULL)))
      grDevices::dev.off()
      grDevices::png("Plots/DensityPlot.png",
                     width = width, height = height, units = 'in', res = 300)
      suppressWarnings(print(gridExtra::marrangeGrob(ls[ls.dens], layout_matrix = lay, top = NULL)))
      grDevices::dev.off()
    }

    if (plot == TRUE){
      suppressWarnings(print(gridExtra::marrangeGrob(ls[ls.dens], layout_matrix = lay, top = NULL)))
    }

  } else {
    if (is.null(se2) && is.null(se3)){
      p <- prot.plot_imputation(se1, export = F, plot = F)
      width = 12
      height = 6
    } else if (!is.null(se2) && is.null(se3)){
      p <- prot.plot_imputation(se1, se2, export = F, plot = F)
      width = 12
      height = 12
    } else if (!is.null(se2) && !is.null(se3)){
      p <- prot.plot_imputation(se1, se2, se3, export = F, plot = F)
      width = 12
      height = 18
    }

    if (export == TRUE){
      dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
      message(paste0("Exporting grouped density plots to:\n\"", getwd(), "/Plots/DensityPlot_grp.pdf\" and \".../DensityPlot_grp.png\""))
      grDevices::pdf("Plots/DensityPlot_grp.pdf", width=width, height=height)
      print(p)
      grDevices::dev.off()
      grDevices::png("Plots/DensityPlot_grp.png",
                     width=width, height=height, units = 'in', res = 300)
      print(p)
      grDevices::dev.off()
    }

    if (plot == TRUE){
      print(p)
    }
    invisible(p)
  }
}



#' Plot missing values in SummarizedExperiment
#'
#' This function takes a SummarizedExperiment object and creates a binary heatmap of missing values in the assay.
#'
#' @param se A SummarizedExperiment object
#' @param plot Logical. Whether to return the plot.
#' @param export Logical. Whether to export the heatmap in pdf and png format.
#' @param fontsize Numeric. Font size for the column names in the heatmap.
#' @return A binary heatmap of missing values in the SummarizedExperiment.
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom SummarizedExperiment assay
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom grid gpar
#' @importFrom grDevices pdf png dev.off
#'
#'
prot.plot_missval <- function (se, plot = TRUE, export = FALSE, fontsize = 12)
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  se_assay <- SummarizedExperiment::assay(se)
  if (!any(is.na(se_assay))) {
    stop("No missing values in '", deparse(substitute(se)),
         "'", call. = FALSE)
  }
  # Make assay data binary (1 = valid value, 0 = missing value)
  df <- se_assay %>% data.frame(check.names = FALSE)
  missval <- df[apply(df, 1, function(x) any(is.na(x))), ]
  missval <- ifelse(is.na(missval), 0, 1)
  # Plot binary heatmap
  ht2 = ComplexHeatmap::Heatmap(missval, col = c("white", "black"),
                                column_names_side = "top", show_row_names = FALSE,
                                show_column_names = TRUE, name = "Missing values pattern",
                                column_names_gp = grid::gpar(fontsize = fontsize), heatmap_legend_param = list(at = c(0,
                                                                                                                      1), labels = c("Missing value", "Valid value")))
  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
    message(paste0("Exporting Missing values pattern to:\n\"", getwd(), "/Plots/MissValPlot.pdf\" and \".../MissValPlot.png\""))

    grDevices::pdf("Plots/MissValPlot.pdf")
    ComplexHeatmap::draw(ht2, heatmap_legend_side = "top")
    grDevices::dev.off()

    grDevices::png("Plots/MissValPlot.png",
                   width = 6, height = 6, units = 'in', res = 300)
    ComplexHeatmap::draw(ht2, heatmap_legend_side = "top")
    grDevices::dev.off()
  }

  if (plot == TRUE){
    ComplexHeatmap::draw(ht2, heatmap_legend_side = "top")
  }
}

#' Plot detectability for a SummarizedExperiment
#'
#' Plots the group-wise density of proteins with missing values
#'
#' @param se A SummarizedExperiment object
#' @param basesize Base font size for the plots
#' @param plot Logical, whether to plot the figure
#' @param export Logical, whether to export the figure
#'
#' @return Plots of the density distributions and cumulative fraction of proteins with and without missing values
#'
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom SummarizedExperiment assay
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_density stat_density geom_line labs guides guide_legend theme_minimal
#' @importFrom grDevices pdf png
#' @importFrom gridExtra grid.arrange
#' @importFrom stringr str_replace
prot.plot_detect <- function (se, basesize = 10, plot = TRUE, export = FALSE)
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  se_assay <- SummarizedExperiment::assay(se)
  if (!any(is.na(se_assay))) {
    stop("No missing values in '", deparse(substitute(se)),
         "'", call. = FALSE)
  }
  df <- se_assay %>% data.frame(check.names = FALSE) %>% tibble::rownames_to_column() %>%
    gather(ID, val, -rowname)
  stat <- df %>% group_by(rowname) %>% summarize(mean = mean(val,
                                                             na.rm = TRUE), missval = any(is.na(val)))
  cumsum <- stat %>% group_by(missval) %>% arrange(mean) %>%
    mutate(num = 1, cs = cumsum(num), cs_frac = cs/n())
  p1 <- ggplot(stat, aes(mean, col = missval)) +
    geom_density(na.rm = TRUE, show.legend = FALSE) +
    stat_density(geom="line", position ="identity") +
    labs(x = expression(log[2] ~ "Intensity"), y = "Density") +
    guides(col = guide_legend(title = paste0(
      "Missing values\n(",
      round(
        100 * nrow(cumsum[cumsum$missval == TRUE, ]) / nrow(cumsum[cumsum$missval == FALSE, ]),
        digits = 1
      ),
      "% missing)"
    ))) +
    theme_minimal(base_size = basesize)
  p2 <- ggplot(cumsum, aes(mean, cs_frac, col = missval)) +
    geom_line() + labs(x = expression(log[2] ~ "Intensity"),
                       y = "Cumulative fraction") +
    guides(col = guide_legend(title = paste0(
      "Missing values\n(",
      round(
        100 * nrow(cumsum[cumsum$missval == TRUE, ]) / nrow(cumsum[cumsum$missval ==
                                                                     FALSE, ]),
        digits = 1
      ),
      "% missing)"
    ))) +
    theme_minimal(base_size = basesize)

  df_per_group <- df
  df_per_group$ID <- str_replace(df_per_group$ID, "_[:digit:]+$", "")
  stat_per_group <- df_per_group %>% group_by(rowname, ID) %>% summarize(mean = mean(val,
                                                                                     na.rm = TRUE), missval = any(is.na(val)))
  cumsum_per_group <- stat_per_group %>% group_by(missval) %>% arrange(mean) %>%
    mutate(num = 1, cs = cumsum(num), cs_frac = cs/n())

  p1_per_group <- ggplot(stat_per_group, aes(mean, col = missval)) +
    geom_density(na.rm = TRUE, show.legend = FALSE) +
    stat_density(geom = "line", position = "identity") +
    labs(x = expression(log[2] ~ "Intensity"), y = "Density") +
    guides(col = guide_legend(title = paste0(
      "Missing values\n(",
      round(
        100 * nrow(cumsum_per_group[cumsum_per_group$missval == TRUE, ]) / nrow(cumsum_per_group[cumsum_per_group$missval ==
                                                                                                   FALSE, ]),
        digits = 1
      ),
      "% missing)"
    ))) +
    theme_DEP1(basesize = basesize)
  p2_per_group <-
    ggplot(cumsum_per_group, aes(mean, cs_frac, col = missval)) +
    geom_line() + labs(x = expression(log[2] ~ "Intensity"),
                       y = "Cumulative fraction") +
    guides(col = guide_legend(title = paste0(
      "Missing values\n(",
      round(
        100 * nrow(cumsum_per_group[cumsum_per_group$missval == TRUE, ]) / nrow(cumsum_per_group[cumsum_per_group$missval ==
                                                                                                   FALSE, ]),
        digits = 1
      ),
      "% missing)"
    ))) +
    theme_DEP1(basesize = basesize)



  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
    message(paste0("Exporting density distributions and cumulative fraction of proteins with and without missing values to:\n\"", getwd(), "/Plots/DensMissPlot.pdf\" and \".../DensMissPlot.png\""))

    grDevices::pdf("Plots/DensMissPlot.pdf")
    gridExtra::grid.arrange(p1, p2, ncol = 1)
    grDevices::dev.off()

    grDevices::png("Plots/DensMissPlot.png",
                   width = 6, height = 6, units = 'in', res = 300)
    gridExtra::grid.arrange(p1, p2, ncol = 1)
    grDevices::dev.off()
  }

  if (plot == TRUE){
    gridExtra::grid.arrange(p1, p2, ncol = 1)
  }
}

#' Plot Results of Principal Component Analysis
#'
#' Performs PCA on a given SummarizedExperiment object with log2-transformed protein abundances and plots the results.
#'
#' @param dep SummarizedExperiment object
#' @param x x-axis PC (default: 1)
#' @param y y-axis PC (default: 2)
#' @param indicate features to indicate in the plot ("condition" or "replicate")
#' @param title title of the plot (default: "PCA plot - top \code{n} variable proteins")
#' @param label whether to label points in the plot (default: FALSE)
#' @param n number of variables to take into account, sorted according to their variance in descending order. Only the \code{n} most variable proteins are used to perform PCA. default: number of columns in \code{dep})
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
#' @return (invisibly) a data frame containing the PCA coordinates
#' @importFrom ggplot2 ggplot  aes labs coord_fixed geom_point facet_wrap scale_shape_manual scale_fill_brewer scale_fill_manual scale_color_manual
#' @export
prot.plot_pca <- function (dep,
                           x = 1,
                           y = 2,
                           indicate = c("condition", "replicate"),
                           title = NULL,
                           label = FALSE,
                           n = ncol(t(SummarizedExperiment::assay(dep))),
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
                           export = FALSE)
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
    inherits(dep, "SummarizedExperiment"),
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
  if (x > ncol(dep) || y > ncol(dep)) {
    stop(
      paste0(
        "'x' and/or 'y' arguments are not valid\n",
        "Run prot.plot_pca() with 'x' and 'y' <= ",
        ncol(dep),
        "."
      ),
      call. = FALSE
    )
  }
  if (n > nrow(dep)) {
    stop(
      paste0(
        "'n' argument is not valid.\n",
        "Run prot.plot_pca() with 'n' <= ",
        nrow(dep),
        "."
      ),
      call. = FALSE
    )
  }
  columns <- colnames(SummarizedExperiment::colData(dep))
  if (!is.null(indicate)) {
    if (length(indicate) > 3) {
      stop(
        "Too many features in 'indicate'\n        Run plot_pca() with a maximum of 3 indicate features"
      )
    }
    if (any(!indicate %in% columns)) {
      stop(
        paste0(
          "'",
          paste0(indicate, collapse = "' and/or '"),
          "' column(s) is/are not present in ",
          deparse(substitute(dep)),
          ".\nValid columns are: '",
          paste(columns, collapse = "', '"),
          "'."
        ),
        call. = FALSE
      )
    }
  }
  var <- apply(SummarizedExperiment::assay(dep), 1, sd)
  df <- SummarizedExperiment::assay(dep)[order(var, decreasing = TRUE)[seq_len(n)],]
  pca <- stats::prcomp(t(df), scale = FALSE)
  pca_df <- pca$x %>% data.frame(check.names = FALSE) %>% tibble::rownames_to_column() %>%
    left_join(., data.frame(SummarizedExperiment::colData(dep), check.names = FALSE), by = c(rowname = "ID"))
  percent <- round(100 * pca$sdev ^ 2 / sum(pca$sdev ^ 2), 1)
  if (!is.null(title)){
    title <- as.character(title)
  } else {
    title <- paste0("PCA plot - top ", n, " variable proteins")
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
      scale_shape_manual(values=c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21, 22, 23, 24, 25))

    if (length(unique(dep$condition)) <= 8) {
      p <- p + scale_fill_brewer(palette = "Dark2", guide = guide_legend(override.aes = list(shape = 21)))+
        scale_color_manual(values=c(rep("black", length(unique(dep$condition)))))
    } else if (length(unique(dep$condition)) <= 12) {
      p <- p + scale_fill_brewer(palette = "Set3", guide = guide_legend(override.aes = list(shape = 21))) +
        scale_color_manual(values=c(rep("black", length(unique(dep$condition)))))
    } else if (length(unique(dep$condition)) <= 25) {
      p <- p + scale_fill_manual(
        values = c(
          "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
          "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
          "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
          "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
          "green1", "yellow4", "yellow3", "darkorange4", "brown"
        )
      )
    } else {
      p <- p + scale_fill_manual(values = rainbow(length(unique(dep$condition))))
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

  df <- pca_df %>% select(rowname, paste0("PC", c(x, y)),
                          match(indicate, colnames(pca_df)))
  colnames(df)[1] <- "sample"
  invisible(df)
}



#' Plots a volcano plot for a given contrast
#'
#' @param dep A SummarizedExperiment object
#' @param contrast The contrast of interest in the format "conditionA_vs_conditionB"
#' @param label_size The size of labels for the protein names
#' @param alpha The alpha value used to determine significance.
#' @param lfc The log2 fold change threshold.
#' @param add_names Logical, whether to add protein names to the plot
#' @param adjusted Logical, whether to plot adjusted p-values. The \code{alpha} threshold is applied to adjusted p values even if \code{adjusted = FALSE}.
#' @param plot Logical, whether to return the volcano plot
#' @param export Logical, whether to export the volcano plot as PND and PDF file
#'
#' @return (invisibly) A data frame with columns for protein, log2 fold change, and -log10 p-value
#'
#' @export
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom assertthat assert_that
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline labs scale_shape scale_colour_manual labs xlab ylab theme_bw theme coord_cartesian
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggh4x force_panelsizes
#' @importFrom scales pretty_breaks
#' @importFrom stringr str_replace
#' @importFrom grDevices pdf png dev.off
#' @importFrom dplyr filter arrange mutate select
prot.plot_volcano <-
  function (dep,
            contrast,
            label_size = 3,
            alpha = 0.05,
            lfc = 2,
            add_names = TRUE,
            adjusted = FALSE,
            plot = TRUE,
            export = FALSE)
  {
    if (is.integer(label_size))
      label_size <- as.numeric(label_size)
    assertthat::assert_that(
      inherits(dep, "SummarizedExperiment"),
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
    row_data <- SummarizedExperiment::rowData(dep, use.names = FALSE)
    if (any(!c("name", "ID") %in% colnames(row_data))) {
      stop(
        paste0(
          "'name' and/or 'ID' columns are not present in '",
          deparse(substitute(dep)),
          "'.\nRun make_unique() to obtain required columns."
        ),
        call. = FALSE
      )
    }
    if (length(grep("_p.adj|_diff", colnames(row_data))) <
        1) {
      stop(
        paste0(
          "'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
          deparse(substitute(dep)),
          "'.\nRun test_diff() to obtain the required columns."
        ),
        call. = FALSE
      )
    }
    if (length(grep("_significant", colnames(row_data))) <
        1) {
      stop(
        paste0(
          "'[contrast]_significant' columns are not present in '",
          deparse(substitute(dep)),
          "'.\nRun add_rejections() to obtain the required columns."
        ),
        call. = FALSE
      )
    }
    if (length(grep(paste(contrast, "_diff", sep = ""),
                    colnames(row_data))) == 0) {
      valid_cntrsts <-
        row_data %>% data.frame(check.names = FALSE) %>% select(ends_with("_diff")) %>%
        colnames(.) %>% gsub("_diff", "", .)
      valid_cntrsts_msg <- paste0("Valid contrasts are: '",
                                  paste0(valid_cntrsts, collapse = "', '"),
                                  "'")
      stop(
        "Not a valid contrast, please run `prot.plot_volcano()` with a valid contrast as argument\n",
        valid_cntrsts_msg,
        call. = FALSE
      )
    }
    diff <- match(paste(contrast, "_diff", sep = ""),
                  colnames(row_data))
    if (adjusted) {
      p_values <- match(paste(contrast, "_p.adj", sep = ""),
                        colnames(row_data))
    }
    else {
      p_values <- match(paste(contrast, "_p.val", sep = ""),
                        colnames(row_data))
    }
    signif <- match(paste(contrast, "_significant", sep = ""),
                    colnames(row_data))

    df <- data.frame(
      x = row_data[, diff],
      y = -log10(row_data[, p_values]),
      significant = row_data[, signif],
      name = row_data$name,
      check.names = FALSE
    ) %>%
      dplyr::filter(!is.na(significant)) %>% arrange(significant)

    # create new columns for shapes "circle" and "triangle" if proteins lie within or outside (for proteins with very small adj.p value < 3.162278e-05) of the y bounds, respectively.
    if (adjusted) {
      df$shape <- ifelse(df$y > 4.9, "triangle", "circle")
    }
    else{
      df$shape <- "circle"
    }
    # change the -Log10(q value) of proteins exceeding the y plot bound at 5.0 to 4.9 so that they are displayed at the axis border
    if (adjusted) {
      if (max(df$y) > 5.0) {
        df$y[df$y > 4.9] <- 4.9
      }
    }

    df_volcano <- df %>%
      # Categorize proteins based on the log2(fold change) and q value thresholds
      mutate(Proteins = factor(
        case_when(
          x >= lfc & significant == "TRUE" ~ "Upregulated",
          x <= -lfc & significant == "TRUE" ~ "Downregulated",
          TRUE ~ "Not significant"
        )
      ))
    # Plot with points colored according to the thresholds
    plot_volcano <- df_volcano %>%
      ggplot(aes(x, y, group = name, colour = Proteins)) +
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
        labels = c("Downregulated", "Not significant", "Upregulated"),
        values = c(
          "Upregulated" = "#e31f26",
          "Not significant" = "#525352",
          "Downregulated" = "#387fb9")) +
      labs(title =  str_replace(contrast, "_vs_", " vs. ")) +
      xlab(expression(log[2]("Fold change"))) + ylab(expression(-Log[10]("q-value"))) +       # Label the axes
      theme_bw(base_size = 20) + #Set the theme, define label font size
      theme(legend.position = "bottom", title = element_text(size = exp(-nchar(str_replace(contrast, "_vs_", " vs. "))/40)*38),
            axis.title = element_text(size = 22),
            legend.title = element_text(size = 22)) +
      ggh4x::force_panelsizes(rows = unit(6.5, "in"),
                              cols = unit(6.5, "in")) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = max(abs(df$x)) + 1)) +
      guides(color = guide_legend(reverse = F))

    if (adjusted) {
      plot_volcano <- plot_volcano + coord_cartesian(
        xlim = c(-max(abs(df$x)) - 0.3,  # Lower x bound [defined dynamically based on the largest absolute log2(fc) value]
                 max(abs(df$x)) + 0.3),
        # Upper x bound [defined dynamically based on the largest absolute log2(fc) value]
        ylim = c((min(df$y) - 0.1), 5.0),
        expand = FALSE
      )     # Plot y bounds (the bottom axis is defined dynamically based on the smallest q value)
    }
    else{
      plot_volcano <- plot_volcano + coord_cartesian(
        xlim = c(-max(abs(df$x)) - 0.3,  # Lower x bound [defined dynamically based on the largest absolute log2(fc) value]
                 max(abs(df$x)) + 0.3),
        # Upper x bound [defined dynamically based on the largest absolute log2(fc) value]
        ylim = c((min(df$y) - 0.1), max(df$y)),
        expand = TRUE
      )     # Plot y bounds (the bottom axis is defined dynamically based on the smallest q value)
    }
    if (adjusted) {
      plot_volcano <- plot_volcano +
        geom_hline(yintercept = -log10(alpha),
                   linetype = "dashed",
                   alpha = 0.5)      # Add dotted line to indicate the adj.p value threshold, semi-transparent
    }
    if (add_names) {
      plot_volcano <- plot_volcano + ggrepel::geom_text_repel(
        data = dplyr::filter(df_volcano, significant),
        aes(label = name),
        size = label_size,
        box.padding = unit(0.25,
                           "lines"),
        point.padding = unit(0.1, "lines"),
        segment.size = 0.5,
        show.legend = FALSE,
        max.overlaps = 25
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
      colnames(df)[c(1, 2, 3)] <- c("protein", "log2_fold_change",
                                    "p_value_-log10")
      if (adjusted) {
        colnames(df)[3] <- "adjusted_p_value_-log10"
      }
    }
    invisible(df)
  }



#' Plot a correlation heatmap
#'
#' This function plots a Pearson correlation heatmap for a given SummarizedExperiment object.
#'
#' @param dep A SummarizedExperiment object
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
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom dplyr filter
#' @importFrom tibble rownames_to_column
#'
prot.plot_corrheatmap <- function (dep, significant = TRUE, lower = 0, upper = 1, pal = "PRGn",
                                   pal_rev = FALSE, indicate = NULL, font_size = 12, plot = TRUE,
                                   ...)
{
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.logical(significant), length(significant) == 1, is.numeric(lower),
                          length(lower) == 1, is.numeric(upper), length(upper) ==
                            1, is.character(pal), length(pal) == 1, is.logical(pal_rev),
                          length(pal_rev) == 1, is.numeric(font_size), length(font_size) ==
                            1, is.logical(plot), length(plot) == 1)
  if (!(lower >= -1 & upper >= -1 & lower <= 1 & upper <= 1)) {
    stop("'lower' and/or 'upper' arguments are not valid\n         Run prot.plot_corrheatmap() with 'lower' and 'upper' between -1 and 1",
         call. = FALSE)
  }
  pals <- RColorBrewer::brewer.pal.info %>% tibble::rownames_to_column() %>%
    dplyr::filter(category != "qual")
  if (!pal %in% pals$rowname) {
    stop("'", pal, "' is not a valid color panel",
         " (qualitative panels also not allowed)\n",
         "Run plot_pca() with one of the following 'pal' options: ",
         paste(pals$rowname, collapse = "', '"), "'",
         call. = FALSE)
  }
  if (!is.null(indicate)) {
    assertthat::assert_that(is.character(indicate))
    col_data <- SummarizedExperiment::colData(dep) %>% as.data.frame()
    columns <- colnames(col_data)
    if (any(!indicate %in% columns)) {
      stop("'", paste0(indicate, collapse = "' and/or '"),
           "' column(s) is/are not present in ", deparse(substitute(dep)),
           ".\nValid columns are: '", paste(columns,
                                            collapse = "', '"), "'.", call. = FALSE)
    }
    anno <- SummarizedExperiment::colData(dep) %>% data.frame() %>% select(all_of(indicate))
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
  if (significant) {
    if (!"significant" %in% colnames(SummarizedExperiment::rowData(dep, use.names = FALSE))) {
      stop("'significant' column is not present in '",
           deparse(substitute(dep)), "'\nRun add_rejections() to obtain the required column",
           call. = FALSE)
    }
    dep <- dep[!is.na(SummarizedExperiment::rowData(dep, use.names = FALSE)$significant), ]
    dep <- dep[SummarizedExperiment::rowData(dep, use.names = FALSE)$significant, ]
  }
  cor_mat <- cor(SummarizedExperiment::assay(dep))
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
  df <- as.data.frame(cor_mat)
  invisible(df)
}


#' Plots a heatmap of the differentially expressed proteins
#'
#' @param dep A SummarizedExperiment object containing the normalized protein abundances and the associated metadata
#' @param type Type of heatmap to plot. Either "centered" or "contrast".
#' @param contrast (String or vector of strings) Analyze only significant proteins contained in the specified contrast(s). Default is NULL. Contrasts must be given in the form "ConditionA_vs_ConditionB"
#' @param show_all Show protein abundances of all conditions or only of those in the specified constrast(s). Default is TRUE
#' @param pal Color palette to use. Default is "RdBu"
#' @param kmeans Perform k-means clustering on the rows. Default is FALSE
#' @param k Number of clusters to use when kmeans is TRUE. Default is 6
#' @param col_limit Specify a custom range for the color scale of the heatmap. Default is NA, in which case the minimum and maximum of the color scale are calculated as the 5% and 95% percentiles.
#' @param indicate Annotate the heatmap with the specified column from the metadata.
#' @param clustering_distance Distance metric used for clustering. Default is "euclidean".
#' @param show_row_names Show the protein names on the left side of the heatmap. Default is FALSE.
#' @param row_font_size Font size for the protein names. Default is 6
#' @param col_font_size Font size for the condition names. Default is 10
#' @param plot Plot the heatmap. Default is TRUE
#' @param export Export the heatmap as pdf and png. Default is TRUE
#' @param ... Other parameters passed to the ComplexHeatmap::Heatmap function
#' @return (Invisibly) returns a data frame with the proteins and the normalized abundances
#' @importFrom assertthat assert_that
#' @importFrom dplyr select filter
#' @importFrom stringr str_detect
#' @export
prot.plot_heatmap <- function (dep, type = c("centered", "contrast"),
                               contrast = NULL, # Analyze only significant proteins contained in the specified contrast(s)
                               show_all = TRUE, # Show protein abundances of all conditions or only of those in the specified constrast(s)
                               pal = "RdBu",
                               kmeans = FALSE, k = 6,
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
                               export = FALSE, ...)
{
  get_annotation_contrast <- function (dep, indicate, contrast = contrast_samples)
  {
    assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                            is.character(indicate))
    col_data <- SummarizedExperiment::colData(dep) %>% data.frame(check.names = FALSE)
    columns <- colnames(col_data)
    if (all(!indicate %in% columns)) {
      stop("'", paste0(indicate, collapse = "' and/or '"),
           "' column(s) is/are not present in ", deparse(substitute(dep)),
           ".\nValid columns are: '", paste(columns, collapse = "', '"),
           "'.", call. = FALSE)
    }
    if (any(!indicate %in% columns)) {
      indicate <- indicate[indicate %in% columns]
      warning("Only used the following indicate column(s): '",
              paste0(indicate, collapse = "', '"), "'")
    }
    anno <- dplyr::select(col_data, all_of(indicate))
    anno <- dplyr::filter(anno, str_detect(condition, paste(contrast, collapse = "|")))
    names <- colnames(anno)
    anno_col <- vector(mode = "list", length = length(names))
    names(anno_col) <- names
    for (i in names) {
      var = anno[[i]] %>% unique() %>% sort()
      if (length(var) == 1)
        cols <- c("black")
      if (length(var) == 2)
        cols <- c("orangered", "cornflowerblue")
      if (length(var) < 7 && length(var) > 2)
        cols <- RColorBrewer::brewer.pal(length(var), "Pastel1")
      if (length(var) > 7 && length(var) <= 12){
        cols <- RColorBrewer::brewer.pal(length(var), "Set3")
      } else {
        pal <- c(
          "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
                       "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
                       "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
                       "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
                       "green1", "yellow4", "yellow3", "darkorange4", "brown"
        )
        cols <-pal[1:length(var)]
      }


      names(cols) <- var
      anno_col[[i]] <- cols
    }
    ComplexHeatmap::HeatmapAnnotation(df = anno, col = anno_col, show_annotation_name = TRUE)
  }
  if(is.null(contrast) & show_all == FALSE){
    contrast <- SummarizedExperiment::rowData(dep, use.names = FALSE) %>% data.frame(check.names = FALSE) %>% select(ends_with("_diff")) %>%
      colnames(.) %>% gsub("_diff", "", .)
  }
  if (length(grep(paste(paste(contrast, collapse = "|"), "_diff", sep = ""),
                  colnames(SummarizedExperiment::rowData(dep, use.names = FALSE)))) == 0) {
    valid_cntrsts <-
      SummarizedExperiment::rowData(dep, use.names = FALSE) %>% data.frame(check.names = FALSE) %>% select(ends_with("_diff")) %>%
      colnames(.) %>% gsub("_diff", "", .)
    valid_cntrsts_msg <- paste0("Valid contrasts are: '",
                                paste0(valid_cntrsts, collapse = "', '"),
                                "'")
    stop(
      "Not a valid contrast, please run `prot.plot_heatmap()` with a valid contrast as argument\n",
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
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(type), is.logical(kmeans), is.numeric(k),
                          length(k) == 1, is.numeric(row_font_size), length(row_font_size) ==
                            1, is.numeric(col_font_size), length(col_font_size) ==
                            1, is.logical(plot), length(plot) == 1)
  type <- match.arg(type)
  clustering_distance <- match.arg(clustering_distance)
  row_data <- SummarizedExperiment::rowData(dep, use.names = FALSE)
  col_data <- SummarizedExperiment::colData(dep) %>% data.frame(check.names = FALSE)
  if (any(!c("label", "condition", "replicate") %in%
          colnames(col_data))) {
    stop(paste0("'label', 'condition' and/or 'replicate' columns are not present in '",
                deparse(substitute(dep)), "'"), call. = FALSE)
  }
  if (length(grep("_diff", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' columns are not present in '",
                deparse(substitute(dep)), "'.\nRun test_diff() to obtain the required columns."),
         call. = FALSE)
  }
  if (!"significant" %in% colnames(row_data)) {
    stop(paste0("'significant' column is not present in '",
                deparse(substitute(dep)), "'.\nRun add_rejections() to obtain the required column."),
         call. = FALSE)
  }
  if (!is.null(indicate) && type == "contrast") {
    warning("Heatmap annotation only applicable for type = 'centered'",
            call. = FALSE)
  }
  if (!is.null(indicate) && type == "centered" && show_all == TRUE) {
    ha1 <- get_annotation(dep, indicate)
  } else if (!is.null(indicate) && type == "centered" && show_all == FALSE) {
    ha1 <- get_annotation_contrast(dep, indicate, contrast = contrast_samples)
  } else {
    ha1 <- NULL
  }
  filtered <- dep[row_data$significant, ]
  if (nrow(SummarizedExperiment::assay(filtered)) == 0){
    stop("No proteins with significantly different abundance were found.")
  }
  if (any(is.na(SummarizedExperiment::assay(filtered)))) {
    warning("Missing values in '", deparse(substitute(dep)),
            "'. ", "Using clustering_distance = 'gower'",
            call. = FALSE)
    clustering_distance <- "gower"
    obs_NA <- TRUE
  } else {
    obs_NA <- FALSE
  }
  if (type == "centered") {
    SummarizedExperiment::rowData(filtered)$mean <- rowMeans(SummarizedExperiment::assay(filtered), na.rm = TRUE)
    df <- SummarizedExperiment::assay(filtered) - SummarizedExperiment::rowData(filtered, use.names = FALSE)$mean
    if (!is.null(contrast) && show_all == FALSE){
      df <- data.frame(df, check.names = FALSE)[,str_detect(colnames(df),paste(contrast_samples, collapse = "|"))]
    }
  }
  if (type == "contrast") {
    df <- SummarizedExperiment::rowData(filtered, use.names = FALSE) %>% data.frame(check.names = FALSE) %>%
      tibble::column_to_rownames(var = "name") %>% select(ends_with("_diff"))
    if(length(grep("_diff", colnames(SummarizedExperiment::rowData(filtered)))) == 1){
      df[[1]] <- df[[1]] %>% sort()
    }
    colnames(df) <- gsub("_diff", "", colnames(df)) %>%
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
        summarize(index = sum(row)/n()) %>% arrange(desc(index)) %>%
        pull(cluster) %>% match(seq_len(k), .)
      df_kmeans$cluster <- order[df_kmeans$cluster]
    }
    if (type == "contrast") {
      order <- data.frame(df, check.names = FALSE) %>% cbind(., cluster = df_kmeans$cluster) %>%
        gather(condition, diff, -cluster) %>% group_by(cluster) %>%
        summarize(row = mean(diff)) %>% arrange(desc(row)) %>%
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
    df_cent <- SummarizedExperiment::assay(filtered) - rowMeans(SummarizedExperiment::assay(filtered))
    col_lim <- ceiling(stats::quantile(df_cent, probs= 0.95)) - ceiling(stats::quantile(df_cent, probs= 0.05))
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
    w <- 6+length(dep$condition)/10
    if (show_row_names == TRUE){
      len = nrow(dep[SummarizedExperiment::rowData(dep)$significant, ]) / 12
    } else {
      len = nrow(dep[SummarizedExperiment::rowData(dep)$significant, ]) / 45
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
  }

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
    return <- df[unlist(ComplexHeatmap::row_order(ht1)), ]), silent = T)
  try(df <- data.frame(protein = row.names(return), return, check.names = FALSE) %>% mutate(order = row_number()), silent = T)

  invisible(df)

}




#' Plots a heatmap of all protein abundances
#'
#' @param se A SummarizedExperiment object
#' @param show_all Logical. Show protein abundances of all conditions or only of those in the specified constrast(s). Default is TRUE
#' @param pal Character. Color palette to use. Default is "RdBu"
#' @param kmeans Logical. Perform kmeans clustering on protein abundances. Default is FALSE
#' @param k Numeric. Number of clusters to use in kmeans clustering. Default is 6
#' @param col_limit Numeric. Maximum intensity for the heatmap. Default is NA
#' @param indicate Character. Column name of the annotation to indicate in the heatmap. Default is "condition"
#' @param clustering_distance Character. Distance metric to use for clustering. Default is "euclidean"
#' @param show_row_names Logical. Show protein names on the left of the heatmap. Default is FALSE
#' @param row_font_size Numeric. Font size for protein names. Default is 6
#' @param col_font_size Numeric. Font size for condition names. Default is 10
#' @param plot Logical. Plot the heatmap or not. Default is TRUE
#' @param export Logical. Export the heatmap to a pdf and png file. Default is TRUE
#' @param scale Logical. Scale the data before plotting. Default is TRUE
#' @param show_row_dend Logical. Show row dendrogram. Default is FALSE
#' @param ... Other parameters to pass to ComplexHeatmap::Heatmap
#'
#' @return A data frame with protein abundances
#'
#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr select filter
#' @importFrom stringr str_detect
prot.plot_heatmap_all <- function (se,
                                   show_all = TRUE, # Show protein abundances of all conditions or only of those in the specified constrast(s)
                                   pal = "RdBu",
                                   kmeans = FALSE,
                                   k = 6,
                                   col_limit = NA,
                                   indicate = "condition",
                                   clustering_distance = "euclidean",
                                   show_row_names = FALSE,
                                   row_font_size = 6,
                                   col_font_size = 10,
                                   plot = TRUE,
                                   export = FALSE,
                                   scale = TRUE,
                                   show_row_dend = FALSE,
                                   ...)
{


  if (is.integer(k)) k <- as.numeric(k)
  if (is.integer(row_font_size)) row_font_size <- as.numeric(row_font_size)
  if (is.integer(col_font_size)) col_font_size <- as.numeric(col_font_size)
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.logical(kmeans), is.numeric(k),
                          length(k) == 1, is.numeric(row_font_size), length(row_font_size) ==
                            1, is.numeric(col_font_size), length(col_font_size) ==
                            1, is.logical(plot), length(plot) == 1)
  clustering_distance <- match.arg(clustering_distance)
  row_data <- SummarizedExperiment::rowData(se, use.names = FALSE)
  col_data <- SummarizedExperiment::colData(se) %>% data.frame(check.names = FALSE)
  if (any(!c("label", "condition", "replicate") %in%
          colnames(col_data))) {
    stop(paste0("'label', 'condition' and/or 'replicate' columns are not present in '",
                deparse(substitute(se)), "'"), call. = FALSE)
  }
  ha1 <- get_annotation(se, indicate)
  obs_NA <- FALSE
  SummarizedExperiment::rowData(se)$mean <- rowMeans(SummarizedExperiment::assay(se), na.rm = TRUE)
  df <- SummarizedExperiment::assay(se) - SummarizedExperiment::rowData(se, use.names = FALSE)$mean


  df <- as.matrix(df)
  if(scale == TRUE){
    df <- scale(t(df)) %>% t()
  }
  if (kmeans && obs_NA) {
    warning("Cannot perform kmeans clustering with missing values",
            call. = FALSE)
    kmeans <- FALSE
  }
  if (kmeans && !obs_NA) {
    set.seed(1)
    df_kmeans <- kmeans(df, k)
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
    df_cent <- SummarizedExperiment::assay(se) - rowMeans(SummarizedExperiment::assay(se))
    col_lim <- ceiling(stats::quantile(df_cent, probs= 0.95, na.rm=TRUE)) - ceiling(stats::quantile(df_cent, probs= 0.05, na.rm=TRUE))
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
  legend <- "log2 Centered intensity"
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
    cluster_columns = row_clust,
    row_names_side = "left",
    column_names_side = "top",
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
    w <- 6+length(se$condition)/10
    if (show_row_names == TRUE){
      len = nrow(se[SummarizedExperiment::rowData(se)$significant, ]) / 12
    } else {
      len = nrow(se[SummarizedExperiment::rowData(se)$significant, ]) / 45
    }
    if (len < 6) { len = 6 }
    message(paste0("Exporting heat map to:\n\"",
                   getwd(),
                   "/Plots/HeatMap_all",
                   ".pdf",
                   " and \"HeatMap_all",
                   ".png\""))

    grDevices::pdf(paste0("Plots/HeatMap_all",
                          ".pdf"),
                   width = w, height = len)
    ComplexHeatmap::draw(ht1, heatmap_legend_side = "top")
    grDevices::dev.off()

    grDevices::png(paste0("Plots/HeatMap_all",
                          ".png"),
                   width = w, height = len, units = 'in', res = 300)
    ComplexHeatmap::draw(ht1, heatmap_legend_side = "top")
    grDevices::dev.off()
  }
  if (plot) {
    ComplexHeatmap::draw(ht1, heatmap_legend_side = "top")
  }
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
    df <- data.frame(protein = row.names(return), return, check.names = FALSE) %>% mutate(order = row_number())
    invisible(df)
}


#' SCREE plot
#'
#' prot.plot_screeplot() creates a SCREE plot of the explained variation in the principal components of a 'princomp' object
#'
#' @param pcaobj An object of class 'princomp'
#' @param components A numeric vector of principal component numbers to include in the plot
#' @param xlim A numeric vector of length 2 specifying the x-axis limits
#' @param ylim A numeric vector of length 2 specifying the y-axis limits. Default is c(0, 100)
#' @param xlab A character string for x-axis label
#' @param xlabAngle Numeric value for x-axis label angle
#' @param xlabhjust Numeric value for x-axis label horizontal justification
#' @param xlabvjust Numeric value for x-axis label vertical justification
#' @param ylab A character string for y-axis label
#' @param ylabAngle Numeric value for y-axis label angle
#' @param ylabhjust Numeric value for y-axis label horizontal justification
#' @param ylabvjust Numeric value for y-axis label vertical justification
#' @param axisLabSize Numeric value for axis label size
#' @param title A character string for the plot title
#' @param subtitle A character string for the plot subtitle
#' @param caption A character string for the plot caption
#' @param titleLabSize Numeric value for title label size
#' @param subtitleLabSize Numeric value for subtitle label size
#' @param captionLabSize Numeric value for caption label size
#' @param colBar A character string for the color of the bar chart
#' @param drawCumulativeSumLine Logical value to indicate if the cumulative sum line should be plotted
#' @param colCumulativeSumLine A character string for the color of the cumulative sum line
#' @param sizeCumulativeSumLine Numeric value for the size of the cumulative sum line
#' @param drawCumulativeSumPoints Logical value to indicate if the cumulative sum points should be plotted
#' @param colCumulativeSumPoints A character string for the color of the cumulative sum points
#' @param sizeCumulativeSumPoints Numeric value for the size of the cumulative sum points
#' @param hline A numeric value for horizontal line
#' @param hlineType A character string for the type of horizontal line
#' @param hlineCol A character string for the color of the horizontal line
#' @param hlineWidth Numeric value for the width of the horizontal line
#' @param vline A numeric value for vertical line
#' @param vlineType A character string for the type of vertical line
#' @param vlineCol A character string for the color of the vertical line
#' @param vlineWidth Numeric value for the width of the vertical line
#' @param gridlines.major Logical value to indicate if major gridlines should be plotted
#' @param gridlines.minor Logical value to indicate if minor gridlines should be plotted
#' @param borderWidth Numeric value for the width of the border
#' @param borderColour A character string for the color of the border
#' @param plot Logical value to indicate if the plot should be plotted
#' @param export Logical value to indicate if the plot should be exported as PNG and PDF
#' @return A SCREE plot
#' @importFrom ggplot2 ggplot aes element_blank element_line element_rect geom_bar geom_hline geom_line geom_point geom_vline labs theme theme_bw xlab xlim ylab ylim
#' @export
prot.plot_screeplot <- function (pcaobj, components, xlim = NULL,
                                 ylim = c(0, 100), xlab = "Principal component", xlabAngle = 90,
                                 xlabhjust = 0.5, xlabvjust = 0.5, ylab = "Explained variation (%)",
                                 ylabAngle = 0, ylabhjust = 0.5, ylabvjust = 0.5, axisLabSize = 16,
                                 title = "SCREE plot", subtitle = "", caption = "",
                                 titleLabSize = 16, subtitleLabSize = 12, captionLabSize = 12,
                                 colBar = "dodgerblue", drawCumulativeSumLine = TRUE,
                                 colCumulativeSumLine = "red2", sizeCumulativeSumLine = 1.5,
                                 drawCumulativeSumPoints = TRUE, colCumulativeSumPoints = "red2",
                                 sizeCumulativeSumPoints = 2, hline = NULL, hlineType = "longdash",
                                 hlineCol = "black", hlineWidth = 0.4, vline = NULL,
                                 vlineType = "longdash", vlineCol = "black", vlineWidth = 0.4,
                                 gridlines.major = TRUE, gridlines.minor = TRUE, borderWidth = 0.8,
                                 borderColour = "black", plot = TRUE, export = FALSE)
{
  PC <- Variance <- NULL
  components <- PCAtools::getComponents(pcaobj)[1:length(pcaobj$components)]
  if(length(PCAtools::getComponents(pcaobj)) > 10){
    components <- components[1:10]
  }
  plotobj <- data.frame(components, pcaobj$variance[components], check.names = FALSE)
  colnames(plotobj) <- c("PC", "Variance")
  plotobj$PC <- factor(plotobj$PC, levels = plotobj$PC[seq_along(plotobj$PC)])
  th <- theme_bw(base_size = 24) + theme(legend.background = element_rect(),
                                         plot.title = element_text(angle = 0, size = titleLabSize,
                                                                   face = "bold", vjust = 1), plot.subtitle = element_text(angle = 0,
                                                                                                                           size = subtitleLabSize, face = "plain", vjust = 0.5),
                                         plot.caption = element_text(angle = 0, size = captionLabSize,
                                                                     face = "plain", vjust = 0.5), axis.text.x = element_text(angle = xlabAngle,
                                                                                                                              size = axisLabSize, hjust = xlabhjust, vjust = xlabvjust),
                                         axis.text.y = element_text(angle = ylabAngle, size = axisLabSize,
                                                                    hjust = ylabhjust, vjust = ylabvjust), axis.title = element_text(size = axisLabSize),
                                         legend.position = "none")
  p <- ggplot(plotobj, aes(x = PC)) + th + geom_bar(aes(y = Variance),
                                                    fill = colBar, width = 0.8, stat = "identity",
                                                    show.legend = FALSE)
  if (drawCumulativeSumLine == TRUE) {
    p <- p + geom_line(aes(y = cumsum(Variance), group = 1),
                       colour = colCumulativeSumLine, size = sizeCumulativeSumLine,
                       show.legend = FALSE)
  }
  if (drawCumulativeSumPoints == TRUE) {
    p <- p + geom_point(aes(y = cumsum(Variance)),
                        colour = colCumulativeSumPoints, size = sizeCumulativeSumPoints,
                        show.legend = FALSE)
  }
  p <- p + xlab(xlab) + ylab(ylab)
  if (!is.null(xlim)) {
    p <- p + xlim(xlim[1], xlim[2])
  }
  if (!is.null(ylim)) {
    p <- p + ylim(ylim[1], ylim[2])
  }
  p <- p + labs(title = title, subtitle = subtitle, caption = caption)
  if (!is.null(vline)) {
    p <- p + geom_vline(xintercept = vline, linetype = vlineType,
                        colour = vlineCol, size = vlineWidth)
  }
  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = hlineType,
                        colour = hlineCol, size = hlineWidth)
  }
  p <- p + theme(panel.border = element_rect(colour = borderColour,
                                             fill = NA, size = borderWidth))
  if (gridlines.major == TRUE) {
    p <- p + theme(panel.grid.major = element_line())
  }
  else {
    p <- p + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor == TRUE) {
    p <- p + theme(panel.grid.minor = element_line())
  }
  else {
    p <- p + theme(panel.grid.minor = element_blank())
  }
  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
    message(paste0("Exporting scree plot to:\n\"", getwd(), "/Plots/ScreePlot.pdf\" and \".../ScreePlot.png\""))

    grDevices::pdf("Plots/ScreePlot.pdf")
    print(p)
    grDevices::dev.off()

    grDevices::png("Plots/ScreePlot.png",
                   width = 6, height = 6, units = 'in', res = 300)
    print(p)
    grDevices::dev.off()
  }
  if (plot == TRUE) {
    return(p)
  }
  else if (plot == FALSE) {
    p
  }
}




#' Plot the loadings of a principal components analysis
#'
#' @param pcaobj An object that contains the results of a PCA analysis.
#' @param components The components to plot. Defaults to the first 5 components.
#' @param rangeRetain a numeric value between 0 and 1 indicating the proportion of loading values to retain. Default is 0.05.
#' @param absolute a logical value indicating whether to plot absolute loading values (TRUE) or signed loading values (FALSE). Default is FALSE.
#' @param col a vector of colors for the points on the plot. Default is c("gold", "white", "royalblue").
#' @param colMidpoint a numeric value indicating the midpoint for the color gradient. Default is 0.
#' @param shape an integer indicating the shape of the points on the plot. Default is 21.
#' @param shapeSizeRange a numeric vector of length 2 indicating the range of point sizes. Default is c(10, 10).
#' @param legendPosition a character string indicating the position of the legend. Default is "top".
#' @param legendLabSize an integer indicating the size of the legend labels. Default is 10.
#' @param legendIconSize a numeric value indicating the size of the legend icons. Default is 3.
#' @param xlim a numeric vector of length 2 indicating the range of the x-axis. Default is NULL.
#' @param ylim a numeric vector of length 2 indicating the range of the y-axis. Default is NULL.
#' @param labSize an integer indicating the size of the point labels. Default is 2.
#' @param labhjust a numeric value indicating the horizontal justification of the point labels. Default is 1.5.
#' @param labvjust a numeric value indicating the vertical justification of the point labels. Default is 0.
#' @param drawConnectors a logical value indicating whether to draw connectors between points. Default is TRUE.
#' @param positionConnectors a character string indicating the position of the connectors. Default is "right".
#' @param widthConnectors a numeric value indicating the width of the connectors. Default is 0.5.
#' @param typeConnectors a character string indicating the type of the connectors. Default is "closed".
#' @param endsConnectors a character string indicating the ends of the connectors. Default is "first".
#' @param lengthConnectors a numeric value indicating the length of the connectors. Default is unit(0.01, "npc").
#' @param colConnectors a character string indicating the color of the connectors. Default is "grey50".
#' @param xlab a character string indicating the label for the x-axis. Default is "Principal component".
#' @param xlabAngle a numeric value indicating the angle of the x-axis label. Default is 0.
#' @param xlabhjust a numeric value indicating the horizontal justification of the x-axis label. Default is 0.5.
#' @param xlabvjust a numeric value indicating the vertical justification of the x-axis label. Default is 0.5.
#' @param ylab a character string indicating the label for the y-axis. Default is "Component loading".
#' @param ylabAngle a numeric value indicating the angle of the y-axis label. Default is 0.
#' @param ylabhjust a numeric value indicating the horizontal justification of the y-axis label. Default is 0.5.
#' @param ylabvjust a numeric value indicating the vertical justification of the y-axis label. Default is 0.5.
#' @param axisLabSize an integer indicating the size of the axis labels. Default is 16.
#' @param title a character string indicating the title of the plot. Default is an empty string.
#' @param subtitle a character string indicating the subtitle of the plot. Default is an empty string.
#' @param caption a character string indicating the caption of the plot. Default is an empty string.
#' @param titleLabSize an integer indicating the size of the title label. Default is 16.
#' @param subtitleLabSize an integer indicating the size of the subtitle label. Default is 12.
#' @param captionLabSize an integer indicating the size of the caption label. Default is 12.
#' @param hline a numeric vector indicating the horizontal lines to be plotted. Default is c(0).
#' @param hlineType a character string indicating the type of the horizontal lines. Default is "longdash".
#' @param hlineCol a character string indicating the color of the horizontal lines. Default is "black".
#' @param hlineWidth a numeric value indicating the width of the horizontal lines. Default is 0.4.
#' @param vline a numeric vector indicating the vertical lines to be plotted. Default is NULL.
#' @param vlineType a character string indicating the type of the vertical lines. Default is "longdash".
#' @param vlineCol a character string indicating the color of the vertical lines. Default is "black".
#' @param vlineWidth a numeric value indicating the width of the vertical lines. Default is 0.4.
#' @param gridlines.major a logical value indicating whether to show major gridlines. Default is TRUE.
#' @param gridlines.minor a logical value indicating whether to show minor gridlines. Default is TRUE.
#' @param borderWidth a numeric value indicating the width of the border. Default is 0.8.
#' @param borderColour a character string indicating the color of the border. Default is "black".
#' @param plot a logical value indicating whether to plot the graph. Default is TRUE.
#' @param export a logical value indicating whether to export the plot. Default is TRUE.
#' @return a scatter plot of the loading values for each gene/protein for the specified principal components in a PCA object.
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_hline geom_vline
#'   theme_bw theme labs guides guide_legend scale_size
#'   scale_fill_continuous scale_fill_gradient2 xlim ylim
#'   element_rect element_line element_text element_blank arrow
#' @importFrom ggrepel geom_text_repel
#' @export
prot.plot_loadings <- function (pcaobj, components = PCAtools::getComponents(pcaobj, seq_len(5)),
                                rangeRetain = 0.05, absolute = FALSE, col = c("gold",
                                                                              "white", "royalblue"), colMidpoint = 0, shape = 21,
                                shapeSizeRange = c(10, 10), legendPosition = "top",
                                legendLabSize = 10, legendIconSize = 3, xlim = NULL, ylim = NULL,
                                labSize = 2, labhjust = 1.5, labvjust = 0, drawConnectors = TRUE,
                                positionConnectors = "right", widthConnectors = 0.5,
                                typeConnectors = "closed", endsConnectors = "first",
                                lengthConnectors = unit(0.01, "npc"), colConnectors = "grey50",
                                xlab = "Principal component", xlabAngle = 0, xlabhjust = 0.5,
                                xlabvjust = 0.5, ylab = "Component loading", ylabAngle = 0,
                                ylabhjust = 0.5, ylabvjust = 0.5, axisLabSize = 16, title = "",
                                subtitle = "", caption = "", titleLabSize = 16,
                                subtitleLabSize = 12, captionLabSize = 12, hline = c(0),
                                hlineType = "longdash", hlineCol = "black", hlineWidth = 0.4,
                                vline = NULL, vlineType = "longdash", vlineCol = "black",
                                vlineWidth = 0.4, gridlines.major = TRUE, gridlines.minor = TRUE,
                                borderWidth = 0.8, borderColour = "black", plot = TRUE, export = FALSE)
{
  x <- pcaobj$loadings[, components, drop = FALSE]
  retain <- c()
  for (i in seq_along(components)) {
    offset <- (max(x[, i]) - min(x[, i])) * rangeRetain
    uppercutoff <- max(x[, i]) - offset
    lowercutoff <- min(x[, i]) + offset
    retain <- unique(c(retain, which(x[, i] >= uppercutoff),
                       which(x[, i] <= lowercutoff)))
  }
  message("-- variables retained:")
  message(paste0(rownames(x)[retain], collapse = ", "))
  x <- x[retain, , drop = FALSE]
  PC <- Loading <- NULL
  x <- data.frame(rownames(x), x[, components, drop = FALSE], check.names = FALSE)
  colnames(x)[1] <- "var"
  plotobj <- reshape::melt(x, id = "var")
  colnames(plotobj) <- c("var", "PC", "Loading")
  if (absolute == TRUE) {
    plotobj$Loading <- abs(plotobj$Loading)
  }
  else if (absolute == FALSE) {
    plotobj$Loading <- plotobj$Loading
  }
  th <- theme_bw(base_size = 24) + theme(legend.background = element_rect(),
                                         plot.title = element_text(angle = 0, size = titleLabSize,
                                                                   face = "bold", vjust = 1), plot.subtitle = element_text(angle = 0,
                                                                                                                           size = subtitleLabSize, face = "plain", vjust = 1),
                                         plot.caption = element_text(angle = 0, size = captionLabSize,
                                                                     face = "plain", vjust = 1), axis.line = element_line(size = 1.5,
                                                                                                                          colour = "black"), axis.text.x = element_text(angle = xlabAngle,
                                                                                                                                                                        size = axisLabSize, hjust = xlabhjust, vjust = xlabvjust),
                                         axis.text.y = element_text(angle = ylabAngle, size = axisLabSize,
                                                                    hjust = ylabhjust, vjust = ylabvjust), axis.title = element_text(size = axisLabSize),
                                         legend.position = legendPosition, legend.direction = "horizontal",
                                         legend.box = "horizontal", legend.key = element_blank(),
                                         legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = legendLabSize),
                                         title = element_text(size = legendLabSize), legend.title = element_blank())

  p <- ggplot(plotobj, aes(x = PC, y = Loading, size = Loading,
                           fill = Loading)) + th + geom_point(shape = shape) + labs(title = title,
                                                                                    subtitle = subtitle, caption = caption) + labs(x = xlab,
                                                                                                                                   y = ylab, size = ylab, fill = ylab) + guides(fill = guide_legend(),
                                                                                                                                                                                size = guide_legend(), colour = guide_legend(override.aes = list(size = legendIconSize)))
  p <- p + scale_size(range = shapeSizeRange)
  if (length(col) == 2) {
    p <- p + scale_fill_continuous(low = col[1], high = col[2])
  }
  else if (length(col) == 3) {
    p <- p + scale_fill_gradient2(low = col[1], mid = col[2],
                                  high = col[3], midpoint = colMidpoint, space = "Lab")
  }
  if (!is.null(xlim)) {
    p <- p + xlim(xlim[1], xlim[2])
  }
  if (!is.null(ylim)) {
    p <- p + ylim(ylim[1], ylim[2])
  }
  if (!is.null(vline)) {
    p <- p + geom_vline(xintercept = vline, linetype = vlineType,
                        colour = vlineCol, size = vlineWidth)
  }
  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = hlineType,
                        colour = hlineCol, size = hlineWidth)
  }
  p <- p + theme(panel.border = element_rect(colour = borderColour,
                                             fill = NA, size = borderWidth))
  if (gridlines.major == TRUE) {
    p <- p + theme(panel.grid.major = element_line())
  }
  else {
    p <- p + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor == TRUE) {
    p <- p + theme(panel.grid.minor = element_line())
  }
  else {
    p <- p + theme(panel.grid.minor = element_blank())
  }
  if (drawConnectors == TRUE) {
    p <- p + ggrepel::geom_text_repel(data = plotobj, aes(label = as.character(var)),
                                      size = labSize, nudge_x = ifelse(positionConnectors ==
                                                                         "left", -0.75, ifelse(positionConnectors ==
                                                                                                 "right", 0.75, 0)), direction = "y",
                                      segment.color = colConnectors, segment.size = widthConnectors,
                                      arrow = ggplot2::arrow(length = lengthConnectors, type = typeConnectors,
                                                    ends = endsConnectors), show.legend = FALSE,
                                      hjust = labhjust, vjust = labvjust, max.overlaps = 25)
  }
  else if (drawConnectors == FALSE) {
    p <- p + geom_text(data = plotobj, aes(label = as.character(var)),
                       size = labSize, check_overlap = TRUE, show.legend = FALSE,
                       hjust = labhjust, vjust = labvjust)
  }
  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
    message(paste0("Exporting loadings plot to:\n\"", getwd(), "/Plots/LoadingsPlot.pdf\" and \".../LoadingsPlot.png\""))

    grDevices::pdf("Plots/LoadingsPlot.pdf")
    print(p)
    grDevices::dev.off()

    grDevices::png("Plots/LoadingsPlot.png",
                   width = 6, height = 6, units = 'in', res = 300)
    print(p)
    grDevices::dev.off()
  }
  if (plot == TRUE) {
    return(p)
  }
  else if (plot == FALSE) {
    p
  }
}



#' Plot pathway enrichment
#'
#' Plots the enrichment of pathways in a given dataset.
#'
#' @param enrichset A data frame containing the enrichment results.
#' @param title The title of the plot.
#' @param subtitle The subtitle of the plot.
#' @param plot Logical; whether to return the plot or not.
#' @param export Logical; whether to export the plot as PDF an PNG files
#' @param kegg Logical; whether the enrichment results are from KEGG or not.
#'
#' @return A ggplot object.
#'
#'
#' @export
#'
#' @importFrom forcats fct_reorder
#' @importFrom ggplot2 ggplot aes geom_segment geom_point
#'   scale_color_viridis_c scale_size_area theme_bw xlab ylab ggtitle guide_colorbar
#' @importFrom grDevices pdf png
#'
#' @seealso
#' \code{\link{pathway_enrich}}
prot.plot_enrichment <- rna.plot_enrichment <- function(enrichset,
                                 title = "Differentially enriched pathways",
                                 subtitle = "",
                                 plot = TRUE,
                                 export = FALSE,
                                 kegg = TRUE) {
  if ( kegg == TRUE ) {
    enrichset@result$Description <- gsub(" - .*", "", enrichset@result$Description)
  }

  p <- ggplot(enrichset, showCategory = 30,
              aes(richFactor, fct_reorder(Description, richFactor))) +
    geom_segment(aes(xend = 0, yend = Description)) +
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_viridis_c(guide = guide_colorbar(reverse = TRUE)) +
    scale_size_area(max_size = 10) +
    theme_bw(base_size = 16) +
    xlab("Rich factor") +
    ylab(NULL) +
    ggtitle(title, subtitle = subtitle)

  if (export == TRUE) {
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
    message(
      paste0(
        "Exporting pathway enrichment plot to:\n\"",
        getwd(),
        "/Plots/EnrichmentPlot_",
        if_else(
          grepl("Diff", title, ignore.case = T),
          "Diff", if_else(
            grepl("Up", title, ignore.case = T),
            "Up", if_else(
              grepl("Down", title, ignore.case = T),
              "Down", ""))),
        subtitle,
        if_else(kegg==TRUE, "-KEGG", "-custom"),
        ".pdf\" and \".../EnrichmentPlot_",
        if_else(
          grepl("Diff", title, ignore.case = T),
          "Diff", if_else(
            grepl("Up", title, ignore.case = T),
            "Up", if_else(
              grepl("Down", title, ignore.case = T),
              "Down", ""))),
        subtitle,
        if_else(kegg==TRUE, "-KEGG", "-custom"),
        ".png\""
      )
    )
    w <- 6+max(str_count(enrichset$Description))/6
    h <- 5+length(enrichset$Description)/6
    grDevices::pdf(paste0("Plots/EnrichmentPlot_",
                          if_else(
                            grepl("Diff", title, ignore.case = T),
                            "Diff_", if_else(
                              grepl("Up", title, ignore.case = T),
                              "Up_", if_else(
                                grepl("Down", title, ignore.case = T),
                                "Down_", ""))),
                          subtitle,
                          if_else(kegg==TRUE, "-KEGG", "-custom"),
                          ".pdf"),
                   width = w,
                   height = h)
    print(p)
    grDevices::dev.off()

    grDevices::png(
      paste0("Plots/EnrichmentPlot_",
             if_else(
               grepl("Diff", title, ignore.case = T),
               "Diff_", if_else(
                 grepl("Up", title, ignore.case = T),
                 "Up_", if_else(
                   grepl("Down", title, ignore.case = T),
                   "Down_", ""))),
             subtitle,
             if_else(kegg==TRUE, "-KEGG", "-custom"),
             ".png"),
      width = w,
      height = h,
      units = 'in',
      res = 300
    )
    print(p)
    grDevices::dev.off()
  }
  if (plot) {
    print(p)
  } else {
    p
  }
}

#' @title Plot an UpSet enrichment plot
#'
#' @description Plot an UpSet chart from an enrichment set
#'
#' @param enrichset An enrichment set object
#' @param order.by A character string indicating how to order the sets
#' @param point.size Numeric indicating the size of the points
#' @param line.size Numeric indicating the size of the lines
#' @param text.scale Numeric indicating the size of the text
#' @param ... Further arguments to pass to the upset function
#'
#' @return An UpSetR plot
#'
#' @details for further arguments, see `?UpSetR::upset`
#'
#' @details For more information about upset plots, see https://upset.app/
#'
#' @export
#'
prot.plot_upset <- function(enrichset, order.by = "freq", point.size = 3,
                            line.size = 1, text.scale = c(2, 2, 2, 2, 2, 1.5), ...)
{
  # Get the number of pathways (sets)
  num_pathways <- nrow(enrichset@result)

  # Dynamically adjust text scale based on the number of pathways
  if (num_pathways > 30) {
    text.scale <- c(1.2, 1.2, 1.2, 1.2, 1, 1)  # Smallest text scale for very large sets
  } else if (num_pathways > 20) {
    text.scale <- c(1.5, 1.5, 1.5, 1.5, 1.2, 1)  # Smaller text scale for large sets
  } else if (num_pathways > 10) {
    text.scale <- c(1.8, 1.8, 1.8, 1.8, 1.5, 1.2)  # Moderate scaling for medium-sized sets
  } else {
    text.scale <- c(2, 2, 2, 2, 2, 1.5)  # Larger text scale for smaller sets
  }

  enrichset@result$Description <- gsub(" - .*", "", enrichset@result$Description)

  # Prepare the upset list (convert gene IDs to list format)
  upsetlist <- enrichset$geneID %>%
    str_replace_all("/", ", ") %>%
    strsplit(", ")

  # Assign pathway descriptions as names for each list element
  names(upsetlist) <- enrichset$Description

  # Generate the UpSet plot using the dynamically adjusted text scaling
  UpSetR::upset(UpSetR::fromList(upsetlist),
                order.by = order.by,
                nsets = length(upsetlist),
                point.size = point.size,
                line.size = line.size,
                text.scale = text.scale,
                ...)
}



#' Plot bar plots for protein abundance or fold change
#'
#' This function generates bar plots for protein abundance or fold change,
#' normalized by either a reference protein or centered to the mean.
#' The user can also plot the contrast of two conditions.
#'
#' @param dep A SummarizedExperiment object containing the data.
#' @param proteins A vector of protein names or IDs.
#' @param combine Logical value indicating whether to combine the data from different proteins in one plot.
#' @param type Type of plot: "abundance", "ReferenceProt", "ReferenceCondition",
#'   "centered", or "contrast".
#' @param ref.prot ID or name of the reference protein for normalization
#'   when type is "ReferenceProt".
#' @param ref.condition Name of the reference condition by which all average
#'   abundance values are divided when type is "ReferenceCondition".
#' @param contrast Name of the contrast for plotting when type is "contrast".
#' @param col.id Name of the column with protein IDs in \code{dep}.
#' @param name_table A data frame containing the mapping of protein names
#'   to reference names.
#' @param match.id The name of a column in \code{name_table} to find matches
#'   with \code{col.id}.
#' @param match.name The name of a column in \code{name_table} to assign new
#'   names based on matches of \code{col.id} and \code{match.id}.
#' @param convert_name Logical value indicating whether to convert the protein
#'   names to reference names based on \code{name_table} (e.g., gene names).
#' @param shape.size Numeric value defining the (replicate) symbol size in the plot.
#' @param y.lim Limits of the y axis (numeric vector of length 2). If NULL,
#'   will be inferred from data (unless `combine=TRUE`, then each chunk is
#'   plotted without that inference).
#' @param basesize Base font size for the plot theme.
#' @param plot Whether to immediately print/return the plot(s).
#' @param export Logical; whether to export the plot(s) as PDF and PNG files.
#' @param export.nm Name of the output PDF and PNG files if \code{export = TRUE}.
#' @param width Numeric value defining the width of the exported plot (inches).
#' @param height Numeric value defining the height of the exported plot (inches).
#'
#' @return (Invisibly) A ggplot object (or list of ggplot objects, if multiple chunks).
#'         If \code{plot=TRUE}, the plot(s) are also shown in the Plots pane.
#'
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom tidyr gather spread
#' @importFrom stats qnorm
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stringr str_split str_count
#' @importFrom ggplot2 ggplot aes geom_hline geom_col geom_errorbar
#'   facet_wrap scale_fill_manual labs theme_bw coord_cartesian
#'   scale_x_discrete theme element_text position_dodge
#' @import dplyr
#' @export
prot.plot_bar <- function(dep,
                          proteins,
                          combine = FALSE,
                          type = c("contrast",
                                   "centered",
                                   "ReferenceProt",
                                   "ReferenceCondition",
                                   "abundance"),
                          ref.prot = NULL,
                          ref.condition = NULL,
                          contrast = NULL,
                          col.id = "ID",
                          name_table = NULL,
                          match.id = "Accession",
                          match.name = "Name",
                          convert_name = FALSE,
                          shape.size = 2.5,
                          y.lim = NULL,
                          basesize = 12,
                          plot = TRUE,
                          export = FALSE,
                          export.nm = NULL,
                          width = NULL,
                          height = NULL)
{
  # If name_table is defined, make name_table[[match.name]] unique by appending "(A)", "(B)", etc. to duplicates
  # if (!is.null(name_table)) {
  #   name_table[[match.name]] <- make.unique(name_table[[match.name]])
  # }
  # Custom function to make names unique by appending (A), (B), (C) ...
  make_unique_letters <- function(x) {
    # 'ave' will apply the given function for each group of identical values in x
    ave(x, x, FUN = function(z) {
      if (length(z) == 1) {
        # only one occurrence -> nothing to change
        return(z)
      } else {
        # for duplicates, keep the first occurrence the same,
        # append (A) to the second, (B) to the third, etc.
        for (i in seq_along(z)) {
          if (i > 1) {
            suffix <- LETTERS[i - 1]            # i=2 -> 'A', i=3 -> 'B', etc.
            if (is.na(suffix)) suffix <- i - 1  # fallback if >26 duplicates
            z[i] <- paste0(z[i], " (", suffix, ")")
          }
        }
        z
      }
    })
  }

  ## --------------------------------------------------------------------------
  ## 1) HELPER FUNCTIONS
  ## --------------------------------------------------------------------------

  # Helper A: Split an integer vector into chunks of size `chunk_size`.
  chunk_indices <- function(n, chunk_size = 8) {
    idx <- seq_len(n)
    group_id <- ceiling(idx / chunk_size)  # e.g. 1,1,1,1,2,2,2,2,...
    split(idx, group_id)
  }

  # Helper B: Convert protein names (vector) using name_table
  convert_protein_names <- function(vec, name_table, from_col, to_col) {
    name_table[[to_col]][ match(vec, name_table[[from_col]]) ]
  }

  # Helper C: Palettes for condition
  condition_palette <- function(n_condition) {
    if (n_condition <= 8) {
      RColorBrewer::brewer.pal(n_condition, "Dark2")
    } else if (n_condition <= 12) {
      RColorBrewer::brewer.pal(n_condition, "Set3")
    } else {
      c(
        "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black",
        "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
        "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1",
        "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3",
        "darkorange4", "brown"
      )
    }
  }

  # Helper D: Palettes for replicate
  replicate_palette <- function(n_reps) {
    if (n_reps <= 2) {
      # Exactly 1 or 2 colors
      c("#dad9d9", "#878787")[1:n_reps]
    } else if (n_reps <= 8) {
      RColorBrewer::brewer.pal(n_reps, "Greys")
    } else {
      # Gradient
      grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(8, "Greys")
      )(n_reps)
    }
  }

  # Helper E: Summarize replicate values into mean, CI.L, CI.R
  summarize_with_ci <- function(df_reps, value_col = "val") {
    df_reps %>%
      dplyr::group_by(condition, rowname) %>%
      dplyr::summarize(
        mean = mean(.data[[value_col]], na.rm = TRUE),
        sd   = sd(.data[[value_col]], na.rm = TRUE),
        n    = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        error = stats::qnorm(0.975) * .data$sd / sqrt(.data$n),
        CI.L  = .data$mean - .data$error,
        CI.R  = .data$mean + .data$error
      )
  }

  # Helper F: Summarize for contrast data
  #   e.g., returns a data frame with columns: "name"/"ID", "diff", "CI.L", "CI.R", "contrast"
  gather_contrasts <- function(row_data_sub, contrast = NULL, convert_name, name_table, match.id, match.name) {
    # Figure out which columns to select
    # if user gave a specific contrast, we pick e.g. "XX_diff", "XX_CI.L", "XX_CI.R".
    # otherwise we pick all columns that end in _diff, _CI.L, _CI.R
    if (!is.null(contrast)) {
      needed_cols <- c("name", paste0(contrast, "_diff"), paste0(contrast, "_CI.L"), paste0(contrast, "_CI.R"))
    } else {
      # pick name plus any that match: ends_with(_diff,_CI.L,_CI.R)
      # use tidyselect approach
      needed_cols <- c("name", grep("_diff$|_CI\\.L$|_CI\\.R$", colnames(row_data_sub), value = TRUE))
    }
    df <- row_data_sub[, needed_cols, drop = FALSE]
    df <- as.data.frame(df, check.names = FALSE)

    # pivot to long form
    df_long <- tidyr::gather(df, var, val, -name)
    df_long$contrast <- gsub("_diff|_CI\\.L|_CI\\.R", "", df_long$var)
    df_long$var      <- gsub(".*_", "", df_long$var)   # e.g. "diff" or "CI.L" or "CI.R"
    df_wide <- tidyr::spread(df_long, var, val)        # columns: name, contrast, diff, CI.L, CI.R

    if (convert_name) {
      #message("Converting protein names to reference names...")
      # add debug breakpoint
      # Example usage inside your if-statement:
      if (!is.null(name_table)) {
        # reduce the name_table to rows where df_reps$rowname matches name_table$match.id
        name_table <- name_table[match(df_reps$rowname, name_table[[match.id]]), ]
        # Suppose 'match.name' is the column whose values you need to deduplicate
        name_table[[match.name]] <- make_unique_letters(name_table[[match.name]])
      }
      else {
        stop("Please provide a 'name_table' when 'convert_name' is TRUE.")
      }
      # Convert all rownames:
      converted_names <- convert_protein_names(
        df_reps$rowname,
        name_table,
        from_col = match.id,
        to_col   = match.name
      )

      # Fallback: keep the original row name if 'converted_names' is NA.
      df_reps$rowname <- ifelse(
        is.na(converted_names),
        df_reps$rowname,   # fallback to the old name
        converted_names    # use the mapped name
      )
    }
    df_wide
  }

  # Helper G: Perform the data transformation for each of the "non-contrast" types
  # Returns a data.frame with columns rowname, condition, replicate, val, etc.
  transform_data_noncontrast <- function(dep_subset,
                                         type,
                                         ref.prot,
                                         ref.condition,
                                         row_data_full) {
    # We'll compute a matrix of values for each row (protein) x column (replicate).
    mat <- SummarizedExperiment::assay(dep_subset)

    if (type == "centered") {
      # subtract rowMeans
      means <- rowMeans(mat, na.rm = TRUE)
      mat   <- sweep(mat, 1, means, FUN = "-")
    } else if (type == "abundance") {
      # do nothing
      # mat stays as it is
    } else if (type == "ReferenceProt") {
      # find the reference's mean
      # note: row_data_full is the rowData for the *full* dep, so we can find the reference
      if (is.null(ref.prot)) {
        stop("Please provide ref.prot when type='ReferenceProt'.")
      }
      # figure out if the user gave ID or 'name'
      # search row_data_full by "ID" or "name"
      ref_idx <- grep(ref.prot, row_data_full[,"name"])
      col_used <- "name"
      if (length(ref_idx) == 0) {
        # try "ID"
        ref_idx <- grep(ref.prot, row_data_full[,"ID"])
        col_used <- "ID"
      }
      if (length(ref_idx) == 0) {
        stop("'ref.prot' not found among row_data names or IDs.")
      }
      # get reference name for labeling
      reference_name <- row_data_full[ref_idx, "name"]

      # compute reference mean
      ref_mean <- mean(SummarizedExperiment::assay(dep)[ref_idx, ], na.rm = TRUE)
      mat <- mat - ref_mean
    } else if (type == "ReferenceCondition") {
      # find columns that match the ref.condition
      if (is.null(ref.condition)) {
        stop("Please provide 'ref.condition' when type='ReferenceCondition'.")
      }
      cond_vec <- SummarizedExperiment::colData(dep_subset)$condition
      # Actually let's look in the *full* colData for a moment, or just do the naive approach:
      # We use `dep` to find columns whose condition == ref.condition, then compute mean(2^values)
      # This is how the original code did it:
      dep_assay <- SummarizedExperiment::assay(dep)
      col_cond  <- SummarizedExperiment::colData(dep)$condition
      if (!ref.condition %in% col_cond) {
        stop("ref.condition is not in colData(dep)$condition.")
      }
      # compute reference mean in linear space (2^assay)
      ref_mean <- mean(2^dep_assay[proteins,
                                   which(col_cond == ref.condition)],
                       na.rm = TRUE)
      # now transform mat: 2^mat / ref_mean
      mat <- 2^mat / ref_mean
    }

    # Now gather into long form
    df_reps <- data.frame(mat, check.names = FALSE) %>%
      tibble::rownames_to_column("rowname") %>%
      tidyr::gather("ID", "val", -rowname)

    # Merge with colData
    coldata_sub <- data.frame(SummarizedExperiment::colData(dep_subset),
                              check.names = FALSE) %>%
      tibble::rownames_to_column("ID")
    df_reps <- dplyr::left_join(df_reps, coldata_sub, by = "ID")
    df_reps$replicate <- as.factor(df_reps$replicate)

    df_reps
  }

  # Helper H: Plot a single chunk (non-contrast) or (contrast)
  make_plot <- function(df_main,
                        type,
                        combine,
                        y.lim,
                        basesize,
                        shape.size,
                        n_conditions,
                        n_reps,
                        reference_name = NULL) {

    library(ggplot2)
    library(ggnewscale)

    if (type == "contrast") {
      # df_main is presumably the wide data with columns: name (or ID), contrast, diff, CI.L, CI.R
      # check if combine or not
      if (!combine) {
        p <- ggplot(df_main, aes(x = contrast, y = diff, fill = contrast)) +
          geom_hline(yintercept = 0) +
          geom_col(colour = "black", position = position_dodge(width = 0.9)) +
          geom_errorbar(aes(ymin = CI.L, ymax = CI.R),
                        width = 0.3, position = position_dodge(width = 0.9)) +
          facet_wrap(~ name, scales = "fixed") +
          scale_fill_manual(values = condition_palette(n_conditions)) +
          labs(
            x = "Contrasts",
            y = expression(log[2] ~ "Fold change" ~ "(\u00B195% CI)")
          ) +
          theme_bw(base_size = basesize)
      } else {
        # combine=TRUE => x-axis is "name" or "ID"
        p <- ggplot(df_main, aes(x = name, y = diff, fill = contrast)) +
          geom_hline(yintercept = 0) +
          geom_col(colour = "black", position = position_dodge(width = 0.9)) +
          geom_errorbar(aes(ymin = CI.L, ymax = CI.R),
                        width = 0.3, position = position_dodge(width = 0.9)) +
          scale_fill_manual(values = condition_palette(n_conditions)) +
          labs(
            x = "Protein",
            y = expression(log[2] ~ "Fold change" ~ "(\u00B195% CI)")
          ) +
          theme_bw(base_size = basesize)
      }
      if (!is.null(y.lim)) {
        p <- p + coord_cartesian(ylim = y.lim)
      }
      return(p+scale_x_discrete(labels = function(x) sub("\\s*\\([^()]*\\)$", "", x)) )
    }

    # Otherwise, we have "abundance", "centered", "ReferenceProt", "ReferenceCondition".
    # Typically we do a bar for the mean plus error bars, plus replicate points.
    # We already summarized with CI if needed outside, or we do it now.
    # We'll see if df_main is the replicate data or the summarized data.
    # Usually the original code did replicate data in "df_reps" and summarized in "df".

    # Let's do that summarizing here:
    df_summ <- summarize_with_ci(df_main, value_col = "val")

    # Depending on combine, we do facet_wrap or side-by-side bar
    if (!combine) {
      # x=condition, facet=protein
      p <- ggplot(df_summ, aes(x = condition, y = mean, fill = condition)) +
        geom_hline(yintercept = 0) +
        geom_col(colour = "black") +
        scale_fill_manual(values = condition_palette(n_conditions)) +
        ggnewscale::new_scale_fill() +
        geom_point(
          data = df_main,
          aes(x = condition, y = val, fill = replicate),
          shape = 23, size = shape.size, color = "black",
          position = position_dodge(width = 0.3)
        ) +
        scale_fill_manual(values = replicate_palette(n_reps)) +
        geom_errorbar(aes(ymin = CI.L, ymax = CI.R), width = 0.3) +
        facet_wrap(~ rowname) +
        theme(
          axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1)
          )+
        theme_DEP2(basesize = basesize)

      # Y-axis label depends on type
      if (type == "centered") {
        p <- p + labs(x = "Baits", y = expression(log[2] ~ "Centered intensity (\u00B195% CI)"))
      } else if (type == "abundance") {
        p <- p + labs(x = "Baits", y = expression(log[2] ~ "Intensity (\u00B195% CI)"))
      } else if (type == "ReferenceCondition") {
        p <- p + labs(x = "Baits", y = "Normalized abundance")
      } else if (type == "ReferenceProt") {
        p <- p + labs(
          x = "Baits",
          y = bquote(log[2](intensity) - log[2](intensity[.(reference_name)]) ~ "(\u00B195% CI)")
        )
      }

    } else {
      # combine=TRUE => x=protein, fill=condition => side-by-side bars
      p <- ggplot(df_summ, aes(x = rowname, y = mean, fill = condition)) +
        geom_hline(yintercept = 0) +
        geom_col(colour = "black", position = position_dodge()) +
        scale_fill_manual(values = condition_palette(n_conditions)) +
        geom_errorbar(aes(ymin = CI.L, ymax = CI.R),
                      width = 0.3, position = position_dodge(width = 0.9)) +
        theme(
          axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1)
        )+
        theme_DEP2(basesize = basesize)

      # Y-axis label depends on type again
      if (type == "centered") {
        p <- p + labs(x = "Baits", y = expression(log[2] ~ "Centered intensity (\u00B195% CI)"))
      } else if (type == "abundance") {
        p <- p + labs(x = "Baits", y = expression(log[2] ~ "Intensity (\u00B195% CI)"))
      } else if (type == "ReferenceCondition") {
        p <- p + labs(x = "Baits", y = "Normalized abundance")
      } else if (type == "ReferenceProt") {
        p <- p + labs(
          x = "Baits",
          y = bquote(log[2](intensity) - log[2](intensity[.(reference_name)]) ~ "(\u00B195% CI)")
        )
      }
    }
    if (!is.null(y.lim)) {
      p <- p + coord_cartesian(ylim = y.lim)
    }
    p + scale_x_discrete(labels = function(x) sub("\\s*\\([^()]*\\)$", "", x))
  }

  ## --------------------------------------------------------------------------
  ## 2) MAIN FUNCTION BODY
  ## --------------------------------------------------------------------------

  # Match type argument
  type <- match.arg(type)

  # Basic checks
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.character(proteins),
    is.logical(plot),
    length(plot) == 1
  )

  # Ensure rownames of `dep` match col.id
  row_data <- SummarizedExperiment::rowData(dep, use.names = FALSE)
  if (!col.id %in% colnames(row_data)) {
    stop("The specified 'col.id' (", col.id, ") is not a column in rowData(dep).")
  }
  rownames(dep) <- row_data[[col.id]]
  # replace any NA rownames with "0"
  rownames(dep)[is.na(rownames(dep))] <- "0"

  # Additional checks
  cd <- SummarizedExperiment::colData(dep)
  names(cd)[names(cd) == "ID"] <- "old_ID" # e.g. rename
  SummarizedExperiment::colData(dep) <- cd
  if (any(!c("label", "condition", "replicate") %in% colnames(cd))) {
    stop("'dep' must contain columns 'label', 'condition', 'replicate' in colData.\n",
         "Make sure you used make_se() or equivalent to create 'dep'.")
  }
  # If type=contrast, check for columns ending in _diff
  if (type == "contrast") {
    any_diff_cols <- any(grepl("_diff$", colnames(row_data)))
    if (!any_diff_cols) {
      stop("No columns named '[contrast]_diff' or similar found.\n",
           "Please run test_diff() or set type != 'contrast'.")
    }
  }
  # Must have "name" column in row_data
  if (!"name" %in% colnames(row_data)) {
    stop("Column 'name' not found in rowData(dep). Make sure you used the standard pipeline (make_se).")
  }

  # Check if user-supplied proteins are actually present
  all_proteins_in_data <- rownames(dep)
  not_found <- setdiff(proteins, all_proteins_in_data)
  if (length(not_found) == length(proteins)) {
    # none found
    stop("None of the specified proteins were found in the column '", col.id, "'.")
  } else if (length(not_found) > 0) {
    warning("These proteins were not found and will be ignored: ", paste(not_found, collapse = ", "))
    proteins <- setdiff(proteins, not_found)
  }

  # Create chunked subsets
  if (combine) {
    subset_list <- list(dep[proteins, ])
  } else {
    idx_list    <- chunk_indices(length(proteins), chunk_size = 8)
    subset_list <- lapply(idx_list, function(idx) {
      dep[proteins[idx], ]
    })
  }

  ## Optional: Pre-compute y.lim if combine=FALSE and user didn't specify it,
  ## because the original code tries to guess a global y.lim from all proteins
  ## in some types.
  ## We replicate the original logic only if y.lim is NULL, type != "contrast",
  ## combine=FALSE.
  if (!combine && is.null(y.lim) && type != "contrast") {
    # We'll do the transformation for the entire set of proteins, then
    # compute min & max from that.
    # The original code tried to do it chunk by chunk, but we can do it once:
    all_dep_subset <- dep[proteins, ]
    df_all <- transform_data_noncontrast(all_dep_subset,
                                         type,
                                         ref.prot,
                                         ref.condition,
                                         row_data)
    # Range of replicate values
    rng_vals <- range(df_all$val, na.rm = TRUE)
    # If needed, we also gather CI from summarized data
    df_summ <- summarize_with_ci(df_all, "val")
    rng_ci  <- range(c(df_summ$CI.L, df_summ$CI.R), na.rm = TRUE)
    overall_min <- min(rng_vals[1], rng_ci[1])
    overall_max <- max(rng_vals[2], rng_ci[2])
    y.lim <- c(overall_min, overall_max)
  }

  # We'll store all resulting ggplots in a list
  plots <- vector("list", length(subset_list))

  n_conditions <- length(unique(cd$condition))
  n_reps       <- as.numeric(max(cd$replicate))

  for (i in seq_along(subset_list)) {
    dep_subset <- subset_list[[i]]

    if (type == "contrast") {
      # We'll gather the relevant rowData columns
      row_data_sub <- SummarizedExperiment::rowData(dep_subset, use.names = FALSE)
      df_contrast  <- gather_contrasts(
        row_data_sub,
        contrast     = contrast,
        convert_name = convert_name,
        name_table   = name_table,
        match.id     = match.id,
        match.name   = match.name
      )
      p <- make_plot(df_contrast,
                     type = "contrast",
                     combine = combine,
                     y.lim = y.lim,
                     basesize = basesize,
                     shape.size = shape.size,
                     n_conditions = n_conditions,
                     n_reps = n_reps)
      plots[[i]] <- p

    } else {
      # For the other types
      df_reps <- transform_data_noncontrast(dep_subset,
                                            type,
                                            ref.prot,
                                            ref.condition,
                                            row_data_full = row_data)
      # If needed, convert protein names
      if (convert_name) {
        message("Converting protein names to reference names...")
        if (!is.null(name_table)) {
          # reduce the name_table to rows where df_reps$rowname matches name_table$match.id
          name_table <- name_table[match(df_reps$rowname, name_table[[match.id]]), ]
          # Suppose 'match.name' is the column whose values you need to deduplicate
          name_table[[match.name]] <- make_unique_letters(name_table[[match.name]])
        }
        else {
          stop("Please provide a 'name_table' when 'convert_name' is TRUE.")
        }
        # Convert all rownames:
        converted_names <- convert_protein_names(
          df_reps$rowname,
          name_table,
          from_col = match.id,
          to_col   = match.name
        )

        # Fallback: keep the original row name if 'converted_names' is NA.
        df_reps$rowname <- ifelse(
          is.na(converted_names),
          df_reps$rowname,   # fallback to the old name
          converted_names    # use the mapped name
        )
      }
      # For referenceProt label, we might want to pass reference_name to the plot
      reference_name <- NULL
      if (type == "ReferenceProt" && !is.null(ref.prot)) {
        # Attempt to find the actual "name" from row_data
        ref_idx <- grep(ref.prot, row_data[,"name"])
        if (length(ref_idx) == 0) {
          ref_idx <- grep(ref.prot, row_data[,"ID"])
        }
        if (length(ref_idx)) {
          reference_name <- row_data[ref_idx[1], "name"]
        }
      }

      p <- make_plot(df_reps,
                     type = type,
                     combine = combine,
                     y.lim = y.lim,
                     basesize = basesize,
                     shape.size = shape.size,
                     n_conditions = n_conditions,
                     n_reps = n_reps,
                     reference_name = reference_name)
      plots[[i]] <- p
    }

    # Handle export if requested
    if (export) {
      # Compute file name
      if (length(proteins) <= 8) {
        chunk_suffix <- ""
      } else {
        chunk_suffix <- paste0("_", i)
      }
      if (!is.null(export.nm)) {
        export_name <- file.path("Plots", paste0(export.nm, chunk_suffix))
      } else {
        # fallback name
        export_name <- file.path("Plots",
                                 paste0("BarPlot_",
                                        paste(proteins, collapse = "_"),
                                        chunk_suffix))
      }
      dir.create("Plots", showWarnings = FALSE)
      message("Exporting bar plot(s) to:\n'",
              paste0(export_name, ".pdf' and '", export_name, ".png'"))

      w <- if (!is.null(width)) width else 8
      h <- if (!is.null(height)) height else 10

      grDevices::pdf(paste0(export_name, ".pdf"), width = w, height = h)
      print(plots[[i]])
      grDevices::dev.off()

      grDevices::png(paste0(export_name, ".png"),
                     width = w, height = h, units = "in", res = 300)
      print(plots[[i]])
      grDevices::dev.off()
    }

    # If plot=TRUE, print to the Plots pane
    if (plot) {
      print(plots[[i]])
    }
  } # end for-loop

  # If only one chunk, return that plot invisibly
  if (length(plots) == 1) {
    invisible(plots[[1]])
  } else {
    invisible(plots)
  }
}

#' Boxplot Intensity
#'
#' Generate boxplot for the intensity of each sample in a
#' SummarizedExperiment object.
#'
#' @param se A SummarizedExperiment object.
#' @param ... Additional SummarizedExperiment objects.
#' @param plot Logical. Should the plot be printed?
#' @param export Logical. Should the plot be exported as PDF and PNG files?
#' @return (Invisibly) a ggplot object.
#' @importFrom ggplot2 ggplot aes geom_boxplot coord_flip facet_wrap
#'   labs scale_fill_brewer scale_fill_manual
#' @export
prot.boxplot_intensity <- function (se, ..., plot = TRUE, export = FALSE)
{
  call <- match.call()
  call$export <- NULL
  call$plot <- NULL
  arglist <- lapply(call[-1], function(x) x)
  var.names <- vapply(arglist, deparse, character(1))
  arglist <- lapply(arglist, eval.parent, n = 2)
  names(arglist) <- var.names
  lapply(arglist, function(x) {
    if (any(!c("label", "ID", "condition",
               "replicate") %in% colnames(SummarizedExperiment::colData(x)))) {
      stop("'label', 'ID', 'condition' and/or 'replicate' ",
           "columns are not present in (one of) the input object(s)",
           "\nRun make_se() or make_se_parse() to obtain the required columns",
           call. = FALSE)
    }
  })
  gather_join <- function(se) {
    SummarizedExperiment::assay(se) %>% data.frame(check.names = FALSE) %>% gather(ID, val) %>%
      left_join(., data.frame(SummarizedExperiment::colData(se), check.names = FALSE), by = "ID")
  }
  df <- purrr::map_df(arglist, gather_join, .id = "var") %>%
    mutate(var = factor(var, levels = names(arglist)))
  p <- ggplot(df, aes(x = ID, y = val, fill = condition)) +
    ggplot2::geom_boxplot(notch = TRUE, na.rm = TRUE) +
    coord_flip() +
    facet_wrap( ~ var, ncol = 1) +
    labs(x = "", y = expression(log[2] ~ "Intensity"))


  if (length(se$condition) > 15) {
    n.cond <- length(se$condition)
    theme <- theme_DEP1(basesize = 10)
    theme[["axis.text.y"]]["size"] <- exp(-n.cond / 38) * 16
    p <- p + theme
  } else{
    theme <- theme_DEP1(basesize = 10)
    p <- p + theme
  }

  if (length(unique(se$condition)) <= 8) {
    p <- p + scale_fill_brewer(palette = "Dark2")
  } else if (length(unique(se$condition)) <= 12) {
    p <- p + scale_fill_brewer(palette = "Set3")
  } else if (length(unique(se$condition)) <= 25) {
    p <- p + scale_fill_manual(
      values = c(
        "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
        "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
        "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
        "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
        "green1", "yellow4", "yellow3", "darkorange4", "brown"
      )
    )
  } else {
    p <- p + scale_fill_manual(values = rainbow(length(unique(se$condition))))
  }

  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
    message(paste0("Exporting normalization plots to:\n\"", getwd(), "/Plots/BoxPlotIntensity.pdf\" and \".../BoxPlotIntensity.png\""))
    grDevices::png("Plots/BoxPlotIntensity.png",
                   width = 8, height = 10, units = 'in', res = 300)
    print(p)
    grDevices::dev.off()
    grDevices::pdf("Plots/BoxPlotIntensity.pdf",
                   width = 8, height = 10)
    print(p)
    grDevices::dev.off()
  }

  if (plot == TRUE){
    print(p)
  }
  invisible(p)
}

#' @title Theme DEP2
#' @description Modify the theme DEP1
#' @param ... Arguments passed to \code{\link{theme_DEP1}}
#' @return A modified theme
#' @export
theme_DEP2 <- function(...)
{
  theme <- theme_DEP1(...)
  theme$axis.text.x$angle <- 90
  theme$axis.text.x$hjust <- 1
  theme$axis.text.x$vjust <- 0.5
  return(theme)
}



#' Plot coverage of a SummarizedExperiment
#'
#' @param se A SummarizedExperiment
#' @param plot A logical indicating whether to return the coverage or not
#'
#' @return A data frame (if plot = FALSE) or a ggplot object
#'
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom SummarizedExperiment assay
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_grey labs
plot_coverage <- function (se, plot = TRUE)
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.logical(plot), length(plot) == 1)
  df <- SummarizedExperiment::assay(se) %>% data.frame() %>% tibble::rownames_to_column() %>%
    gather(ID, bin, -rowname) %>% mutate(bin = ifelse(is.na(bin),
                                                      0, 1))
  stat <- df %>% group_by(rowname) %>% summarize(sum = sum(bin))
  table <- table(stat$sum) %>% data.frame()
  p <- ggplot(table, aes(x = "all", y = Freq, fill = Var1)) +
    geom_col(col = "white") + scale_fill_grey(start = 0.8,
                                              end = 0.2) + labs(title = "Protein coverage", x = "",
                                                                y = "Number of proteins", fill = "Samples") + theme_DEP1()
  if (plot) {
    return(p)
  }
  else {
    df <- as.data.frame(table)
    colnames(df) <- c("samples", "proteins")
    return(df)
  }
}



#' @title Plot row standard deviations versus row means
#'
#' @description
#' This function creates a mean-sd plot for a given SummarizedExperiment object.
#'
#' @param x A SummarizedExperiment object.
#' @param ranks A logical value indicating whether to plot the rank of the mean instead of the mean.
#' @param xlab A character string for the x-axis label.
#' @param ylab A character string for the y-axis label.
#' @param pch A vector of plotting characters.
#' @param plot A logical value indicating whether to plot the mean-sd plot.
#' @param bins An integer indicating the number of bins in the plot.
#' @param ... Other arguments to be passed to the underlying plotting function.
#'
#' @return
#' If plot is TRUE, a mean-sd plot is returned. In any case returns a named list with five components: its elements px and py are the x- and y-coordinates of the individual data points in the plot; its first and second element are the x-coordinates and values of the running median estimator (the red line in the plot). Its element gg is the plot object (see examples). Depending on the value of plot, the method can (and by default does) have a side effect, which is to print gg on the active graphics device.
#'
#' @export
#'
#' @details See \code{\link[vsn]{meanSdPlot}} for further details
meanSdPlot <- function (x, ranks = TRUE, xlab = ifelse(ranks, "rank(mean)",
                                                       "mean"), ylab = "sd", pch, plot = TRUE, bins = 50, ...)
{
  vsn::meanSdPlot(SummarizedExperiment::assay(x), ranks = ranks, xlab = xlab, ylab = ylab,
                  pch = pch, plot = plot, ...)
}

