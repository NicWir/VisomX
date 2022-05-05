####____theme_DEP1____(internal)####
# Adjusted aesthetics for theme_DEP1; adjustable base size
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

####____prot.plot_numbers____####
prot.plot_numbers <- function (se, plot = TRUE, export = FALSE)
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.logical(plot), length(plot) == 1)
  df <- assay(se) %>% data.frame() %>% rownames_to_column() %>%
    gather(ID, bin, -rowname) %>% mutate(bin = ifelse(is.na(bin),
                                                      0, 1))
  stat <- df %>% group_by(ID) %>% summarize(n = n(), sum = sum(bin)) %>%
    left_join(., data.frame(colData(se)), by = "ID")
  p <- ggplot(stat, aes(x = ID, y = sum, fill = condition)) +
    geom_col() + geom_hline(yintercept = unique(stat$n)) +
    labs(title = "Proteins per sample", x = "",
         y = "Number of proteins") + theme_DEP2()

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
####____prot.plot_imputation____####
prot.plot_imputation <- function(se, ..., plot = TRUE, basesize = 12, export = TRUE)
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
               "replicate") %in% colnames(colData(x)))) {
      stop("'label', 'ID', 'condition' and/or 'replicate' ",
           "columns are not present in (one of) the input object(s)",
           "\nRun make_se() or make_se_parse() to obtain the required columns",
           call. = FALSE)
    }
  })
  gather_join <- function(se) {
    assay(se) %>% data.frame(check.names = FALSE) %>% gather(ID, val) %>% left_join(.,
                                                                                    data.frame(colData(se)), by = "ID", check.names = FALSE)
  }
  df <- map_df(arglist, gather_join, .id = "var") %>%
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
  } else {
    return(p)
  }
}

####____prot.plot_density____####
prot.plot_density <- function(se1,
                              se2 = NULL,
                              se3 = NULL,
                              title1 = deparse(substitute(se1)),
                              title2 = deparse(substitute(se2)),
                              title3 = deparse(substitute(se3)),
                              group = FALSE,
                              basesize = 15,
                              plot = plot,
                              export = export) {
  ls <- list()
  if (group == F) {


    # Create data frame from SummarizedExperiment with normalized data
    se.df1 <- data.frame(assay(se1), check.names = FALSE) %>% rownames_to_column()
    if (!is.null(se2)){
      se.df2 <- data.frame(assay(se2), check.names = FALSE) %>% rownames_to_column()
    }
    if (!is.null(se3)){
      se.df3 <- data.frame(assay(se3), check.names = FALSE) %>% rownames_to_column()
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
      pal <- brewer.pal(n = 8, name = "Dark2")
    } else if (length(unique(se1$condition)) <= 12) {
      pal <- brewer.pal(n = 12, name = "Set3")
    } else {
      pal <- c(
        "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
        "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
        "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
        "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
        "green1", "yellow4", "yellow3", "darkorange4", "brown")
    }
    pal <- brewer.pal(n = 8, name = "Dark2")
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
      suppressWarnings(print(marrangeGrob(ls[ls.dens], layout_matrix = lay, top = NULL)))
      grDevices::dev.off()
      grDevices::png("Plots/DensityPlot.png",
          width = width, height = height, units = 'in', res = 300)
      suppressWarnings(print(marrangeGrob(ls[ls.dens], layout_matrix = lay, top = NULL)))
      grDevices::dev.off()
    }

    if (plot == TRUE){
      suppressWarnings(print(marrangeGrob(ls[ls.dens], layout_matrix = lay, top = NULL)))
    }

  } else {
    if (is.null(se2) && is.null(se3)){
      p <- plot_imputation(se1, export = F, plot = F)
      width = 12
      height = 6
    } else if (!is.null(se2) && is.null(se3)){
      p <- plot_imputation(se1, se2, export = F, plot = F)
      width = 12
      height = 12
    } else if (!is.null(se2) && !is.null(se3)){
      p <- plot_imputation(se1, se2, se3, export = F, plot = F)
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
    } else {
      return(p)
    }
  }
}

####____prot.plot_missval____####
# Modified plot_missval function from package DEP with export as PDF and PNG
prot.plot_missval <- function (se, plot = TRUE, export = TRUE, fontsize = 12)
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  se_assay <- assay(se)
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
####____prot.plot_detect____####
# Modified plot_detect function from package DEP with export as PDF and PNG. Plots the group-wise density of proteins with missing values
prot.plot_detect <- function (se, basesize = 10, plot = TRUE, export = TRUE)
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  se_assay <- assay(se)
  if (!any(is.na(se_assay))) {
    stop("No missing values in '", deparse(substitute(se)),
         "'", call. = FALSE)
  }
  df <- se_assay %>% data.frame(check.names = FALSE) %>% rownames_to_column() %>%
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

####____prot.plot_pca____####
# Edited PCA plot function from package DEP including vertical and horizontal reference lines at x = 0 and y = 0
prot.plot_pca <- function (dep,
                           x = 1,
                           y = 2,
                           indicate = c("condition", "replicate"),
                           title = NULL,
                           label = FALSE,
                           n = ncol(t(assay(dep))),
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
        "Run plot_pca() with 'x' and 'y' <= ",
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
        "Run plot_pca() with 'n' <= ",
        nrow(dep),
        "."
      ),
      call. = FALSE
    )
  }
  columns <- colnames(colData(dep))
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
  var <- apply(assay(dep), 1, sd)
  df <- assay(dep)[order(var, decreasing = TRUE)[seq_len(n)],]
  pca <- stats::prcomp(t(df), scale = FALSE)
  pca_df <- pca$x %>% data.frame(check.names = FALSE) %>% rownames_to_column() %>%
    left_join(., data.frame(colData(dep), check.names = FALSE), by = c(rowname = "ID"))
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
      scale_shape_manual(values=c(21, 22, 23, 24, 25, 0, 1, 2, 5, 6))

    if (length(unique(dep$condition)) <= 8) {
      p <- p + scale_fill_brewer(palette = "Dark2", guide = guide_legend(override.aes = list(shape = 21)))+
        scale_color_manual(values=c(rep("black", length(unique(dep$condition)))))
    } else if (length(unique(dep$condition)) <= 12) {
      p <- p + scale_fill_brewer(palette = "Set3", guide = guide_legend(override.aes = list(shape = 21))) +
        scale_color_manual(values=c(rep("black", length(unique(dep$condition)))))
    } else {
      pal <- c(
        "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
        "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
        "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
        "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
        "green1", "yellow4", "yellow3", "darkorange4", "brown")
      p <- p + scale_fill_manual(values = pal[1:length(unique(dep$condition))],
                                 guide = guide_legend(override.aes = list(shape = 21))) +
        scale_color_manual(values=c(rep("black", length(unique(dep$condition)))))
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

####____prot.plot_volcano____####
#Revised volcano plot function from the DEP package with enhanced aesthetics.

prot.plot_volcano <-
  function (dep,
            contrast,
            label_size = 3,
            alpha = 0.05,
            lfc = 2,
            add_names = TRUE,
            adjusted = FALSE,
            plot = TRUE,
            export = TRUE)
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
    row_data <- rowData(dep, use.names = FALSE)
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
      filter(!is.na(significant)) %>% arrange(significant)

    # create new columns for shapes "circle" and "triangle" if proteins lie within or outside (for proteins with very small adj.p value < 3.162278e-05) of the y bounds, respectively.
    if (adjusted) {
      df$shape <- ifelse(df$y > 4.9, "triangle", "circle")
    }
    else{
      df$shape <- "circle"
    }
    # change the -Log10(q value) of proteins exceeding the y plot bound at 5.0 to 4.9 so that they are displayed at the axis border
    #df$y[df$y>4.9] <- 4.9
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
        labels = c("Upregulated", "Not significant", "Downregulated"),
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
      force_panelsizes(rows = unit(6.5, "in"),
                       cols = unit(6.5, "in")) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = max(abs(df$x)) + 1))

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
        data = filter(df_volcano, significant),
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
      colnames(df)[c(1, 2, 3)] <- c("protein", "log2_fold_change",
                                    "p_value_-log10")
      if (adjusted) {
        colnames(df)[3] <- "adjusted_p_value_-log10"
      }
      return(df)
    }
  }

####____prot.plot_corrheatmap____####
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
    stop("'lower' and/or 'upper' arguments are not valid\n         Run plot_pca() with 'lower' and 'upper' between -1 and 1",
         call. = FALSE)
  }
  pals <- RColorBrewer::brewer.pal.info %>% rownames_to_column() %>%
    filter(category != "qual")
  if (!pal %in% pals$rowname) {
    stop("'", pal, "' is not a valid color panel",
         " (qualitative panels also not allowed)\n",
         "Run plot_pca() with one of the following 'pal' options: ",
         paste(pals$rowname, collapse = "', '"), "'",
         call. = FALSE)
  }
  if (any(is.na(assay(dep)))) {
    stop("Missing values in '", deparse(substitute(dep)),
         "'. Use plot_dist() instead")
  }
  if (!is.null(indicate)) {
    assertthat::assert_that(is.character(indicate))
    col_data <- colData(dep) %>% as.data.frame()
    columns <- colnames(col_data)
    if (any(!indicate %in% columns)) {
      stop("'", paste0(indicate, collapse = "' and/or '"),
           "' column(s) is/are not present in ", deparse(substitute(dep)),
           ".\nValid columns are: '", paste(columns,
                                            collapse = "', '"), "'.", call. = FALSE)
    }
    anno <- colData(dep) %>% data.frame() %>% select(indicate)
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
    if (!"significant" %in% colnames(rowData(dep, use.names = FALSE))) {
      stop("'significant' column is not present in '",
           deparse(substitute(dep)), "'\nRun add_rejections() to obtain the required column",
           call. = FALSE)
    }
    dep <- dep[rowData(dep, use.names = FALSE)$significant,
    ]
  }
  cor_mat <- cor(assay(dep))
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
####____prot.plot_heatmap____####
# Modified plot_heatmap function with adjusted aesthetics and legend breaks
prot.plot_heatmap <- function (dep, type = c("contrast", "centered"),
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
                               export = TRUE, ...)
{
  get_annotation_contrast <- function (dep, indicate, contrast = contrast_samples)
  {
    assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                            is.character(indicate))
    col_data <- colData(dep) %>% data.frame(check.names = FALSE)
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
    anno <- select(col_data, indicate)
    anno <- filter(anno, str_detect(condition, paste(contrast, collapse = "|")))
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
    contrast <- rowData(dep, use.names = FALSE) %>% data.frame(check.names = FALSE) %>% select(ends_with("_diff")) %>%
      colnames(.) %>% gsub("_diff", "", .)
  }
  if (length(grep(paste(paste(contrast, collapse = "|"), "_diff", sep = ""),
                  colnames(rowData(dep, use.names = FALSE)))) == 0) {
    valid_cntrsts <-
      rowData(dep, use.names = FALSE) %>% data.frame(check.names = FALSE) %>% select(ends_with("_diff")) %>%
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
  row_data <- rowData(dep, use.names = FALSE)
  col_data <- colData(dep) %>% data.frame(check.names = FALSE)
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
  if (nrow(assay(filtered)) == 0){
    stop("No proteins with significantly different abundance were found.")
  }
  if (any(is.na(assay(filtered)))) {
    warning("Missing values in '", deparse(substitute(dep)),
            "'. ", "Using clustering_distance = 'gower'",
            call. = FALSE)
    clustering_distance <- "gower"
    obs_NA <- TRUE
  } else {
    obs_NA <- FALSE
  }
  if (type == "centered") {
    rowData(filtered)$mean <- rowMeans(assay(filtered), na.rm = TRUE)
    df <- assay(filtered) - rowData(filtered, use.names = FALSE)$mean
    if (!is.null(contrast) && show_all == FALSE){
      df <- data.frame(df, check.names = FALSE)[,str_detect(colnames(df),paste(contrast_samples, collapse = "|"))]
    }
  }
  if (type == "contrast") {
    df <- rowData(filtered, use.names = FALSE) %>% data.frame(check.names = FALSE) %>%
      column_to_rownames(var = "name") %>% select(ends_with("_diff"))
    if(length(grep("_diff", colnames(rowData(filtered)))) == 1){
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
    df_cent <- assay(filtered) - rowMeans(assay(filtered))
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
      len = nrow(dep[rowData(dep)$significant, ]) / 12
    } else {
      len = nrow(dep[rowData(dep)$significant, ]) / 45
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
    data.frame(protein = row.names(return), return, check.names = FALSE) %>% mutate(order = row_number())
  }
}

####____prot.plot_heatmap_all____####
# Modified plot_heatmap function with adjusted aesthetics and legend breaks
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
                                   export = TRUE,
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
  row_data <- rowData(se, use.names = FALSE)
  col_data <- colData(se) %>% data.frame(check.names = FALSE)
  if (any(!c("label", "condition", "replicate") %in%
          colnames(col_data))) {
    stop(paste0("'label', 'condition' and/or 'replicate' columns are not present in '",
                deparse(substitute(se)), "'"), call. = FALSE)
  }
  ha1 <- get_annotation(se, indicate)
  obs_NA <- FALSE
  rowData(se)$mean <- rowMeans(assay(se), na.rm = TRUE)
  df <- assay(se) - rowData(se, use.names = FALSE)$mean


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
    df_cent <- assay(se) - rowMeans(assay(se))
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
      len = nrow(se[rowData(se)$significant, ]) / 12
    } else {
      len = nrow(se[rowData(se)$significant, ]) / 45
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
    data.frame(protein = row.names(return), return, check.names = FALSE) %>% mutate(order = row_number())
  }
}


####____prot.plot_screeplot____####
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
                                 borderColour = "black", plot = TRUE, export = TRUE)
{
  PC <- Variance <- NULL
  components <- getComponents(pcaobj)[1:length(pcaobj$components)]
  if(length(getComponents(pcaobj)) > 10){
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
####____prot.plot_loadings____####
prot.plot_loadings <- function (pcaobj, components = getComponents(pcaobj, seq_len(5)),
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
                                borderWidth = 0.8, borderColour = "black", plot = TRUE, export = TRUE)
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
    p <- p + geom_text_repel(data = plotobj, aes(label = as.character(var)),
                             size = labSize, nudge_x = ifelse(positionConnectors ==
                                                                "left", -0.75, ifelse(positionConnectors ==
                                                                                        "right", 0.75, 0)), direction = "y",
                             segment.color = colConnectors, segment.size = widthConnectors,
                             arrow = arrow(length = lengthConnectors, type = typeConnectors,
                                           ends = endsConnectors), show.legend = FALSE,
                             hjust = labhjust, vjust = labvjust)
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
####____prot.plot_enrichment____####
prot.plot_enrichment <- function(enrichset,
                                 title = "Differentially enriched pathways",
                                 subtitle = "",
                                 plot = TRUE,
                                 export = FALSE,
                                 kegg = TRUE) {

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

####____prot.plot_upset____####
prot.plot_upset <- function(enrichset, order.by = "freq", point.size = 3,
                            line.size = 1, text.scale = c(2, 2, 2, 2, 2, 1.5), ...)
{
  upsetlist <- enrichset$geneID %>% str_replace_all(., "/", ", ") %>% strsplit(., ", ")
  names(upsetlist) <- enrichset$Description
  UpSetR::upset(fromList(upsetlist), order.by = order.by, nsets = length(upsetlist),
                point.size = point.size, line.size = line.size, text.scale = text.scale, ...)
}

####____prot.plot_bar____####
# Plot bar plots of single proteins based on their Ensembl gene ID
prot.plot_bar <- function (dep,
                           proteins,
                           type = c("contrast", "centered"),
                           contrast = NULL,
                           col.id = "ID",
                           match.id = "Accession",
                           match.name = "Name",
                           convert_name = FALSE,
                           shape.size = 2.5,
                           name_table = NULL,
                           plot = TRUE,
                           export = TRUE,
                           export.nm = NULL)
{
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(proteins), is.character(type), is.logical(plot),
                          length(plot) == 1)
  require("ggnewscale")
  type <- match.arg(type)
  # Replace rownames of dep with column of original data frame defined in "col.id" argument
  row_data <- rowData(dep, use.names = F)
  GeneID <- match(col.id,
                  colnames(row_data))
  rownames(dep) <- row_data[, GeneID]
  rownames(dep)[is.na(rownames(dep))] <- 0
  row_data <- rowData(dep, use.names = FALSE)

  if (any(!c("label", "condition", "replicate") %in%
          colnames(colData(dep)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(dep)), "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1 && type == "contrast") {
    stop("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(dep)), "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }
  if (!"name" %in% colnames(row_data)) {
    stop("'name' column not present in '", deparse(substitute(dep)),
         "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  if (all(!proteins %in% row_data[,col.id])) {
    if (length(proteins) == 1) {
      rows <- grep(substr(proteins, 1, nchar(proteins) -
                            1), row_data[,col.id])
      possibilities <- row_data[,col.id][rows]
    }
    else {
      rows <- lapply(proteins, function(x) grep(substr(x,
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
    stop(paste0("The proteins ", paste(proteins, collapse=", "), " were not found in the column '", col.id,
                "'. Please run `prot.plot_bar()` with a valid protein names in the 'proteins' argument or provide the valid column ID in the 'col.id' argument.\n"),
         possibilities_msg, call. = FALSE)
  }
  if (any(!proteins %in% row_data[,col.id])) {
    proteins <- proteins[proteins %in% row_data[,col.id]]
    warning("Only used the following protein(s): '",
            paste0(proteins, collapse = "', '"), "'")
  }
  subset_list <- list()
  if(length(proteins)<=8){
    subset_list[[1]] <- dep[proteins]
  } else if(length(proteins)<=16){
    subset_list[[1]] <- dep[proteins[1:8]]
    subset_list[[2]] <- dep[proteins[9:length(proteins)]]
  } else if(length(proteins)<=24){
    subset_list[[1]] <- dep[proteins[1:8]]
    subset_list[[2]] <- dep[proteins[9:16]]
    subset_list[[3]] <- dep[proteins[17:length(proteins)]]
  }else if(length(proteins)<=32){
    subset_list[[1]] <- dep[proteins[1:8]]
    subset_list[[2]] <- dep[proteins[9:16]]
    subset_list[[3]] <- dep[proteins[17:24]]
    subset_list[[4]] <- dep[proteins[25:length(proteins)]]
  } else if(length(proteins)<=40){
    subset_list[[1]] <- dep[proteins[1:8]]
    subset_list[[2]] <- dep[proteins[9:16]]
    subset_list[[3]] <- dep[proteins[17:24]]
    subset_list[[4]] <- dep[proteins[25:32]]
    subset_list[[5]] <- dep[proteins[33:length(proteins)]]
  } else if(length(proteins)<=48){
    subset_list[[1]] <- dep[proteins[1:8]]
    subset_list[[2]] <- dep[proteins[9:16]]
    subset_list[[3]] <- dep[proteins[17:24]]
    subset_list[[4]] <- dep[proteins[25:32]]
    subset_list[[5]] <- dep[proteins[33:40]]
    subset_list[[6]] <- dep[proteins[41:length(proteins)]]
  } else if(length(proteins)<=56){
    subset_list[[1]] <- dep[proteins[1:8]]
    subset_list[[2]] <- dep[proteins[9:16]]
    subset_list[[3]] <- dep[proteins[17:24]]
    subset_list[[4]] <- dep[proteins[25:32]]
    subset_list[[5]] <- dep[proteins[33:40]]
    subset_list[[6]] <- dep[proteins[41:48]]
    subset_list[[7]] <- dep[proteins[49:length(proteins)]]
  } else {
    subset_list[[1]] <- dep[proteins[1:8]]
    subset_list[[2]] <- dep[proteins[9:16]]
    subset_list[[3]] <- dep[proteins[17:24]]
    subset_list[[4]] <- dep[proteins[25:32]]
    subset_list[[5]] <- dep[proteins[33:40]]
    subset_list[[6]] <- dep[proteins[41:48]]
    subset_list[[7]] <- dep[proteins[49:56]]
    subset_list[[8]] <- dep[proteins[57:length(proteins)]]
  }

  for(i in 1:length(subset_list)) {
    if (type == "centered") {
      means <- rowMeans(assay(subset_list[[i]]), na.rm = TRUE)
      df_reps <-
        data.frame(assay(subset_list[[i]]) - means, check.names = FALSE) %>% rownames_to_column() %>%
        gather(ID, val, -rowname) %>% left_join(., data.frame(colData(subset_list[[i]]), check.names = FALSE),
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
        df_reps %>% group_by(condition, rowname) %>% summarize(
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
      if (length(unique(dep$condition)) <= 8) {
        p <- p + scale_fill_manual(values = brewer.pal(n = length(unique(dep$condition)), name = "Dark2"))
      } else if (length(unique(dep$condition)) <= 12) {
        p <- p + scale_fill_manual(values = brewer.pal(n = length(unique(dep$condition)), name = "Set3"))
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
        new_scale_fill() +
        geom_point(
          data = df_reps,
          aes(condition, val, fill = replicate),
          shape = 23,
          size = shape.size,
          color = "black",
          position = position_dodge(width = 0.3)
        ) +
        scale_fill_manual(values = case_when(
          as.numeric(max(dep$replicate))<=8       ~ brewer.pal(n=as.numeric(max(dep$replicate)), name="Greys"),
          as.numeric(max(dep$replicate))>8        ~ colorRampPalette(brewer.pal(n=8, name="Greys"))(as.numeric(max(dep$replicate))))
        ) +
        geom_errorbar(aes(ymin = CI.L, ymax = CI.R), width = 0.3) +
        labs(
          x = "Baits",
          y = expression(log[2] ~ "Centered intensity" ~
                           "(?95% CI)"),
          col = "Rep"
        ) + facet_wrap( ~ rowname) +
        theme(basesize = 12) + theme_DEP2()
      w <- 8+log(max(str_count(df[,"condition"]))-3, base = 1.6)
      h <- 10
    }
    if (type == "contrast") {
      if (convert_name == TRUE) {
        if (!is.null(contrast)) {
          df <-
            rowData(subset_list[[i]], use.names = FALSE) %>% data.frame(check.names = FALSE) %>%
            select(
              name,
              paste0(contrast, "_diff"),
              paste0(contrast, "_CI.L"),
              paste0(contrast, "_CI.R")
            ) %>% gather(var, val, -name) %>%
            mutate(
              contrast = gsub("_diff|_CI.L|_CI.R", "", var),
              var = gsub(".*_", "", var)
            ) %>%
            spread(var, val)
        } else {
          df <-
            rowData(subset_list[[i]], use.names = FALSE) %>% data.frame(check.names = FALSE) %>%
            select(name,
                   ends_with("_diff"),
                   ends_with("_CI.L"),
                   ends_with("_CI.R")) %>% gather(var, val, -name) %>%
            mutate(
              contrast = gsub("_diff|_CI.L|_CI.R", "", var),
              var = gsub(".*_", "", var)
            ) %>%
            spread(var, val)
        }
        p <-
          ggplot(df, aes(contrast, diff)) + geom_hline(yintercept = 0) +
          geom_col(aes(y = diff, fill = contrast), colour = "black")
        if (length(unique(dep$condition)) <= 8) {
          p <- p + scale_fill_manual(values = brewer.pal(n = length(unique(dep$condition)), name = "Dark2"))
        } else if (length(unique(dep$condition)) <= 12) {
          p <- p + scale_fill_manual(values = brewer.pal(n = length(unique(dep$condition)), name = "Set3"))
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
          geom_errorbar(aes(ymin = CI.L, ymax = CI.R), width = 0.3) +
          labs(
            x = "Contrasts",
            y = expression(log[2] ~ "Fold change" ~
                             "(?95% CI)")) + facet_wrap( ~ name) +
          theme(basesize = 12) + theme_DEP2()


      } else {
        if (!is.null(contrast)) {
          df <-
            rowData(subset_list[[i]], use.names = FALSE) %>% data.frame(check.names = FALSE) %>%
            select(
              ID,
              paste0(contrast, "_diff"),
              paste0(contrast, "_CI.L"),
              paste0(contrast, "_CI.R")
            ) %>% gather(var, val, -ID) %>%
            mutate(
              contrast = gsub("_diff|_CI.L|_CI.R", "", var),
              var = gsub(".*_", "", var)
            ) %>%
            spread(var, val)
        } else {
          df <-
            rowData(subset_list[[i]], use.names = FALSE) %>% data.frame(check.names = FALSE) %>%
            select(ID,
                   ends_with("_diff"),
                   ends_with("_CI.L"),
                   ends_with("_CI.R")) %>% gather(var, val, -ID) %>%
            mutate(
              contrast = gsub("_diff|_CI.L|_CI.R", "", var),
              var = gsub(".*_", "", var)
            ) %>%
            spread(var, val)
        }

        p <-
          ggplot(df, aes(contrast, diff)) + geom_hline(yintercept = 0) +
          geom_col(aes(y = diff, fill = contrast), colour = "black")
        if (length(unique(dep$condition)) <= 8) {
          p <- p + scale_fill_manual(values = brewer.pal(n = length(unique(dep$condition)), name = "Dark2"))
        } else if (length(unique(dep$condition)) <= 12) {
          p <- p + scale_fill_manual(values = brewer.pal(n = length(unique(dep$condition)), name = "Set3"))
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
          geom_errorbar(aes(ymin = CI.L, ymax = CI.R), width = 0.3) +
          labs(
            x = "Contrasts",
            y = expression(log[2] ~ "Fold change" ~
                             "(?95% CI)")) + facet_wrap( ~ ID) +
          theme(basesize = 12) + theme_DEP2()


      }
      w <- 8+log(max(str_count(df[,"contrast"]))-3, base = 1.6)
      h <- 10
    }

    if (export == TRUE) {
      if (!is.null(export.nm)) {
        if (length(proteins) <= 8) {
          export_name <- paste0("Plots/", export.nm)
        } else{
          export_name <- paste0("Plots/", export.nm, "_", i)
        }
      } else {
        if (length(proteins) <= 8) {
          export_name <- paste0("Plots/BarPlot_",
                                paste(proteins, collapse = "_"))
        } else{
          export_name <- paste0("Plots/BarPlot_",
                                paste(proteins, collapse = "_"), "_", i)
        }
      }
      dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
      message(
        paste0(
          "Exporting bar plot(s) to:\n\"",
          getwd(),
          export_name,
          ".pdf",
          "and \"...",
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
    if(length(proteins)>8){
      plot = FALSE
    }
    if (plot) {
      return(p)
    } else {
      if (type == "centered") {
        df <- df %>% select(rowname, condition, mean, CI.L,
                            CI.R)
        colnames(df) <- c("protein", "condition",
                          "log2_intensity", "CI.L", "CI.R")
      }
      if (type == "contrast") {
        df <- df %>% select(name, contrast, diff, CI.L, CI.R) %>%
          mutate(contrast = paste0(contrast))
        colnames(df) <- c("protein", "contrast",
                          "log2_fold_change", "CI.L", "CI.R")
      }
      if (length(proteins) <= 8) {
        return(df)
      }
    }
  }
}


####____prot.boxplot_intensity____####
# Modified plot_normalization function from package DEP with adjusted colors and
# function to export plot as PDF and PNG
prot.boxplot_intensity <- function (se, ..., plot = TRUE, export = TRUE) {
  call <- match.call()
  call$export <- NULL
  call$plot <- NULL
  arglist <- lapply(call[-1], function(x) x)
  var.names <- vapply(arglist, deparse, character(1))
  arglist <- lapply(arglist, eval.parent, n = 2)
  names(arglist) <- var.names
  lapply(arglist, function(x) {
    assertthat::assert_that(inherits(x, "SummarizedExperiment"),
                            msg = "input objects need to be of class 'SummarizedExperiment'")
    if (any(!c("label", "ID", "condition",
               "replicate") %in% colnames(colData(x)))) {
      stop("'label', 'ID', 'condition' and/or 'replicate' ",
           "columns are not present in (one of) the input object(s)",
           "\nRun make_se() or make_se_parse() to obtain the required columns",
           call. = FALSE)
    }
  })
  gather_join <- function(se) {
    assay(se) %>% data.frame(check.names = FALSE) %>% gather(ID, val) %>%
      left_join(., data.frame(colData(se), check.names = FALSE), by = "ID")
  }
  df <- map_df(arglist, gather_join, .id = "var") %>%
    mutate(var = factor(var, levels = names(arglist)))
  p <- ggplot(df, aes(x = ID, y = val, fill = condition)) +
    geom_boxplot(notch = TRUE, na.rm = TRUE) +
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
}
