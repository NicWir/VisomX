#' Run upon attaching package VisomX
#'
#' Changes debug option for package \code{rgl} to avoid Rstudio crashing upon attaching it and prints welcome message
#'
#' @param libname library name
#' @param pkgname package name
.onAttach <- function (libname, pkgname){
  options(rgl.debug = TRUE)
  k1 <- paste("VisomX",utils::packageVersion( "VisomX"),"initialized Successfully !")
  k0 <- "\n"
  k2 <- paste("https://github.com/NicWir/VisomX")
  packageStartupMessage(c(k1,k0,k2))
}

#' Format font color for Markdown reports
#'
#' \code{colFmt} formats the input depending on PDF or HTML output to give colored text in reports.
#'
#' @param x A character string. The text to be colored.
#' @param color (Character) A color.
colFmt <- function(x, color) {
  outputFormat <- knitr::opts_knit$get("rmarkdown.pandoc.to")

  if (outputFormat == "latex") {
    ret <- paste("\\textcolor{", color, "}{", gsub("_", "\\\\_", x), "}", sep = "")
  } else if (outputFormat == "html") {
    ret <- paste("<font color='", color, "'>", x, "</font>", sep = "")
  } else {
    ret <- x
  }
  return(ret)
}

#' Write object in CSV file
#'
#' @param dat An R data object (e.g., list, data frame)
#' @param file (Character) The name of the CSV file.
#' @param row.names (Logical) Add row names as column (\code{TRUE}) or not (\code{FALSE}).
fast.write.csv <- function(dat, file, row.names = TRUE) {
  tryCatch(
    {
      if (is.data.frame(dat)) {
        # there is a rare bug in data.table (R 3.6) which kill the R process in some cases
        data.table::fwrite(dat, file, row.names = row.names)
      } else {
        utils::write.csv(dat, file, row.names = row.names)
      }
    },
    error = function(e) {
      print(e)
      write.csv(dat, file, row.names = row.names)
    },
    warning = function(w) {
      print(w)
      write.csv(dat, file, row.names = row.names)
    }
  )
}

GetShapeSchema <- function(mSetObj=NA, show.name, grey.scale){
  if(exists("shapeVec") && all(shapeVec >= 0)){
    sps <- rep(0, length=length(mSetObj$dataSet$cls));
    clsVec <- as.character(mSetObj$dataSet$cls)
    grpnms <- names(shapeVec);
    for(i in 1:length(grpnms)){
      sps[clsVec == grpnms[i]] <- shapeVec[i];
    }
    shapes <- sps;
  }else{
    if(show.name | grey.scale){
      shapes <- as.numeric(mSetObj$dataSet$cls)+1;
    }else{
      shapes <- rep(21, length(mSetObj$dataSet$cls));
    }
  }
  return(shapes);
}

#' Multiple set version of intersect
#'
#' @param x A list
#' @export
Intersect <- function (x) {
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

#' Remove the union of the y's from the common x's.
#'
#' @param x A list of characters
#' @param y A list of characters
#' @export
Setdiff <- function (x, y) {
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}

#' Multiple set version of union
#'
#' @param x A list
#' @export
Union <- function (x) {
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

move.file <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir)
  file.rename(from = from,  to = to)
}

#' Multiple set version of union
#'
#' @param x A list
#' @export
prot.get_kegg_pathways <- function(organism){
  assertthat::assert_that(is.character(organism))
  # Get gene sets from KEGG
  kegg <- download_KEGG(organism)
  # replace KEGG pathway numbers with names
  kegg$KEGGPATHID2EXTID$from <- kegg$KEGGPATHID2NAME[match(kegg$KEGGPATHID2EXTID$from,
                                                           kegg$KEGGPATHID2NAME$from), 2]
  # Create list of KEGG gene sets
  kegg.gs <- split(kegg$KEGGPATHID2EXTID[, 2], kegg$KEGGPATHID2EXTID[, 1])

  return(kegg.gs)
}

####____prot.get_pathway_genes____####
prot.get_pathway_genes <- function(pathway_name, pathway_table, colid_pathways, colid_genes, gene_sep = ", "){
  genes <- unlist(str_split(pathway_table[match(pathway_name, pathway_table[[colid_pathways]]),
                                          colid_genes], gene_sep))

  return(genes)
}

GetDefaultPLSCVComp <- function(mSetObj = NA){
  return(min(5, dim(mSetObj$dataSet$norm)[1] - 2, dim(mSetObj$dataSet$norm)[2],
             mSetObj$dataSet$min.grp.size))
}

.get.ttest.res <- function (data, inx1, inx2, paired = FALSE, equal.var = TRUE,
                            nonpar = F){
  print("Performing regular t-tests ....")
  univ.test <- function(x) {
    t.test(x[inx1], x[inx2], paired = paired, var.equal = equal.var)
  }
  if (nonpar) {
    univ.test <- function(x) {
      wilcox.test(x[inx1], x[inx2], paired = paired)
    }
  }
  my.fun <- function(x) {
    tmp <- try(univ.test(x))
    if (class(tmp) == "try-error") {
      return(c(NA, NA))
    }
    else {
      return(c(tmp$statistic, tmp$p.value))
    }
  }
  res <- apply(data, 2, my.fun)
  return(t(res))
}

ReplaceMissingByLoD <- function (int.mat){
  int.mat <- as.matrix(int.mat)
  rowNms <- rownames(int.mat)
  colNms <- colnames(int.mat)
  int.mat <- apply(int.mat, 2, .replace.by.lod)
  rownames(int.mat) <- rowNms
  colnames(int.mat) <- colNms
  return(int.mat)
}

CleanDataMatrix <- function (ndata){
  varCol <- apply(data.frame(ndata), 2, var, na.rm = T)
  constCol <- (varCol == 0 | is.na(varCol))
  return(ndata[, !constCol, drop = FALSE])
}

CleanData <- function (bdata, removeNA = T, removeNeg = T, removeConst = T){
  if (sum(bdata == Inf, na.rm = TRUE) > 0) {
    inx <- bdata == Inf
    bdata[inx] <- NA
    bdata[inx] <- max(bdata, na.rm = T) * 2
  }
  if (sum(bdata == -Inf, na.rm = TRUE) > 0) {
    inx <- bdata == -Inf
    bdata[inx] <- NA
    bdata[inx] <- min(bdata, na.rm = T)/2
  }
  if (removeNA) {
    if (sum(is.na(bdata)) > 0) {
      bdata[is.na(bdata)] <- min(bdata, na.rm = T)/2
    }
  }
  if (removeNeg) {
    if (sum(as.numeric(bdata <= 0)) > 0) {
      inx <- bdata <= 0
      bdata[inx] <- NA
      bdata[inx] <- min(bdata, na.rm = T)/2
    }
  }
  if (removeConst) {
    varCol <- apply(data.frame(bdata), 2, var, na.rm = T)
    constCol <- (varCol == 0 | is.na(varCol))
    constNum <- sum(constCol, na.rm = T)
    if (constNum > 0) {
      bdata <- data.frame(bdata[, !constCol, drop = FALSE],
                          check.names = F)
    }
  }
  bdata
}

CalculatePairwiseDiff <- function (mat){
  f <- function(i, mat) {
    z <- mat[, i - 1] - mat[, i:ncol(mat), drop = FALSE]
    colnames(z) <- paste(colnames(mat)[i - 1], colnames(z),
                         sep = "/")
    z
  }
  res <- do.call("cbind", sapply(2:ncol(mat), f, mat))
  round(res, 5)
}

Get.Fstat <- function (x, fac, var.equal = TRUE) {
  x = t(x)
  sqr = function(x) x * x
  stopifnot(length(fac) == ncol(x), is.factor(fac), is.matrix(x))
  x <- x[, !is.na(fac), drop = FALSE]
  fac <- fac[!is.na(fac)]
  k <- nlevels(fac)
  xm <- matrix(sapply(levels(fac), function(fl) rowMeans(x[,
                                                           which(fac == fl), drop = FALSE])), nrow = nrow(x), ncol = nlevels(fac))
  x1 <- xm[, fac, drop = FALSE]
  dff <- k - 1
  x0 <- matrix(rowMeans(x), ncol = ncol(x), nrow = nrow(x))
  dfr <- ncol(x) - dff - 1
  mssf <- rowSums(sqr(x1 - x0))/dff
  mssr <- rowSums(sqr(x - x1))/dfr
  fstat <- mssf/mssr
  return(fstat)
}

GetColorSchema <- function (my.cls, grayscale = F){
  lvs <- levels(my.cls)
  grp.num <- length(lvs)
  if (grayscale) {
    dist.cols <- colorRampPalette(c("grey90", "grey30"))(grp.num)
  }
  else if (exists("colVec") && !any(colVec == "#NA") && length(colVec) ==
           length(levels(my.cls))) {
    dist.cols <- colVec
  }
  else {
    if (grp.num <= 18) {
      dist.cols <- pal_18[1:grp.num]
    }
    else {
      dist.cols <- colorRampPalette(pal_18)(grp.num)
    }
  }
  colors <- vector(mode = "character", length = length(my.cls))
  for (i in 1:length(lvs)) {
    colors[my.cls == lvs[i]] <- dist.cols[i]
  }
  return(colors)
}

get_pheatmap_dims <- function (dat, annotation, view.type, width, cellheight = 15,
                               cellwidth = 15) {
  png("NUL")
  heat_map <- pheatmap::pheatmap(dat, annotation = annotation,
                                 cellheight = cellheight, cellwidth = cellwidth)
  h <- sum(sapply(heat_map$gtable$heights, grid::convertHeight,
                  "in"))
  w <- sum(sapply(heat_map$gtable$widths, grid::convertWidth,
                  "in"))
  dev.off()
  myW <- ncol(dat) * 20 + 200
  if (myW < 650) {
    myW <- 650
  }
  myW <- round(myW/72, 2)
  if (w < myW) {
    w <- myW
  }
  if (view.type == "overview") {
    if (is.na(width)) {
      if (w > 9) {
        w <- 9
      }
    }
    else if (width == 0) {
      if (w > 7.2) {
        w <- 7.2
      }
    }
    else {
      w <- 7.2
    }
    if (h > w) {
      h <- w
    }
  }
  return(list(height = h, width = w))
}

GetRGBColorGradient <- function (vals) {
  library(RColorBrewer)
  seed.cols <- c("#FCF5DF", "#FFEDA0", "#F03B20")
  cols <- colorRampPalette(seed.cols)(length(vals))
  my.alpha <- signif(seq(from = 0.3, to = 0.8, length.out = length(vals)),
                     2)
  rgb.cols <- my.col2rgba(cols, alpha = my.alpha)
  nms.orig <- names(vals)
  names(rgb.cols) <- names(sort(vals))
  ord.cols <- rgb.cols[nms.orig]
  return(as.vector(ord.cols))
}

all.numeric <- function (x, what = c("test", "vector"), extras = c(".", "NA")) {
  what <- match.arg(what)
  old <- options(warn = -1)
  on.exit(options(old))
  x <- sub("[[:space:]]+$", "", x)
  x <- sub("^[[:space:]]+", "", x)
  inx <- x %in% c("", extras)
  xs <- x[!inx]
  isnum <- !any(is.na(as.numeric(xs)))
  if (what == "test")
    isnum
  else if (isnum)
    as.numeric(x)
  else x
}
