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
