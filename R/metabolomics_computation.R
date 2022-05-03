#' Perform ANOVA analysis
#'
#' ANOVA analysis
#'
#' @param mSetObj Input name of the created mSet object. (see \code{\link[VisomX]{met.read_data}})
#' @param nonpar (Logical) Use a non-parametric test (\code{TRUE}) or not (\code{FALSE}).
#' @param thresh (Numeric) From 0 to 1, indicate the p-value threshold.
#' @param post.hoc (Character) Enter the name of the post-hoc test, "fisher" or "tukey".
#' @param all_results (Logical) If \code{TRUE}, it the ANOVA results for all compounds with no post-hoc tests performed will be written as CSV file "anova_all_results.csv".
#' @param silent (Logical) Suppress message with number of significant features found (\code{TRUE}) or not (\code{FALSE}).
#' @return The input mSet object with ANOVA results added at mSetObj$analSet$aov.
#' @export
met.ANOVA.Anal <- function (mSetObj = NA, nonpar = FALSE, thresh = 0.05, post.hoc = "fisher",
          all_results = FALSE, silent = FALSE)
{
  sig.num <- 0
  if (nonpar) {
    aov.nm <- "Kruskal Wallis Test"
    anova.res <- apply(as.matrix(mSetObj$dataSet$norm), 2,
                       kwtest, cls = mSetObj$dataSet$cls)
    res <- unlist(lapply(anova.res, function(x) {
      c(x$statistic, x$p.value)
    }))
    res <- data.frame(matrix(res, nrow = length(anova.res),
                             byrow = T), stringsAsFactors = FALSE)
    fstat <- res[, 1]
    p.value <- res[, 2]
    names(fstat) <- names(p.value) <- colnames(mSetObj$dataSet$norm)
    fdr.p <- p.adjust(p.value, "fdr")
    inx.imp <- fdr.p <= thresh
    sig.num <- sum(inx.imp)
    if (sig.num > 0) {
      sig.f <- fstat[inx.imp]
      sig.p <- p.value[inx.imp]
      fdr.p <- fdr.p[inx.imp]
      sig.mat <- data.frame(signif(sig.f, 5), signif(sig.p,
                                                     5), signif(-log10(sig.p), 5), signif(fdr.p, 5),
                            "NA")
      rownames(sig.mat) <- names(sig.p)
      colnames(sig.mat) <- c("chi.squared", "p.value",
                             "-log10(p)", "FDR", "Post-Hoc")
      ord.inx <- order(sig.p, decreasing = FALSE)
      sig.mat <- sig.mat[ord.inx, , drop = F]
      fileName <- "kw_posthoc.csv"
      my.mat <- sig.mat[, 1:4]
      colnames(my.mat) <- c("chi_squared", "pval_KW",
                            "-log10(p)", "FDR")
    }
  }
  else {
    aov.nm <- "One-way ANOVA"
    aov.res <- apply(as.matrix(mSetObj$dataSet$norm),
                     2, MetaboAnalystR:::aof, cls = mSetObj$dataSet$cls)
    anova.res <- lapply(aov.res, anova)
    res <- unlist(lapply(anova.res, function(x) {
      c(x["F value"][1, ], x["Pr(>F)"][1,
      ])
    }))
    res <- data.frame(matrix(res, nrow = length(aov.res),
                             byrow = T), stringsAsFactors = FALSE)
    fstat <- res[, 1]
    p.value <- res[, 2]
    names(fstat) <- names(p.value) <- colnames(mSetObj$dataSet$norm)
    fdr.p <- p.adjust(p.value, "fdr")
    if (all_results == TRUE) {
      all.mat <- data.frame(signif(p.value, 5), signif(-log10(p.value),
                                                       5), signif(fdr.p, 5))
      rownames(all.mat) <- names(p.value)
      colnames(all.mat) <- c("p.value", "-log10(p)",
                             "FDR")
      fast.write.csv(all.mat, "anova_all_results.csv")
    }
    inx.imp <- fdr.p <= thresh
    sig.num <- sum(inx.imp)
    if (sig.num > 0) {
      aov.imp <- aov.res[inx.imp]
      sig.f <- fstat[inx.imp]
      sig.p <- p.value[inx.imp]
      fdr.p <- fdr.p[inx.imp]
      cmp.res <- NULL
      post.nm <- NULL
      if (post.hoc == "tukey") {
        tukey.res.imp <- lapply(aov.imp, TukeyHSD, conf.level = 1 -
                                  thresh)
        cmp.res.imp <- unlist(lapply(tukey.res.imp, MetaboAnalystR:::parseTukey,
                                     cut.off = thresh))
        posthoc_p.adj.sig.imp <-
          unlist(lapply(tukey.res.imp, function (tukey, cut.off) {
            inx <- tukey$cls[, "p adj"] <= cut.off
            paste(signif(data.frame(tukey["cls"])[, 4][inx], 5), collapse = "; ")
          },
          cut.off = thresh))

        tukey.res.all <- lapply(aov.res, TukeyHSD, conf.level = 1 -
                                  thresh)
        cmp.res.all <- unlist(lapply(tukey.res.all, MetaboAnalystR:::parseTukey,
                                     cut.off = thresh))
        posthoc_p.adj.all <- unlist(lapply(tukey.res.all, function (tukey)
          paste(signif(data.frame(tukey["cls"])[, 4], 5), collapse = "; ")))
        post.nm = "Tukey's HSD"
      }
      else {
        fisher.res <- lapply(aov.imp, MetaboAnalystR:::FisherLSD, thresh)
        cmp.res.imp <- unlist(lapply(fisher.res, MetaboAnalystR:::parseFisher,
                                     cut.off = thresh))
        posthoc_p.adj.sig.imp <-
          unlist(lapply(fisher.res, function (fisher, cut.off) {
            inx <- fisher[, "pvalue"] <= cut.off
            paste(signif(data.frame(fisher)[, 2][inx], 5), collapse = "; ")
          },
          cut.off = thresh))
        fisher.res.all <- lapply(aov.res, MetaboAnalystR:::FisherLSD, thresh)
        posthoc_p.adj.all <- unlist(lapply(fisher.res, function (fisher)
          paste(signif(data.frame(fisher)[, 2], 5), collapse = "; ")))
        post.nm = "Fisher's LSD"
      }
      sig.mat <- data.frame(
        signif(sig.f, 5),
        signif(sig.p, 5),
        signif(-log10(sig.p), 5),
        signif(fdr.p, 5),
        cmp.res.imp,
        posthoc_p.adj.sig.imp
      )
      rownames(sig.mat) <- names(sig.p)
      colnames(sig.mat) <- c("f.value", "p.value",
                             "-log10(p)", "FDR", post.nm, "posthoc_p.adj")
      ord.inx <- order(sig.p, decreasing = FALSE)
      sig.mat <- sig.mat[ord.inx, , drop = F]
      fileName <- "anova_posthoc.csv"
    }
  }
  if(!silent){
    print(paste(c("ANOVA analysis: A total of", sum(inx.imp), "significant features were found."),
                collapse = " "))
  }
  if (sig.num > 0) {
    res <- 1
    fast.write.csv(sig.mat, file = fileName)
    conditions <- levels(mSetObj$dataSet$cls)
    cntrst <- rev(apply(utils::combn(rev(conditions), 2), 2, paste,
                        collapse = "-"))
    aov <- list(aov.nm = aov.nm, sig.num = sig.num, sig.nm = fileName,
                raw.thresh = thresh, thresh = -log10(thresh), p.value = p.value, p.adj = p.adjust(p.value, "fdr"),
                p.log = -log10(p.value), inx.imp = inx.imp, post.hoc = post.hoc,
                sig.mat = sig.mat, post.hoc.all = cbind(contrasts=paste(cntrst, collapse= "; "), data.frame(posthoc_signif=cmp.res.all, p.adj=posthoc_p.adj.all)))
  }
  else {
    res <- 0
    aov <- list(aov.nm = aov.nm, sig.num = sig.num, raw.thresh = thresh,
                thresh = -log10(thresh), p.value = p.value, p.log = -log10(p.value),
                inx.imp = inx.imp)
  }
  mSetObj$analSet$aov <- aov
  return(mSetObj)
}

#' Fold change analysis, unpaired
#'
#' Perform fold change analysis, method can be mean or median
#'
#' @param mSetObj Input name of the created mSet object. (see \code{\link[VisomX]{met.read_data}})
#' @param log2fc.thresh (Numeric) Log2(fold-change) threshold that is considered relevant.
#' @param grp1 (Character) Enter name of the first group for the contrast \code{grp1 vs. grp2}. If both group arguments are empty, the first two names in the list of groups are selected.
#' @param grp2 (Character) Enter name of the second group for the contrast \code{grp1 vs. grp2}. If both group arguments are empty, the first two names in the list of groups are selected.
#' @param paired (Logical) Are the data in both groups paired (\code{TRUE}) or not (\code{FALSE}).
#' @return The input mSet object with results of fold-change analysis added at mSetObj$analSet$fc$\code{grp1}_vs_\code{grp2}.
#' @export
met.FC.Anal <- function (mSetObj = NA, log2fc.thresh = 1, grp1 = NULL, grp2 = NULL, paired = FALSE)
{

  max.thresh = log2fc.thresh
  min.thresh = -log2fc.thresh
  if(is.null(grp1) && is.null(grp2)){
    grp1 <- unique(levels(mSetObj$dataSet$cls))[1]
    grp2 <- unique(levels(mSetObj$dataSet$cls))[2]
  }
  res <- met.GetFC(mSetObj, paired, grp1, grp2)
  fc.all <- res$fc.all
  fc.log <- res$fc.log
  inx.up <- fc.log > max.thresh
  inx.down <- fc.log < min.thresh
  names(inx.up) <- names(inx.down) <- names(fc.all)
  imp.inx <- inx.up | inx.down
  sig.mat <- cbind(fc.all[imp.inx, drop = F], fc.log[imp.inx,
                                                     drop = F])
  colnames(sig.mat) <- c("Fold Change", "log2(FC)")
  inx.ord <- order(abs(sig.mat[, 2]), decreasing = T)
  sig.mat <- sig.mat[inx.ord, , drop = F]
  fileName <- "fold_change.csv"
  fast.write.csv(sig.mat, file = fileName)
  if(!exists("fc", mSetObj$analSet)){
    mSetObj$analSet$fc <- list()
  }
  assign(paste0("fc_",grp1,"_vs_",grp2), list(paired = paired, raw.thresh = log2fc.thresh,
                                              max.thresh = max.thresh, min.thresh = min.thresh, fc.all = fc.all,
                                              fc.log = fc.log, inx.up = inx.up, inx.down = inx.down,
                                              inx.imp = imp.inx, sig.mat = sig.mat, grp1 = grp1, grp2 = grp2))
  mSetObj$analSet$fc[[paste0(grp1,"_vs_",grp2)]] <- get(paste0("fc_",grp1,"_vs_",grp2))

  return(mSetObj)
}

#' Used by higher functions to calculate fold change
#'
#' Utility method to calculate FC, used in higher function
#'
#' @param mSetObj Input name of the created mSet object. (see \code{\link[VisomX]{met.read_data}})
#' @param paired (Logical) Are the data in both groups paired (\code{TRUE}) or not (\code{FALSE}).
#' @param grp1 (Character) Enter name of the first group for the contrast \code{grp1 vs. grp2}. If both group arguments are empty, the first two names in the list of groups are selected.
#' @param grp2 (Character) Enter name of the second group for the contrast \code{grp1 vs. grp2}. If both group arguments are empty, the first two names in the list of groups are selected.
#' @return The input mSet object with results of fold-change analysis added at mSetObj$analSet$fc$\code{grp1}_vs_\code{grp2}.
#' @export
met.GetFC <- function (mSetObj = NA, paired = FALSE, grp1, grp2)
{
  assertthat::assert_that(grp1 %in% levels(mSetObj$dataSet$cls),
                          grp2 %in% levels(mSetObj$dataSet$cls),
                          msg = paste0("'",grp1, "' or '", grp2, "' are not valid conditions in the dataset. Valid conditions are:\n",
                                       paste(levels(mSetObj$dataSet$cls), collapse = ", ")))

  if (paired) {
    if (mSetObj$dataSet$combined.method) {
      data <- mSetObj$dataSet$norm
    }
    else {
      if(!is.null(mSetObj$dataSet$row_norm)){
        row.norm <- mSetObj$dataSet$row_norm
      } else {
        row.norm <-  qs::qread("row_norm.qs")
      }
      data <- log2(row.norm)
    }
    G1 <- data[which(mSetObj$dataSet$cls == levels(mSetObj$dataSet$cls)[charmatch(grp1, levels(mSetObj$dataSet$cls) ) ]), ]
    G2 <- data[which(mSetObj$dataSet$cls == levels(mSetObj$dataSet$cls)[charmatch(grp2, levels(mSetObj$dataSet$cls) )]), ]
    fc.mat <- G1 - G2
    fc.log <- apply(fc.mat, 2, mean)
    fc.all <- signif(2^fc.log, 5)
  }
  else {
    data <- NULL
    if (mSetObj$dataSet$combined.method) {
      data <- mSetObj$dataSet$norm
      m1 <- colMeans(data[which(mSetObj$dataSet$cls ==
                                  levels(mSetObj$dataSet$cls)[charmatch(grp1, levels(mSetObj$dataSet$cls) ) ]), ])
      m2 <- colMeans(data[which(mSetObj$dataSet$cls ==
                                  levels(mSetObj$dataSet$cls)[charmatch(grp2, levels(mSetObj$dataSet$cls) ) ]), ])
      fc.log <- signif(m1 - m2, 5)
      fc.all <- signif(2^fc.log, 5)
    }
    else {
      if(!is.null(mSetObj$dataSet$row_norm)){
        data <- mSetObj$dataSet$row_norm
      } else {
        data <-  qs::qread("row_norm.qs")
      }
      m1 <- colMeans(data[which(mSetObj$dataSet$cls ==
                                  levels(mSetObj$dataSet$cls)[charmatch(grp1, levels(mSetObj$dataSet$cls) ) ]), ])
      m2 <- colMeans(data[which(mSetObj$dataSet$cls ==
                                  levels(mSetObj$dataSet$cls)[charmatch(grp2, levels(mSetObj$dataSet$cls) ) ]), ])

      ratio <- m1/m2
      fc.all <- signif(ratio, 5)
      ratio[ratio < 0] <- 0
      fc.log <- signif(log2(ratio), 5)
      fc.log[is.infinite(fc.log) & fc.log < 0] <- -99
      fc.log[is.infinite(fc.log) & fc.log > 0] <- 99
    }
  }
  names(fc.all) <- names(fc.log) <- colnames(data)
  return(list(fc.all = fc.all, fc.log = fc.log))
}

#' Methods for non-specific filtering of variables
#'
#' \code{met.FilterVariable} filters non-informative variables (i.e., features with very small values, near-constant values, or low repeatability) from the data set, dependent on the user-specified method for filtering. The function applies a filtering method, ranks the variables within the data set, and removes variables based on its rank. The final data set should contain no more than than 5000 variables for effective computing. If more features are present, the IQR filter will be applied to keep only a number of 5000, even if \code{filter = "none"}. Data filtering is performed as part of \code{\link[VisomX]{met.read_data}}.
#'
#' @param mSetObj Input name of the created mSet object (see \code{\link[MetaboAnalystR]{InitDataObjects}} and \code{\link[MetaboAnalystR]{Read.TextData}}).
#' @param filter (Character) Select the filter option. \code{"none"} applies no filtering based on the following ranking criteria:
#' \itemize{
#'  \item \code{"rsd"} filters features with low relative standard deviation across the data set.
#'  \item \code{"nrsd}" is the non-parametric relative standard deviation.
#'  \item \code{"mean"} filters features with low mean intensity value across the data set.
#'  \item \code{"median"} filters features with low median intensity value across the data set.
#'  \item \code{"sd"} filters features with low absolute standard deviation across the data set.
#'  \item \code{"mad"} filters features with low median absolute deviation across the data set.
#'  \item \code{"iqr"} filters features with a low inter-quartile range across the data set.
#' @param qcFilter (Logical) Filter the variables based on the relative standard deviation of features in QC samples (\code{TRUE}), or not (\code{FALSE}). This filter can be applied in addition to other, unspecific filtering methods.
#' @param qc.rsd (Numeric) Define the relative standard deviation cut-off in %. Variables with a RSD greater than this number will be removed from the data set. It is only necessary to specify this argument if \code{qcFilter} is \code{TRUE}. Otherwise, it will not be used in the function.
#' @param remain.num (Numerical) Enter the number of variables to keep in your data set. If \code{NULL}, the following empirical rules are applied during data filtering with the methods specified in \code{filter = ""}:
#' \itemize{
#'   \item \strong{Less than 250 variables:} 5% will be filtered
#'   \item \strong{250 - 500 variables:} 10% will be filtered
#'   \item \strong{500 - 1000 variables:} 25% will be filtered
#'   \item \strong{More than 1000 variables:} 40% will be filtered
#' }
#' @param all.rsd (Numeric or \code{NULL}) Apply a filter based on the in-group relative standard deviation (RSD, in %) or not \code{NULL}. Therefore, the RSD of every feature is calculated for every group in the data set. If the RSD of a variable in any group exceeds the indicated threshold, it is removed from the data set. This filter can be applied in addition to other filtering methods and is especially useful to perform on data with technical replicates.
#' @return The input mSet object with filtered data added at mSetObj$dataSet$filt.
#' @export
met.FilterVariable <- function (mSetObj = NA, filter = "none", qcFilter="F", qc.rsd=0.25, remain.num = NULL, all.rsd = NULL)
{
  mSetObj$dataSet$filt <- NULL
  if (is.null(mSetObj$dataSet$proc)) {
    int.mat <- as.matrix(mSetObj[["dataSet"]][["data_proc"]])
  }
  else {
    int.mat <- as.matrix(mSetObj$dataSet$proc)
  }
  cls <- mSetObj$dataSet$proc.cls
  mSetObj$dataSet$filt.cls <- cls
  if (substring(mSetObj$dataSet$format, 4, 5) == "ts") {
    mSetObj$dataSet$filt.facA <- mSetObj$dataSet$proc.facA
    mSetObj$dataSet$filt.facB <- mSetObj$dataSet$proc.facB
  }
  msg <- ""
  if (qcFilter == "T") {
    rsd <- qc.rsd/100
    qc.hits <- tolower(as.character(cls)) %in% "qc"
    if (sum(qc.hits) > 1) {
      qc.mat <- int.mat[qc.hits, ]
      sds <- apply(qc.mat, 2, sd, na.rm = T)
      mns <- apply(qc.mat, 2, mean, na.rm = T)
      rsd.vals <- abs(sds/mns)
      gd.inx <- rsd.vals < rsd
      int.mat <- int.mat[, gd.inx]
      if (mSetObj$analSet$type == "mummichog") {
        msg <- paste("Removed ", sum(!gd.inx),
                     " features based on QC RSD values. QC samples are excluded from downstream functional analysis.")
      }
      else {
        msg <- paste("Removed ", sum(!gd.inx),
                     " features based on QC RSD values. QC samples are still kept. You can remove them later.")
      }
    }
    else if (sum(qc.hits) > 0) {
      AddErrMsg("RSD requires at least 2 QC samples, and only non-QC based filtering can be applied.")
      return(0)
    }
    else {
      AddErrMsg("No QC Samples (with class label: QC) found.  Please use non-QC based filtering.")
      return(0)
    }
  }
  filt.res <- met.PerformFeatureFilter(int.mat, mSetObj, filter, NULL, mSetObj$analSet$type, all.rsd = all.rsd)
  mSetObj$dataSet$filt <- filt.res$data
  msg <- paste(msg, filt.res$msg)
  mSetObj$msgSet$filter.msg <- msg
  return(mSetObj)
}

#' Used by higher functions for non-specific filtering of variables
#'
#' Utility method to perform data filtering, used in higher function
#'
#' @param int.mat Data matrix with samples in rows and features in columns, generated by higher function.
#' @param mSetObj Input name of the created mSet object (see \code{\link[MetaboAnalystR]{InitDataObjects}} and \code{\link[MetaboAnalystR]{Read.TextData}}).
#' @param filter (Character) Select the filter option. \code{"none"} applies no filtering based on the following ranking criteria:
#' \itemize{
#'  \item \code{"rsd"} filters features with low relative standard deviation across the data set.
#'  \item \code{"nrsd}" is the non-parametric relative standard deviation.
#'  \item \code{"mean"} filters features with low mean intensity value across the data set.
#'  \item \code{"median"} filters features with low median intensity value across the data set.
#'  \item \code{"sd"} filters features with low absolute standard deviation across the data set.
#'  \item \code{"mad"} filters features with low median absolute deviation across the data set.
#'  \item \code{"iqr"} filters features with a low inter-quartile range across the data set.
#' @param remain.num (Numerical) Enter the number of variables to keep in your data set. If \code{NULL}, the following empirical rules are applied during data filtering with the methods specified in \code{filter = ""}:
#' \itemize{
#'   \item \strong{Less than 250 variables:} 5% will be filtered
#'   \item \strong{250 - 500 variables:} 10% will be filtered
#'   \item \strong{500 - 1000 variables:} 25% will be filtered
#'   \item \strong{More than 1000 variables:} 40% will be filtered
#' }
#' @param anal.type (Character) Type of analysis. Extracted from mSetObj by higher function at mSetObj$analSet$type.
#' @param all.rsd (Numeric or \code{NULL}) Apply a filter based on the in-group relative standard deviation (RSD, in %) or not \code{NULL}. Therefore, the RSD of every feature is calculated for every group in the data set. If the RSD of a variable in any group exceeds the indicated threshold, it is removed from the data set. This filter can be applied in addition to other filtering methods and is especially useful to perform on data with technical replicates.
#' @return A list with filtered data and a message to inform about the chosen filtering conditions
met.PerformFeatureFilter <- function (int.mat, mSetObj, filter, remain.num = remain.num, anal.type = NULL, all.rsd = NULL)
{
  feat.num <- ncol(int.mat)
  feat.nms <- colnames(int.mat)
  nm <- NULL
  msg <- ""
  if (filter == "none" && feat.num < 5000) {
    remain <- rep(TRUE, feat.num)
    msg <- paste(msg, "No unspecific filtering was applied")
  }
  else if(!is.null(all.rsd)){
    all.rsd <- all.rsd/100
    smaller_0.05quant <- int.mat[,apply(int.mat<unname(quantile(int.mat, probs=seq(0, 1 , 1/20))[2]), 2, any)]
    vec.smaller_0.05quant <- replace(smaller_0.05quant[1,], !is.null(smaller_0.05quant[1,]), 0)
    greater_0.05quant <- int.mat[,apply(int.mat<unname(quantile(int.mat, probs=seq(0, 1 , 1/20))[2]), 2, any)!=TRUE]
    ls.sd <-
      lapply(1:length(levels(mSetObj[["dataSet"]][["cls"]])), function (x)
        apply(greater_0.05quant[grep(levels(mSetObj[["dataSet"]][["cls"]])[x], as.character(mSetObj[["dataSet"]][["cls"]])),], 2, sd, na.rm = T))
    ls.mns <-
      lapply(1:length(levels(mSetObj[["dataSet"]][["cls"]])), function (x)
        apply(greater_0.05quant[grep(levels(mSetObj[["dataSet"]][["cls"]])[x], as.character(mSetObj[["dataSet"]][["cls"]])), ], 2, mean, na.rm = T))
    ls.rsd <- mapply("/", ls.sd, ls.mns, SIMPLIFY = FALSE)
    mat.rsd <- do.call(rbind, ls.rsd)
    filter.val <-
      c(vec.smaller_0.05quant, apply(mat.rsd, 2, max, na.rm = T))
    names(filter.val) <- str_replace(names(filter.val), "^X", "")
    filter.val<- filter.val[order(names(filter.val))]
    names(filter.val) <- gsub("^([0-9]+)", "X\\1", names(filter.val))
    remain <- filter.val <= all.rsd
    nm <- paste0("Relative standard deviation (", length(remain)-length(remain[remain]), " features with RSD > ", all.rsd, "). Only metabolites with values above the 5% quantile in all samples were considered for filtering.")
    msg <- paste0(msg, "Feature filtering based on: ",
                  nm)
  }
  else {
    if (filter == "rsd") {
      sds <- apply(int.mat, 2, sd, na.rm = T)
      mns <- apply(int.mat, 2, mean, na.rm = T)
      filter.val <- abs(sds/mns)
      nm <- "Relative standard deviation"
    }
    else if (filter == "nrsd") {
      mads <- apply(int.mat, 2, mad, na.rm = T)
      meds <- apply(int.mat, 2, median, na.rm = T)
      filter.val <- abs(mads/meds)
      nm <- "Non-paramatric relative standard deviation"
    }
    else if (filter == "mean") {
      filter.val <- apply(int.mat, 2, mean, na.rm = T)
      nm <- "mean"
    }
    else if (filter == "sd") {
      filter.val <- apply(int.mat, 2, sd, na.rm = T)
      nm <- "standard deviation"
    }
    else if (filter == "mad") {
      filter.val <- apply(int.mat, 2, mad, na.rm = T)
      nm <- "Median absolute deviation"
    }
    else if (filter == "median") {
      filter.val <- apply(int.mat, 2, median, na.rm = T)
      nm <- "median"
    }
    else if (filter == "iqr"){
      filter.val <- apply(int.mat, 2, IQR, na.rm = T)
      nm <- "Interquantile Range"
    }
    else {
      stop(paste0("\"", filter, "\" is not a valid filtering condition. See ?met.FilterVariable for more information about suitable options."))
    }
    rk <- rank(-filter.val, ties.method = "random")
    if (is.null(remain.num)) {
      if (feat.num < 250) {
        remain <- rk < feat.num * 0.95
        msg <- paste(msg, "Further feature filtering based on",
                     nm)
      }
      else if (feat.num < 500) {
        remain <- rk < feat.num * 0.9
        msg <- paste(msg, "Further feature filtering based on",
                     nm)
      }
      else if (feat.num < 1000) {
        remain <- rk < feat.num * 0.75
        msg <- paste(msg, "Further feature filtering based on",
                     nm)
      }
      else {
        remain <- rk < feat.num * 0.6
        msg <- paste(msg, "Further feature filtering based on",
                     nm)
        if (anal.type == "mummichog") {
          max.allow <- 7500
        }
        else if (anal.type == "power" || anal.type ==
                 "ts") {
          max.allow <- 5000
        }
        else {
          max.allow <- 2500
        }
        if (sum(remain) > max.allow) {
          remain <- rk < max.allow
          msg <- paste(msg, paste("Reduced to",
                                  max.allow, "features based on", nm))
        }
      }
    }
    else {
      remain <- rk < remain.num
    }
  }
  print(msg)
  return(list(data = int.mat[, remain], msg = msg))
}

#' Perform T-test analysis
#'
#' \code{met.Ttests.Anal} performs the Student's t-test on chosen sample groups. For large data set (> 1000 variables), both the paired information and the group variance will be ignored, and the default parameters will be used for t-tests to save computational time. If you choose non-parametric tests (Wilcoxon rank-sum test), the group variance will be ignored.
#'
#' @param mSetObj Input name of the created mSet object (see \code{\link[VisomX]{met.read_data}}).
#' @param grp1 (Character) Enter name of the first group for the contrast \code{grp1 vs. grp2}. If both group arguments are empty, the first two names in the list of groups are selected.
#' @param grp2 (Character) Enter name of the second group for the contrast \code{grp1 vs. grp2}. If both group arguments are empty, the first two names in the list of groups are selected.
#' @param nonpar (Logical) Use a non-parametric test (\code{TRUE}) or not (\code{FALSE}).
#' @param threshp (Numerical) Enter the adjusted p-value (FDR) cutoff.
#' @param paired (Logical) Is the data paired (\code{TRUE}) or not (\code{FALSE}).
#' @param equal.var (Logical) Is the group variance equal (\code{TRUE}) or not (\code{FALSE}).
#' @param pvalType (Character) p value type used to apply significances based on the chosen threshold \code{threshp}. Can be \code{"fdr"} for adjusted p values or \code{"raw"} for raw p values.
#' @param all_results (Logical) Create a CSV file with T-Test results for all compounds (\code{TRUE}) or not (\code{FALSE}).
#' @param silent (Logical) Suppress message about the number of significant features found in the console (\code{TRUE}) or not (\code{FALSE}).
#' @return The input mSet object with T-Test results added at mSetObj$analSet$tt$grp1_vs_grp2.
#' @export
$met.Ttests.Anal
function (mSetObj = NA, grp1, grp2, nonpar = FALSE, threshp = 0.05, paired = FALSE,
          equal.var = TRUE, pvalType = "fdr", all_results = TRUE, silent = FALSE)
{
  if(!is.null(grp1) && !is.null(grp2)){
    assertthat::assert_that(grp1 %in% levels(mSetObj$dataSet$cls),
                            grp2 %in% levels(mSetObj$dataSet$cls),
                            msg = paste0("'",grp1, "' or '", grp2, "' are not valid conditions in the dataset. Valid conditions are:\n",
                                         paste(levels(mSetObj$dataSet$cls), collapse = ", ")))
  }
  if(is.null(grp1) && is.null(grp2)){
    grp1 <- unique(levels(mSetObj$dataSet$cls))[1]
    grp2 <- unique(levels(mSetObj$dataSet$cls))[2]
  }
  res <- met.GetTtestRes(mSetObj, grp1, grp2, paired, equal.var, nonpar)
  t.stat <- res[, 1]
  p.value <- res[, 2]
  names(t.stat) <- names(p.value) <- colnames(mSetObj$dataSet$norm)
  p.log <- -log10(p.value)
  fdr.p <- p.adjust(p.value, "fdr")
  if (all_results == TRUE) {
    all.mat <- data.frame(signif(t.stat, 5), signif(p.value,
                                                    5), signif(p.log, 5), signif(fdr.p, 5))
    if (nonpar) {
      tt.nm = "Wilcoxon Rank Test"
      file.nm <- "wilcox_rank_all.csv"
      colnames(all.mat) <- c("V", "p.value",
                             "-log10(p)", "FDR")
    }
    else {
      tt.nm = "T-Tests"
      file.nm <- paste0("t_test_all_", grp1, "_vs_", grp2, ".csv")
      colnames(all.mat) <- c("t.stat", "p.value",
                             "-log10(p)", "FDR")
    }
    fast.write.csv(all.mat, file = file.nm)
  }
  if (pvalType == "fdr") {
    inx.imp <- fdr.p <= threshp
  }
  else {
    inx.imp <- p.value <= threshp
  }
  sig.num <- sum(inx.imp, na.rm = TRUE)
  if(silent==FALSE){
    print(paste("T-test: A total of", sig.num, "significant features were found for contrast:", grp1, "vs.", grp2))
  }
  if (sig.num > 0) {
    sig.t <- t.stat[inx.imp]
    sig.p <- p.value[inx.imp]
    lod <- -log10(sig.p)
    sig.q <- fdr.p[inx.imp]
    sig.mat <- cbind(sig.t, sig.p, lod, sig.q)
    colnames(sig.mat) <- c("t.stat", "p.value",
                           "-log10(p)", "FDR")
    ord.inx <- order(sig.p)
    sig.mat <- sig.mat[ord.inx, , drop = F]
    sig.mat <- signif(sig.mat, 5)
    if (nonpar) {
      tt.nm = "Wilcoxon Rank Test"
      file.nm <- "wilcox_rank.csv"
      colnames(sig.mat) <- c("V", "p.value",
                             "-log10(p)", "FDR")
    }
    else {
      tt.nm = "T-Tests"
      file.nm <- paste0("t_test_", grp1, "_vs_", grp2, ".csv")
      colnames(sig.mat) <- c("t.stat", "p.value",
                             "-log10(p)", "FDR")
    }
    fast.write.csv(sig.mat, file = file.nm)
    tt <- list(tt.nm = tt.nm, sig.nm = file.nm, sig.num = sig.num,
               paired = paired, pval.type = pvalType, raw.thresh = threshp,
               t.score = t.stat, p.value = p.value, p.log = p.log,
               thresh = -log10(threshp), inx.imp = inx.imp, sig.mat = sig.mat)
  }
  else {
    tt <- list(sig.num = sig.num, paired = paired, pval.type = pvalType,
               raw.thresh = threshp, t.score = t.stat, p.value = p.value,
               p.log = p.log, thresh = -log10(threshp), inx.imp = inx.imp)
  }
  if(!exists("tt", mSetObj$analSet)){
    mSetObj$analSet$tt <- list()
  }
  mSetObj$analSet$tt[[paste0(grp1,"_vs_",grp2)]] <- tt
  return(mSetObj)
}

#' Used by higher functions to retrieve T-test p-values
#'
#' Utility method to to get p values via T-test.
#'
#' @param mSetObj Input name of the created mSet object (see \code{\link[VisomX]{met.read_data}}).
#' @param grp1 (Character) Enter name of the first group for the contrast \code{grp1 vs. grp2}.
#' @param grp2 (Character) Enter name of the second group for the contrast \code{grp1 vs. grp2}.
#' @param paired (Logical) Is the data paired (\code{TRUE}) or not (\code{FALSE}).
#' @param equal.var (Logical) Is the group variance equal (\code{TRUE}) or not (\code{FALSE}).
#' @param nonpar (Logical) Use a non-parametric test (\code{TRUE}) or not (\code{FALSE}).
#' @return A data frame with T-test results.
met.GetTtestRes <- function (mSetObj = NA, grp1, grp2, paired = FALSE, equal.var = TRUE, nonpar = F)
{

  inx1 <- which(mSetObj$dataSet$cls == levels(mSetObj$dataSet$cls)[charmatch(grp1, levels(mSetObj$dataSet$cls) ) ])
  inx2 <- which(mSetObj$dataSet$cls == levels(mSetObj$dataSet$cls)[charmatch(grp2, levels(mSetObj$dataSet$cls) ) ])
  if (length(inx1) == 1 || length(inx2) == 1) {
    equal.var <- TRUE # overwrite option if one does not have enough replicates
  }
  data <- as.matrix(mSetObj$dataSet$norm)
  if (!exists("mem.tt")) {
    mem.tt <<- memoise::memoise(MetaboAnalystR:::.get.ttest.res)
  }
  return(mem.tt(data, inx1, inx2, paired, equal.var, nonpar))
}

#' Impute missing variables
#'
#' Replace missing variables via a chosen method. Data need to be re-calibrated after this step, including \code{\link[VisomX::met.PerformFeatureFilter]{filtering}} as well as \code{\link[VisomX::met.normalize]{normalization, transformation, and scaling}}. Data imputation is performed as part of \code{\link[VisomX]{met.read_data}}.
#'
#' @param mSetObj Input name of the created mSet object ((see \code{\link[VisomX]{$met.initialize}} and \code{\link[MetaboAnalystR]{Read.TextData}}).
#' @param method (Character) Is the data paired (\code{TRUE}) or not (\code{FALSE}).
#' @return The input mSet object with imputed data at mSetObj$dataSet$data_proc.
#' @export
met.impute <- function (mSetObj = NA, method = "min") {

  mSetObj$dataSet$filt <- mSetObj$dataSet$edit <- NULL
  if(!is.null(mSetObj$dataSet$preproc)){
    int.mat <- mSetObj$dataSet$preproc
  } else {
    int.mat <- qs::qread("preproc.qs")
  }
  new.mat <- NULL
  msg <- mSetObj$msgSet$replace.msg
  if (method == "exclude") {
    good.inx <- apply(is.na(int.mat), 2, sum) == 0
    new.mat <- int.mat[, good.inx, drop = FALSE]
    msg <- c(msg, "Variables with missing values were excluded.")
  }
  else if (method == "min") {
    new.mat <- MetaboAnalystR:::ReplaceMissingByLoD(int.mat)
    msg <- c(msg, "Missing variables were replaced by LoDs (1/5 of the min positive value for each variable)")
  }
  else if (method == "colmin") {
    new.mat <- apply(int.mat, 2, function(x) {
      if (sum(is.na(x)) > 0) {
        x[is.na(x)] <- min(x, na.rm = T)/2
      }
      x
    })
    msg <- c(msg, "Missing variables were replaced by 1/2 of min values for each feature column.")
  }
  else if (method == "mean") {
    new.mat <- apply(int.mat, 2, function(x) {
      if (sum(is.na(x)) > 0) {
        x[is.na(x)] <- mean(x, na.rm = T)
      }
      x
    })
    msg <- c(msg, "Missing variables were replaced with the mean value for each feature column.")
  }
  else if (method == "median") {
    new.mat <- apply(int.mat, 2, function(x) {
      if (sum(is.na(x)) > 0) {
        x[is.na(x)] <- median(x, na.rm = T)
      }
      x
    })
    msg <- c(msg, "Missing variables were replaced with the median for each feature column.")
  }
  else {
    if (method == "knn_var") {
      new.mat <- t(impute::impute.knn(t(int.mat))$data)
    }
    else if (method == "knn_smp") {
      new.mat <- impute::impute.knn(data.matrix(int.mat))$data
    }
    else {
      if (method == "bpca") {
        new.mat <- pcaMethods::pca(int.mat, nPcs = 5,
                                   method = "bpca", center = T)@completeObs
      }
      else if (method == "ppca") {
        new.mat <- pcaMethods::pca(int.mat, nPcs = 5,
                                   method = "ppca", center = T)@completeObs
      }
      else if (method == "svdImpute") {
        new.mat <- pcaMethods::pca(int.mat, nPcs = 5,
                                   method = "svdImpute", center = T)@completeObs
      }
    }
    msg <- c(msg, paste("Missing variables were imputated using",
                        toupper(method)))
  }
  mSetObj$dataSet$proc.feat.num <- ncol(new.mat)
  mSetObj$dataSet$data_proc <- as.data.frame(new.mat)
  message("Writing data set with imputed values to 'data_proc.qs'")
  qs::qsave(as.data.frame(new.mat), file = "data_proc.qs")
  mSetObj$msgSet$replace.msg <- msg
  print(mSetObj$msgSet$replace.msg)
  return(mSetObj)
}

$met.initialize
function (data.type, anal.type, paired = FALSE)
{
  dataSet <- list()
  dataSet$type <- data.type
  dataSet$design.type <- "regular"
  dataSet$cls.type <- "disc"
  dataSet$format <- "rowu"
  dataSet$paired <- paired
  analSet <- list()
  analSet$type <- anal.type
  mSetObj <- list()
  mSetObj$dataSet <- dataSet
  mSetObj$analSet <- analSet
  mSetObj$imgSet <- list()
  mSetObj$msgSet <- list()
  mSetObj$msgSet$msg.vec <- vector(mode = "character")
  mSetObj$cmdSet <- vector(mode = "character")
  if (anal.type == "mummichog") {
    mSetObj$paramSet$mumRT <- NA
    mSetObj$paramSet$mumRT.type <- NA
    mSetObj$paramSet$version <- NA
    mSetObj$paramSet$mumDataContainsPval <- 1
    mSetObj$paramSet$mode <- NA
    mSetObj$paramSet$adducts <- NA
    mSetObj$paramSet$peakFormat <- "mpt"
  }
  else if (anal.type == "metapaths") {
    paramSet <- list()
    paramSet$mumRT <- NA
    paramSet$mumRT.type <- NA
    paramSet$version <- NA
    paramSet$mumDataContainsPval <- 1
    paramSet$mode <- NA
    paramSet$adducts <- NA
    paramSet$peakFormat <- "mpt"
    paramSet$metaNum <- 0
    mSetObj$paramSet <- paramSet
    dataNMs <- names(mSetObj)[grepl("MetaData", names(mSetObj))]
    if (length(dataNMs) > 0) {
      for (n in dataNMs) {
        mSetObj[[n]] <- NULL
      }
    }
  }
  module.count <<- 0
  smpdbpw.count <<- 0
  mdata.all <<- list()
  mdata.siggenes <<- vector("list")
  meta.selected <<- TRUE
  anal.type <<- anal.type
  Cairo::CairoFonts(regular = "Arial:style=Regular",
                    bold = "Arial:style=Bold", italic = "Arial:style=Italic",
                    bolditalic = "Arial:style=Bold Italic", symbol = "Symbol")
  print("MetaboAnalyst R objects initialized ...")
  return(mSetObj)
}

$met.normalize
function (mSetObj = NA, rowNorm, transNorm, scaleNorm, ref = NULL, norm.vec = NULL,
          ratio = FALSE, ratioNum = 20)
{
  if(!is.null(mSetObj$dataSet$prenorm)){
    data <- mSetObj$dataSet$prenorm
  } else {
    data <- qs::qread("prenorm.qs")
  }
  cls <- mSetObj$dataSet$prenorm.cls
  if (substring(mSetObj$dataSet$format, 4, 5) == "ts") {
    if (is.null(mSetObj$dataSet$prenorm.facA)) {
      nfacA <- mSetObj$dataSet$facA
      nfacB <- mSetObj$dataSet$facB
    }
    else {
      nfacA <- mSetObj$dataSet$prenorm.facA
      nfacB <- mSetObj$dataSet$prenorm.facB
    }
    mSetObj$dataSet$facA <- nfacA
    mSetObj$dataSet$facB <- nfacB
    if (mSetObj$dataSet$design.type == "time" | mSetObj$dataSet$design.type ==
        "time0") {
      if (tolower(mSetObj$dataSet$facA.lbl) == "time") {
        time.fac <- nfacA
        exp.fac <- nfacB
      }
      else {
        time.fac <- nfacB
        exp.fac <- nfacA
      }
      lvls <- levels(time.fac)
      time.points <- as.numeric(as.character(lvls))
      ord.lvls <- lvls[order(time.points)]
      time.fac <- ordered(time.fac, levels = ord.lvls)
      mSetObj$dataSet$time.fac <- time.fac
      mSetObj$dataSet$exp.fac <- exp.fac
    }
  }
  colNames <- colnames(data)
  rowNames <- rownames(data)
  if (rowNorm == "QuantileNorm") {
    data <- QuantileNormalize(data)
    varCol <- apply(data, 2, var, na.rm = T)
    constCol <- (varCol == 0 | is.na(varCol))
    constNum <- sum(constCol, na.rm = T)
    if (constNum > 0) {
      print(paste("After quantile normalization",
                  constNum, "features with a constant value were found and deleted."))
      data <- data[, !constCol, drop = FALSE]
      colNames <- colnames(data)
      rowNames <- rownames(data)
    }
    rownm <- "Quantile Normalization"
  }
  else if (rowNorm == "GroupPQN") {
    grp.inx <- cls == ref
    ref.smpl <- apply(data[grp.inx, , drop = FALSE], 2, mean)
    data <- t(apply(data, 1, ProbNorm, ref.smpl))
    rownm <- "Probabilistic Quotient Normalization by a reference group"
  }
  else if (rowNorm == "SamplePQN") {
    ref.smpl <- data[ref, , drop = FALSE]
    data <- t(apply(data, 1, ProbNorm, ref.smpl))
    rownm <- "Probabilistic Quotient Normalization by a reference sample"
  }
  else if (rowNorm == "CompNorm") {
    data <- t(apply(data, 1, CompNorm, ref))
    rownm <- "Normalization by a reference feature"
  }
  else if (rowNorm == "SumNorm") {
    data <- t(apply(data, 1, SumNorm))
    rownm <- "Normalization to constant sum"
  }
  else if (rowNorm == "MedianNorm") {
    data <- t(apply(data, 1, MedianNorm))
    rownm <- "Normalization to sample median"
  }
  else if (rowNorm == "SpecNorm") {
    if (is.null(norm.vec)) {
      norm.vec <- rep(1, nrow(data))
      print("No sample specific information were given, all set to 1.0")
    }
    rownm <- "Normalization by sample-specific factor"
    data <- data/norm.vec
  }
  else {
    rownm <- "N/A"
  }
  rownames(data) <- rowNames
  colnames(data) <- colNames
  if (rowNorm == "CompNorm" && !is.null(ref)) {
    inx <- match(ref, colnames(data))
    data <- data[, -inx, drop = FALSE]
    colNames <- colNames[-inx]
  }
  row.norm <- as.data.frame(MetaboAnalystR:::CleanData(data, T, T))
  qs::qsave(row.norm, file = "row_norm.qs")
  mSetObj$dataSet$row_norm <- row.norm
  if (ratio) {
    min.val <- min(abs(data[data != 0]))/2
    norm.data <- log2((data + sqrt(data^2 + min.val))/2)
    transnm <- "Log2 Normalization"
    ratio.mat <- CalculatePairwiseDiff(norm.data)
    fstats <- Get.Fstat(ratio.mat, cls)
    hit.inx <- rank(-fstats) < ratioNum
    ratio.mat <- ratio.mat[, hit.inx, drop = FALSE]
    data <- cbind(norm.data, ratio.mat)
    colNames <- colnames(data)
    rowNames <- rownames(data)
    mSetObj$dataSet$use.ratio <- TRUE
    mSetObj$dataSet$proc.ratio <- data
  }
  else {
    mSetObj$dataSet$use.ratio <- FALSE
    if (transNorm == "LogNorm") {
      min.val <- min(abs(data[data != 0]))/10
      data <- apply(data, 2, MetaboAnalystR:::LogNorm, min.val)
      transnm <- "Log10 Transformation"
    }
    else if (transNorm == "SrNorm") {
      min.val <- min(abs(data[data != 0]))/10
      data <- apply(data, 2, MetaboAnalystR:::SquareRootNorm, min.val)
      transnm <- "Square Root Transformation"
    }
    else if (transNorm == "CrNorm") {
      norm.data <- abs(data)^(1/3)
      norm.data[data < 0] <- -norm.data[data < 0]
      data <- norm.data
      transnm <- "Cubic Root Transformation"
    }
    else {
      transnm <- "N/A"
    }
  }
  if (scaleNorm == "MeanCenter") {
    data <- apply(data, 2, MetaboAnalystR:::MeanCenter)
    scalenm <- "Mean Centering"
  }
  else if (scaleNorm == "AutoNorm") {
    data <- apply(data, 2, MetaboAnalystR:::AutoNorm)
    scalenm <- "Autoscaling"
  }
  else if (scaleNorm == "ParetoNorm") {
    data <- apply(data, 2, MetaboAnalystR:::ParetoNorm)
    scalenm <- "Pareto Scaling"
  }
  else if (scaleNorm == "RangeNorm") {
    data <- apply(data, 2, MetaboAnalystR:::RangeNorm)
    scalenm <- "Range Scaling"
  }
  else {
    scalenm <- "N/A"
  }
  rownames(data) <- rowNames
  colnames(data) <- colNames
  data <- MetaboAnalystR:::CleanData(data, T, F)
  if (ratio) {
    mSetObj$dataSet$ratio <- MetaboAnalystR:::CleanData(ratio.mat, T, F)
  }
  mSetObj$dataSet$norm <- as.data.frame(data)
  if (substring(mSetObj$dataSet$format, 4, 5) == "ts") {
    if (rownames(mSetObj$dataSet$norm) != rownames(mSetObj$dataSet$meta.info)) {
      print("Metadata and data norm are not synchronized.")
    }
    mSetObj$dataSet$meta.info <- mSetObj$dataSet$meta.info[rownames(data),
    ]
  }
  qs::qsave(mSetObj$dataSet$norm, file = "complete_norm.qs")
  mSetObj$dataSet$cls <- cls
  mSetObj$dataSet$rownorm.method <- rownm
  mSetObj$dataSet$trans.method <- transnm
  mSetObj$dataSet$scale.method <- scalenm
  mSetObj$dataSet$combined.method <- FALSE
  mSetObj$dataSet$norm.all <- NULL
  return(mSetObj)
}



$met.PLS_ImpScatter
function(mSetObj, imgName=NULL, format = "png", dpi = 300, width = NA,
         feat.nm="coef.mean", vip.nm = c("Comp. 1", "Comp. 2"), plot = TRUE, export = FALSE, vip.thresh = 1,
         show.title = FALSE, title){
  if(!is.character(imgName)){
    imgName <- paste0("Plots/PLS_ImpScatter_", feat.nm)
  }
  imgName = paste(imgName, ".", format,
                  sep = "")
  if (is.na(width)) {
    w <- 9
  }
  else if (width == 0) {
    w <- 8
  }
  else {
    w <- width
  }
  h <- w * if_else(show.title==TRUE, 1.1, 1.0)
  mSetObj[["imgSet"]][[paste0("pls.ImpScatter_", feat.nm)]] <- imgName
  df <- data.frame(mSetObj$analSet$plsda$vip.mat, mSetObj$analSet$plsda$coef.mat[rownames(mSetObj$analSet$plsda$vip.mat),], check.names = TRUE)
  feat <- make.names(feat.nm)
  if(length(vip.nm)>1){
    df$vip_avg <- rowMeans(df[,grep(paste0(make.names(vip.nm), collapse="|"), colnames(df))])
    vip <- "vip_avg"
  } else {
    vip <- make.names(vip.nm)
  }
  p <-
    ggplot(df, aes_string(x = vip, y = feat)) +
    geom_point(color = "Black", fill="Gray", shape=21, size = 3.5, alpha = 0.6) +
    ggrepel::geom_text_repel(
      aes(label = rownames(df)),
      box.padding = unit(0.25,
                         "lines"),
      point.padding = unit(0.5, "lines")
    ) +
    geom_vline(xintercept = vip.thresh,
               linetype = "longdash", alpha = 0.4) +
    ylab(if_else(feat.nm=="coef.mean", feat.nm, paste0("Coef. ", feat.nm))) +
    theme_bw(base_size = 18) +
    theme(axis.text=element_text(size = 20),
          axis.title = element_text(size=24,face="bold"))
  if(length(vip.nm)>1){
    p <- p + xlab(paste0("average VIP (", paste(vip.nm, collapse = ", "), ")"))
  } else {
    p <- p + xlab(paste0("VIP - ", vip.nm))
  }
  if(show.title == TRUE){
    p <- p + ggtitle(title) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  }
  if (export == TRUE){
    Cairo::Cairo(file = imgName, unit = "in", dpi = dpi,
                 width = w, height = h, type = format, bg = "white")
    print(p)
    dev.off()
  }
  mSetObj[["imgSet"]][[paste0("pls.ImpScatter_plot_",feat.nm)]] <-  p
  if (plot == TRUE){
    print(p)
  }
  return(mSetObj)
}

  $met.PLSDA.CV
function (mSetObj = NA, methodName = "T", compNum = MetaboAnalystR:::GetDefaultPLSCVComp(mSetObj),
          choice = "Q2", data = "all")
{
  if(data=="anova.sig") {
    plsda.cls <-
      caret::train(
        mSetObj$dataSet$norm[, mSetObj$analSet$aov[["inx.imp"]]],
        mSetObj$dataSet$cls,
        "pls",
        trControl = caret::trainControl(method = ifelse(methodName ==
                                                          "L", "LOOCV", "CV")),
        tuneLength = compNum
      )
  } else {
    plsda.cls <-
      caret::train(
        mSetObj$dataSet$norm,
        mSetObj$dataSet$cls,
        "pls",
        trControl = caret::trainControl(method = ifelse(methodName ==
                                                          "L", "LOOCV", "CV")),
        tuneLength = compNum
      )
  }

  if (mSetObj$analSet$plsr$reg) {
    cls <- cls <- scale(as.numeric(mSetObj$dataSet$cls))[,
                                                         1]
  }
  else {
    cls <- model.matrix(~mSetObj$dataSet$cls - 1)
  }
  datmat <- as.matrix(mSetObj$dataSet$norm)
  plsda.reg <- pls::plsr(cls ~ datmat, method = "oscorespls",
                         ncomp = compNum, validation = ifelse(methodName == "L",
                                                              "LOO", "CV"))
  fit.info <- pls::R2(plsda.reg, estimate = "all")$val[,
                                                       1, ]
  accu <- plsda.cls$results[, 2]
  all.info <- rbind(accu, fit.info[, -1])
  rownames(all.info) <- c("Accuracy", "R2", "Q2")
  if (choice == "Q2") {
    best.num <- which(all.info[3, ] == max(all.info[3, ]))
  }
  else if (choice == "R2") {
    best.num <- which(all.info[2, ] == max(all.info[2, ]))
  }
  else {
    best.num <- which(all.info[1, ] == max(all.info[1, ]))
  }
  coef.mat <- try(caret::varImp(plsda.cls, scale = T)$importance)
  if (class(coef.mat) == "try-error") {
    coef.mat <- NULL
  }
  else {
    if (mSetObj$dataSet$cls.num > 2) {
      coef.mean <- apply(coef.mat, 1, mean)
      coef.mat <- cbind(coef.mean = coef.mean, coef.mat)
    }
    inx.ord <- order(coef.mat[, 1], decreasing = T)
    coef.mat <- data.matrix(coef.mat[inx.ord, , drop = FALSE])
    fast.write.csv(signif(coef.mat, 5), file = "plsda_coef.csv")
  }
  pls <- mSetObj$analSet$plsr
  b <- c(pls$Yloadings)[1:compNum]
  T <- pls$scores[, 1:compNum, drop = FALSE]
  SS <- b^2 * colSums(T^2)
  W <- pls$loading.weights[, 1:compNum, drop = FALSE]
  Wnorm2 <- colSums(W^2)
  SSW <- sweep(W^2, 2, SS/Wnorm2, "*")
  vips <- sqrt(nrow(SSW) * apply(SSW, 1, cumsum)/cumsum(SS))
  if (compNum > 1) {
    vip.mat <- as.matrix(t(vips))
    ord.inx <- order(-abs(vip.mat[, 1]), -abs(vip.mat[, 2]))
  }
  else {
    vip.mat <- as.matrix(vips)
    ord.inx <- order(-abs(vip.mat[, 1]))
  }
  vip.mat <- vip.mat[ord.inx, ]
  colnames(vip.mat) <- paste("Comp.", 1:ncol(vip.mat))
  fast.write.csv(signif(vip.mat, 5), file = "plsda_vip.csv")

  if(data=="anova.sig") {
    plsda.cls <-
      caret::train(
        mSetObj$dataSet$norm[, mSetObj$analSet$aov[["inx.imp"]]],
        mSetObj$dataSet$cls,
        "pls",
        trControl = caret::trainControl(method = ifelse(methodName ==
                                                          "L", "LOOCV", "CV")),
        tuneLength = compNum
      )
  } else {
    plsda.cls <-
      caret::train(
        mSetObj$dataSet$norm,
        mSetObj$dataSet$cls,
        "pls",
        trControl = caret::trainControl(method = ifelse(methodName ==
                                                          "L", "LOOCV", "CV")),
        tuneLength = compNum
      )
  }
  mSetObj$analSet$plsda <- list(best.num = best.num, choice = choice,
                                coef.mat = coef.mat, vip.mat = vip.mat, fit.info = all.info, data.type = data)
  return(mSetObj)
}

  $met.PLSR.Anal
function (mSetObj = NA, reg = FALSE, data = "all")
{
  comp.num <- dim(mSetObj$dataSet$norm)[1] - 1
  if (comp.num > 8) {
    comp.num <- min(dim(mSetObj$dataSet$norm)[2], 8)
  }
  else if (comp.num < 8) {
    AddMsg(paste("Too few features in your data to do PLS-DA analysis! "))
  }
  if (reg) {
    cls <- scale(as.numeric(mSetObj$dataSet$cls))[, 1]
  }
  else {
    cls <- model.matrix(~mSetObj$dataSet$cls - 1)
  }
  if(data=="anova.sig"&&is.null(mSetObj$analSet$aov)) {
    stop("Please perform 'met.ANOVA.Anal' on your mSet before running PLS analysis with 'data = \"anova\"'")
  }
  if(data=="anova.sig") {
    datmat <- as.matrix(mSetObj$dataSet$norm[,mSetObj$analSet$aov[["inx.imp"]]])
  } else {
    datmat <- as.matrix(mSetObj$dataSet$norm)

  }
  mSetObj$analSet$plsr <- pls::plsr(cls ~ datmat, method = "oscorespls",
                                    ncomp = comp.num)
  mSetObj$analSet$plsr$reg <- reg
  mSetObj$analSet$plsr$loading.type <- "all"
  mSetObj$analSet$plsr$data.type <- data
  mSetObj$custom.cmpds <- c()
  fast.write.csv(signif(mSetObj$analSet$plsr$scores, 5), row.names = rownames(mSetObj$dataSet$norm),
                 file = "plsda_score.csv")
  fast.write.csv(signif(mSetObj$analSet$plsr$loadings, 5),
                 file = "plsda_loadings.csv")
  return(mSetObj)
}

  $met.PreparePrenormData
function (mSetObj = NA)
{
  if (!is.null(mSetObj$dataSet$edit)) {
    mydata <- mSetObj$dataSet$edit
    if (!is.null(mSetObj$dataSet$filt)) {
      hit.inx <- colnames(mydata) %in% colnames(mSetObj$dataSet$filt)
      mydata <- mydata[, hit.inx, drop = FALSE]
    }
    prenorm <- mydata
    mSetObj$dataSet$prenorm.cls <- mSetObj$dataSet$edit.cls
    if (substring(mSetObj$dataSet$format, 4, 5) == "ts") {
      mSetObj$dataSet$prenorm.facA <- mSetObj$dataSet$edit.facA
      mSetObj$dataSet$prenorm.facB <- mSetObj$dataSet$edit.facB
    }
  } else if (!is.null(mSetObj$dataSet$filt)) {
    prenorm <- mSetObj$dataSet$filt
    mSetObj$dataSet$prenorm.cls <- mSetObj$dataSet$filt.cls
    if (substring(mSetObj$dataSet$format, 4, 5) == "ts") {
      mSetObj$dataSet$prenorm.facA <- mSetObj$dataSet$filt.facA
      mSetObj$dataSet$prenorm.facB <- mSetObj$dataSet$filt.facB
    }
  } else {
    prenorm <- mSetObj$dataSet$data_proc
    mSetObj$dataSet$prenorm.cls <- mSetObj$dataSet$proc.cls
    if (substring(mSetObj$dataSet$format, 4, 5) == "ts") {
      mSetObj$dataSet$prenorm.facA <- mSetObj$dataSet$proc.facA
      mSetObj$dataSet$prenorm.facB <- mSetObj$dataSet$proc.facB
    }
  }
  mSetObj$dataSet$prenorm <- prenorm
  qs::qsave(prenorm, "prenorm.qs")
  mSetObj$dataSet$prenorm.smpl.nms <- rownames(prenorm)
  mSetObj$dataSet$prenorm.feat.nms <- colnames(prenorm)
  return(mSetObj)
}

$met.read_data
function (data,
          data.type = "conc", anal.type = "stat", paired = FALSE, # Parameters used to initialize dataSet object
          csvsep = ";", # optional: delimiter if reading CSV file
          format = "rowu", lbl.type = "disc",
          filt.metab = c(""),
          filt.smpl = c(""),
          filt.grp = c(""),
          filt.method = NULL,
          imp.method = "min",
          qcFilter = "F",
          qc.rsd = 25,
          all.rsd = NULL,
          export = TRUE
)
{
  mSetObj <- suppressWarnings(met.initialize(data.type = data.type, anal.type = anal.type, paired = paired))
  mSetObj$dataSet$cls.type <- lbl.type
  mSetObj$dataSet$format <- format
  if (!is.character(data)) {
    dat <- data
  } else {
    # Read table file
    if (file.exists(data)) {
      if (str_replace_all(data, ".{1,}\\.", "") == "csv") {
        dat <-
          read.csv(
            data,
            sep = csvsep,
            header = T,
            stringsAsFactors = F,
            fill = T,
            na.strings = "",
            quote = "",
            comment.char = "",
            check.names = F
          )
      } else if (str_replace_all(data, ".{1,}\\.", "") == "xls" |
                 str_replace(data, ".{1,}\\.", "") == "xlsx") {
        dat <- read_excel(data)
      } else if (str_replace_all(data, ".{1,}\\.", "") == "tsv") {
        dat <-
          read.csv(
            data,
            sep = "\t",
            header = T,
            stringsAsFactors = F,
            fill = T,
            na.strings = "",
            quote = "",
            comment.char = "",
            check.names = F
          )
      } else if (str_replace_all(data, ".{1,}\\.", "") == "txt") {
        dat <-
          read.table(
            data,
            sep = "\t",
            header = T,
            stringsAsFactors = F,
            fill = T,
            na.strings = "",
            quote = "",
            comment.char = "",
            check.names = F
          )
      } else {
        stop(
          "No compatible file format provided.
             Supported formats are: \\.txt (tab delimited), \\.csv (delimiters can be specified with the argument \"csvsep = \", \\.tsv, \\.xls, and \\.xlsx"
        )
      }
    } else {
      stop(paste0("File \"", data, "\" does not exist."), call. = F)
    }
  }
  if (class(dat) == "try-error" || ncol(dat) == 1) {
    AddErrMsg("Data format error. Failed to read in the data!")
    AddErrMsg("Make sure the data table is saved as comma separated values (.csv), XLS/XLSX, or tab separated values (.txt) format!")
    AddErrMsg("Please also check the followings: ")
    AddErrMsg("Either sample or feature names must in UTF-8 encoding; Latin, Greek letters are not allowed.")
    AddErrMsg("We recommend to use a combination of English letters, underscore, and numbers for naming purpose.")
    AddErrMsg("Make sure sample names and feature (peak, compound) names are unique.")
    AddErrMsg("Missing values should be blank or NA without quote.")
    AddErrMsg("Make sure the file delimeters are indicated correctly (e.g., csvsep = \",\")")
    return(0)
  }
  msg <- NULL

  if (substring(format, 1, 3) == "row") {
    msg <- c(msg, "Samples are in rows and features in columns")
    smpl.nms <- dat[, 1]
    dat[, 1] <- NULL
    if (lbl.type == "qc") {
      rownames(dat) <- smpl.nms
      qs::qsave(dat, file = "data_orig.qs")
      mSetObj$dataSet$data_orig <- dat
      mSetObj$dataSet$cmpd <- colnames(dat)
      return(1)
    }
    cls.lbl <- dat[, 1]
    conc <- dat[, -1, drop = FALSE]
    var.nms <- colnames(conc)
    if (lbl.type == "no") {
      cls.lbl <- rep(1, nrow(dat))
      conc <- dat[, , drop = FALSE]
      var.nms <- colnames(conc)
    }
  } else {
    msg <- c(msg, "Samples are in columns and features in rows.")
    if (lbl.type == "no") {
      cls.lbl <- rep(1, ncol(dat))
      conc <- t(dat[, -1])
      var.nms <- dat[, 1]
      smpl.nms <- colnames(dat[, -1])
    } else {
      var.nms <- dat[-1, 1]
      dat[, 1] <- NULL
      smpl.nms <- colnames(dat)
      cls.lbl <- dat[1, ]
      conc <- t(dat[-1, ])
    }
  }
  conc <- conc[,order(colnames(conc))]
  colnames(conc) <- make.names(colnames(conc))
  var.nms <- make.names(var.nms[order(var.nms)])
  mSetObj$dataSet$type.cls.lbl <- class(cls.lbl)
  empty.inx <- is.na(smpl.nms) | smpl.nms == ""
  if (sum(empty.inx) > 0) {
    msg <- c(msg, paste("<font color=\"red\">", sum(empty.inx),
                        "empty rows</font> were detected and excluded from your data."))
    smpl.nms <- smpl.nms[!empty.inx]
    cls.lbl <- cls.lbl[!empty.inx]
    conc <- conc[!empty.inx, ]
  }
  empty.inx <- is.na(cls.lbl) | cls.lbl == ""
  if (sum(empty.inx) > 0) {
    if (mSetObj$analSet$type != "roc") {
      msg <- c(msg, paste("<font color=\"red\">",
                          sum(empty.inx), "empty labels</font> were detected and excluded from your data."))
      smpl.nms <- smpl.nms[!empty.inx]
      cls.lbl <- cls.lbl[!empty.inx]
      conc <- conc[!empty.inx, ]
    }
    else {
      cls.lbl[is.na(cls.lbl)] <- ""
      msg <- c(msg, paste("<font color=\"orange\">",
                          sum(empty.inx), "new samples</font> were detected from your data."))
    }
  }
  if (anal.type == "roc") {
    if (length(unique(cls.lbl[!empty.inx])) > 2) {
      AddErrMsg("ROC analysis is only defined for two-group comparisions!")
      return(0)
    }
  }
  if (length(unique(smpl.nms)) != length(smpl.nms)) {
    dup.nm <- paste(smpl.nms[duplicated(smpl.nms)], collapse = " ")
    AddErrMsg("Duplicate sample names are not allowed!")
    AddErrMsg(dup.nm)
    return(0)
  }
  empty.inx <- is.na(var.nms) | var.nms == ""
  if (sum(empty.inx) > 0) {
    msg <- c(msg, paste("<font color=\"red\">", sum(empty.inx),
                        "empty features</font> were detected and excluded from your data."))
    var.nms <- var.nms[!empty.inx]
    conc <- conc[, !empty.inx]
  }
  if (length(unique(var.nms)) != length(var.nms)) {
    dup.nm <- paste(var.nms[duplicated(var.nms)], collapse = " ")
    AddErrMsg("Duplicate feature names are not allowed!")
    AddErrMsg(dup.nm)
    return(0)
  }
  if (anal.type == "mummichog") {
    is.rt <- mSetObj$paramSet$mumRT
    if (!is.rt) {
      mzs <- as.numeric(var.nms)
      if (sum(is.na(mzs) > 0)) {
        AddErrMsg("Make sure that feature names are numeric values (mass or m/z)!")
        return(0)
      }
    }
  }
  if (sum(is.na(iconv(smpl.nms))) > 0) {
    na.inx <- is.na(iconv(smpl.nms))
    nms <- paste(smpl.nms[na.inx], collapse = "; ")
    AddErrMsg(paste("No special letters (i.e. Latin, Greek) are allowed in sample names!",
                    nms, collapse = " "))
    return(0)
  }
  if (sum(is.na(iconv(var.nms))) > 0) {
    na.inx <- is.na(iconv(var.nms))
    nms <- paste(var.nms[na.inx], collapse = "; ")
    AddErrMsg(paste("No special letters (i.e. Latin, Greek) are allowed in feature names!",
                    nms, collapse = " "))
    return(0)
  }
  smpl.nms <- gsub("\\\\", "-", smpl.nms)
  url.smp.nms <- MetaboAnalystR:::CleanNames(smpl.nms)
  names(url.smp.nms) <- smpl.nms
  var.nms <- gsub("\\\\", "-", var.nms)
  url.var.nms <- MetaboAnalystR:::CleanNames(var.nms)
  names(url.var.nms) <- var.nms
  cls.lbl <- MetaboAnalystR:::ClearStrings(as.vector(cls.lbl))
  rownames(conc) <- smpl.nms
  colnames(conc) <- var.nms
  if (mSetObj$dataSet$paired) {
    mSetObj$dataSet$orig.cls <- mSetObj$dataSet$pairs <- cls.lbl
  } else {
    if (lbl.type == "disc") {
      mSetObj$dataSet$orig.cls <- mSetObj$dataSet$cls <- as.factor(as.character(cls.lbl))
      if (substring(format, 4, 5) == "ts") {
        mSetObj$dataSet$facA.type <- is.numeric(facA)
        mSetObj$dataSet$orig.facA <- mSetObj$dataSet$facA <- as.factor(as.character(facA))
        mSetObj$dataSet$facA.lbl <- facA.lbl
        mSetObj$dataSet$facB.type <- is.numeric(facB)
        mSetObj$dataSet$orig.facB <- mSetObj$dataSet$facB <- as.factor(as.character(facB))
        mSetObj$dataSet$facB.lbl <- facB.lbl
      }
    }
    else {
      mSetObj$dataSet$orig.cls <- mSetObj$dataSet$cls <- tryCatch({
        as.numeric(cls.lbl)
      }, warning = function(na) {
        print("Class labels must be numeric and continuous!")
        return(0)
      })
      if (mSetObj$dataSet$cls == 0) {
        AddErrMsg("Class labels must be numeric and continuous!")
        return(0)
      }
    }
  }
  if (mSetObj$dataSet$type == "conc") {
    mSetObj$dataSet$cmpd <- var.nms
  }
  mSetObj$dataSet$mumType <- "table"
  mSetObj$dataSet$url.var.nms <- url.var.nms
  mSetObj$dataSet$url.smp.nms <- url.smp.nms
  mSetObj$dataSet$data_orig <- conc
  qs::qsave(conc, file = "data_orig.qs")
  mSetObj$msgSet$read.msg <- c(msg, paste("The uploaded data file contains a ",
                                          nrow(conc), " (samples) by ", ncol(conc), " (",
                                          tolower(GetVariableLabel(mSetObj$dataSet$type)), ") data matrix.",
                                          sep = ""))

  mSetObj <- suppressWarnings(met.SanityCheck(mSetObj))
  mSetObj <- met.impute(mSetObj, method = imp.method)
  mSetObj <- met.PreparePrenormData(mSetObj)

  if(!filt.meta b== "" || !filt.smpl == "" || !filt.grp == "") {
    mSetObj <- GetGroupNames(mSetObj, "")
    mSetObj <- met.UpdateData(mSetObj, filt.metab = filt.metab,
                              filt.smpl = filt.smpl,
                              filt.grp = filt.grp)
    mSetObj<- met.PreparePrenormData(mSetObj)
  }
  if(!is.null(filt.method)){
    mSetObj <- met.FilterVariable(mSetObj, filt.method, qcFilter, rsd = qc.rsd, all.rsd=all.rsd)
    mSetObj<- met.PreparePrenormData(mSetObj)
  }
  dir.create(paste0(getwd(), "/met.ProcessedData"), showWarnings = F)
  mSetObj <- met.plot_missval(mSetObj, plot=F, export=export)
  mSetObj <- met.plot_detect(mSetObj, plot=F,  export=export)

  return(mSetObj)
}

$met.report
function(mSetObj, report.nm = NULL, ...){
  args <- list(...)
  for(i in 1:length(args)){
    assign(names(args)[i], args[[i]])
  }
  message("Render reports...")
  if(!is.null(report.nm)){
    wd <- paste0(getwd(), "/", report.nm)
  } else {
    wd <- paste(getwd(), "/Report_", format(Sys.time(),
                                            "%Y%m%d_%H%M%S"), sep = "")
  }
  dir.create(wd, showWarnings = F)
  file <- paste("C:/Users/nicwir/Documents/DTU_Biosustain/Scripts_and_Modelling/fluctuator/220111/R_package", "/Report_Met.Rmd",
                sep = "")
  rmarkdown::render(file, output_format = "all", output_dir = wd,
                    quiet = TRUE)
  message("Save RData object...")
  save(mSetObj, file = paste(wd, "results.RData", sep = "/"))
  message(paste0("Files saved in: '", wd, "'"))
  unlink("C:/Users/nicwir/Documents/DTU_Biosustain/Scripts_and_Modelling/fluctuator/220111/R_package/Plots", recursive = TRUE)
}

  $met.report_test_normalization
function(mSet_list, report.nm = NULL, ...){
  assertthat::assert_that(is.list(mSet_list))
  args <- list(...)
  for(i in 1:length(args)){
    assign(names(args)[i], args[[i]])
  }
  message("Render reports...")
  if(!is.null(report.nm)){
    wd <- paste0(getwd(), "/", report.nm)
  } else {
    wd <- paste(getwd(), "/met.Test_Normalization/Report_", format(Sys.time(),
                                                                   "%Y%m%d_%H%M%S"), sep = "")
  }
  dir.create(wd, showWarnings = F)
  file <- paste("C:/Users/nicwir/Documents/DTU_Biosustain/Scripts_and_Modelling/fluctuator/220111/R_package", "/Report_Met_TestNorm.Rmd",
                sep = "")
  rmarkdown::render(file, output_format = "all", output_dir = wd,
                    quiet = TRUE)
  message("Save RData object...")
  save(mSet_list, file = paste(wd, "results_test_normalizaton.RData", sep = "/"))
  message(paste0("Files saved in: '", wd, "'"))
  unlink(paste0(str_replace_all(tempdir(), "\\\\", "/"), "/Plots"), recursive = TRUE)
}

$met.SanityCheck
function (mSetObj = NA) {
  if(!is.null(mSetObj$dataSet$data_orig)){
    orig.data <- mSetObj$dataSet$data_orig
  } else if (file.exists("data_orig.qs")) {
    orig.data <- qs::qread("data_orig.qs")
  }
  else {
    return(0)
  }
  msg <- NULL
  cls <- mSetObj$dataSet$orig.cls
  mSetObj$dataSet$small.smpl.size <- 0
  if (mSetObj$dataSet$cls.type == "disc") {
    if (substring(mSetObj$dataSet$format, 4, 5) == "ts") {
      metadata <- mSetObj$dataSet$meta.info
      if (mSetObj$dataSet$design.type == "time") {
        msg <- c(msg, "The data is time-series data.")
      }
      else if (mSetObj$dataSet$design.type == "time0") {
        msg <- c(msg, "The data is time-series only data.")
      }
      else {
        msg <- c(msg, "The data is not time-series data.")
      }
      clsA.num <- length(levels(metadata[, 1]))
      clsB.num <- length(levels(metadata[, 2]))
      if (ncol(metadata) == 2) {
        msg <- c(msg, paste(clsA.num, "groups were detected in samples for factor",
                            colnames(metadata)[1]))
        msg <- c(msg, paste(clsB.num, "groups were detected in samples for factor",
                            colnames(metadata)[2]))
      }
      else {
        msg <- c(msg, paste0(clsA.num, " groups were detected from primary meta-data factor: ",
                             colnames(mSetObj$dataSet$meta.info)[1], "."))
      }
    }
    else {
      if (mSetObj$dataSet$paired) {
        msg <- c(msg, "Samples are paired.")
        if (!(mSetObj$dataSet$type == "conc" |
              mSetObj$dataSet$type == "specbin" | mSetObj$dataSet$type ==
              "pktable")) {
          pairs <- ReadPairFile()
          if (length(pairs) != length(mSetObj$dataSet$url.smp.nms)) {
            AddErrMsg("Error: the total paired names are not equal to sample names.")
            return(0)
          }
          else {
            inx <- match(rownames(orig.data), names(pairs))
            if (sum(is.na(inx)) > 0) {
              AddErrMsg("Error: some paired names not match the sample names.")
              return(0)
            }
            else {
              mSetObj$dataSet$pairs <- pairs[inx]
            }
          }
        }
        pairs <- mSetObj$dataSet$pairs
        qc.hits <- tolower(as.character(cls)) %in% "qc"
        if (sum(qc.hits) > 0) {
          AddErrMsg("<font color='red'>Error: QC samples not supported in paired analysis mode.</font>")
          AddErrMsg("You can perform QC filtering using regular two-group labels.")
          AddErrMsg("Then re-upload your data (without QC samples) for paired analysis.")
          return(0)
        }
        else {
          pairs <- as.numeric(pairs)
        }
        label <- as.numeric(pairs)
        cls <- as.factor(ifelse(label > 0, 1, 0))
        mSetObj$dataSet$pairs <- label
        lev <- unique(pairs)
        uni.cl <- length(lev)
        uni.cl.abs <- uni.cl/2
        sorted.pairs <- sort(pairs, index = TRUE)
        if (!all(sorted.pairs$x == c(-uni.cl.abs:-1,
                                     1:uni.cl.abs))) {
          AddErrMsg("There are some problems in paired sample labels! ")
          if (uni.cl.abs != round(uni.cl.abs)) {
            duplicates <- pairs[duplicated(pairs)]
            dup.msg <- paste0("Duplicated labels:",
                              duplicates)
            AddErrMsg(paste("The total samples must be of even number!",
                            dup.msg))
          }
          else {
            AddErrMsg(paste("And class labels between ",
                            -uni.cl.abs, " and 1, and between 1 and ",
                            uni.cl.abs, ".", sep = ""))
          }
          return(0)
        }
        else {
          msg <- c(msg, "The labels of paired samples passed sanity check.")
          msg <- c(msg, paste("A total of", uni.cl.abs,
                              "pairs were detected."))
          x <- sorted.pairs$ix[(uni.cl.abs + 1):uni.cl]
          y <- sorted.pairs$ix[uni.cl.abs:1]
          index <- as.vector(cbind(x, y))
          cls <- cls[index]
          pairs <- pairs[index]
          mSetObj$dataSet$pairs <- pairs
          mSetObj$dataSet$orig.cls <- cls
          orig.data <- orig.data[index, ]
          qs::qsave(orig.data, file = "data_orig.qs")
          mSetObj$dataSet$data_orig <- orig.data
        }
      }
      else {
        cls.lbl <- mSetObj$dataSet$orig.cls
        qb.inx <- tolower(cls.lbl) %in% c("qc",
                                          "blank")
        if (sum(qb.inx) > 0) {
          cls.Clean <- as.factor(as.character(cls.lbl[!qb.inx]))
        }
        else {
          cls.Clean <- cls.lbl
        }
        if (anal.type != "network") {
          if (min(table(cls.Clean)) < 3 | length(levels(cls.Clean)) <
              2) {
            AddErrMsg(paste("A total of", length(levels(cls.Clean)),
                            "groups found with", length(cls.Clean),
                            "samples."))
            AddErrMsg("<font color='red'>At least <b>two</b> groups and <b>three replicates</b> per group are required for analysis</font>!")
            AddErrMsg("You can click the <b>Edit Groups</b> button below to see the group labels for each sample and make corrections.")
            return(-1)
          }
        }
        msg <- c(msg, "Samples are not paired.")
      }
      cls.lbl <- mSetObj$dataSet$orig.cls
      qb.inx <- tolower(cls.lbl) %in% c("qc", "blank")
      if (sum(qb.inx) > 0) {
        cls.lbl <- as.factor(as.character(cls.lbl[!qb.inx]))
      }
      min.grp.size <- min(table(cls.lbl))
      cls.num <- length(levels(cls.lbl))
      if (cls.num/min.grp.size > 3) {
        mSetObj$dataSet$small.smpl.size <- 1
        msg <- c(msg, "<font color='red'>Too many groups with very small number of replicates!</font>")
        msg <- c(msg, "<font color='red'>Only a subset of methods will be available for analysis!</font>")
      }
      if (mSetObj$analSet$type == "ts") {
        msg <- c(msg, paste0(cls.num, " groups were detected from primary meta-data factor: ",
                             colnames(mSetObj$dataSet$meta.info)[1], "."))
        cls.vec <- vector()
        meta.info <- mSetObj$dataSet$meta.info
        meta.types <- mSetObj$dataSet$meta.types
        for (i in 1:ncol(meta.info)) {
          if (meta.types[i] == "disc") {
            cls.lbl <- meta.info[, i]
            qb.inx <- tolower(cls.lbl) %in% c("qc",
                                              "blank")
            if (sum(qb.inx) > 0) {
              cls.Clean <- as.factor(as.character(cls.lbl[!qb.inx]))
            }
            else {
              cls.Clean <- cls.lbl
            }
            meta.name <- colnames(meta.info)[i]
            if (min(table(cls.Clean)) < 3 | length(levels(cls.Clean)) <
                2) {
              cls.vec <- c(cls.vec, meta.name)
            }
          }
        }
      }
      else {
        msg <- c(msg, paste0(cls.num, " groups were detected in samples: ",
                             paste(levels(mSetObj$dataSet$cls), collapse = ", ")))
      }
      mSetObj$dataSet$cls.num <- cls.num
      mSetObj$dataSet$min.grp.size <- min.grp.size
    }
    if (substring(mSetObj$dataSet$format, 4, 5) == "ts") {
      metadata <- mSetObj$dataSet$meta.info
      nfacA <- metadata[, 1]
      nfacB <- metadata[, 2]
      if (mSetObj$dataSet$design.type == "time" |
          mSetObj$dataSet$design.type == "time0") {
        if (tolower(colnames(metadata)[1]) == "time") {
          time.fac <- nfacA
          exp.fac <- nfacB
        }
        else {
          time.fac <- nfacB
          exp.fac <- nfacA
        }
        ord.inx <- order(exp.fac)
      }
      else {
        ord.inx <- order(nfacA)
      }
      mSetObj$dataSet$orig.cls <- mSetObj$dataSet$orig.cls[ord.inx]
      mSetObj$dataSet$url.smp.nms <- mSetObj$dataSet$url.smp.nms[ord.inx]
      mSetObj$dataSet$facA <- mSetObj$dataSet$orig.facA <- metadata[,
                                                                    1][ord.inx]
      mSetObj$dataSet$facB <- mSetObj$dataSet$orig.facB <- metadata[,
                                                                    2][ord.inx]
      orig.data <- orig.data[ord.inx, ]
      mSetObj$dataSet$meta.info <- mSetObj$dataSet$meta.info[rownames(orig.data),
      ]
      qs::qsave(orig.data, file = "data_orig.qs")
      mSetObj$dataSet$data_orig <- orig.data
    }
    else {
      ord.inx <- order(mSetObj$dataSet$orig.cls)
      mSetObj$dataSet$orig.cls <- cls[ord.inx]
      mSetObj$dataSet$url.smp.nms <- mSetObj$dataSet$url.smp.nms[ord.inx]
      orig.data <- orig.data[ord.inx, , drop = FALSE]
      mSetObj$dataSet$data_orig <- orig.data
      qs::qsave(orig.data, file = "data_orig.qs")
      if (mSetObj$dataSet$paired) {
        mSetObj$dataSet$pairs <- mSetObj$dataSet$pairs[ord.inx]
      }
    }
  }
  msg <- c(msg, "Only English letters, numbers, underscore, hyphen and forward slash (/) are allowed.")
  msg <- c(msg, "Other special characters or punctuations (if any) will be stripped off.")
  int.mat <- orig.data
  if (ncol(int.mat) == 1) {
    if (anal.type == "roc") {
      mSetObj$dataSet$roc_cols <- 1
    }
    else {
      AddErrMsg("One-column data is only supported for biomarker analysis.")
      return(0)
    }
  }
  else {
    mSetObj$dataSet$roc_cols <- 2
  }
  rowNms <- rownames(int.mat)
  colNms <- colnames(int.mat)
  naNms <- sum(is.na(int.mat))
  for (c in 1:ncol(int.mat)) {
    if (class(int.mat[, c]) == "integer64") {
      int.mat[, c] <- as.double(int.mat[, c])
    }
  }
  num.mat <- apply(int.mat, 2, as.numeric)
  if (sum(is.na(num.mat)) > naNms) {
    num.mat <- apply(int.mat, 2, function(x) as.numeric(gsub(",",
                                                             "", x)))
    if (sum(is.na(num.mat)) > naNms) {
      msg <- c(msg, "Non-numeric values were found and replaced by NA.")
    }
    else {
      msg <- c(msg, "All data values are numeric.")
    }
  }
  else {
    msg <- c(msg, "All data values are numeric.")
  }
  int.mat <- num.mat
  rownames(int.mat) <- rowNms
  colnames(int.mat) <- colNms
  varCol <- apply(int.mat, 2, var, na.rm = T)
  constCol <- (varCol == 0 | is.na(varCol))
  constNum <- sum(constCol, na.rm = T)
  if (constNum > 0) {
    msg <- c(msg, paste("<font color=\"red\">", constNum,
                        "features with a constant or single value across samples were found and deleted.</font>"))
    int.mat <- int.mat[, !constCol, drop = FALSE]
  }
  totalCount <- nrow(int.mat) * ncol(int.mat)
  naCount <- sum(is.na(int.mat)|int.mat==0)
  naPercent <- round(100 * naCount/totalCount, 1)
  mSetObj$dataSet$missingCount <- naCount
  msg <- c(msg, paste("A total of ", naCount, " (",
                      naPercent, "%) missing or zero values were detected.",
                      sep = ""))
  mSetObj$msgSet$missing.msg <- msg[length(msg)]

  qs::qsave(as.data.frame(int.mat), "preproc.qs")
  mSetObj$dataSet$preproc <- as.data.frame(int.mat)
  mSetObj$dataSet$proc.cls <- mSetObj$dataSet$cls <- mSetObj$dataSet$orig.cls
  if (is.null(mSetObj$dataSet$meta.info)) {
    mSetObj$dataSet$meta.info <- data.frame(mSetObj$dataSet$cls)
    colnames(mSetObj$dataSet$meta.info) = "Class"
  }
  if (substring(mSetObj$dataSet$format, 4, 5) == "ts") {
    mSetObj$dataSet$proc.facA <- mSetObj$dataSet$orig.facA
    mSetObj$dataSet$proc.facB <- mSetObj$dataSet$orig.facB
  }
  mSetObj$msgSet$check.msg <- c(mSetObj$msgSet$read.msg, msg)

  print(mSetObj$msgSet$check.msg)

  return(mSetObj)
}

$met.test_normalization
function(mSetObj,
         test_conditions = NULL,
         ref = NULL, # Input the name of the reference sample or the reference feature, use " " around the name (for rowNorm="CompNorm" or "SpecNorm").
         class_order = FALSE,
         alpha = 0.05,
         lfc = 2,
         posthoc_method = "tukey",
         nonpar = FALSE,
         export = FALSE, # Shall plots be exported as PDF and PNG files?
         export.format = "pdf",
         dpi = 300, # resolution of PNG images
         pls.data = "all",
         permut.num = 500,
         volcano=TRUE,
         vip.thresh = 1,
         volcano.contrast=NULL,
         volcano.add_names = FALSE, #show labels next to circles in volcano plots
         report = FALSE, # Shall a report (HTML and PDF) be created?
         report.nm = NULL, # Folder name for created report (if report = TRUE)
         export.dir = "met.Test_Normalization"){

  dir.create(paste0(getwd(), "/", export.dir), showWarnings = F)
  wd <- getwd()
  setwd(paste0(wd, "/", export.dir))

  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
  }
  conditions <- levels(mSetObj$dataSet$cls)
  if(!is.null(mSetObj$dataSet$prenorm)){
    proc.data <- data.frame(mSetObj$dataSet$prenorm)
  } else if(!is.null(mSetObj$dataSet$data_proc)){
    proc.data <- mSetObj$dataSet$data_proc
  } else {
    proc.data <- qs::qread("data_proc.qs")
  }
  cntrst <- apply(utils::combn(conditions, 2), 2, paste,
                  collapse = " - ")
  test.cond_tmp <- list()
  test.cond_list <- list()
  if(is.null(test_conditions)){
    test_conditions <- c("NULL", "MeanCenter", "LogNorm", "LogNorm/AutoNorm", "CrNorm", "CrNorm/AutoNorm", "AutoNorm", "RangeNorm", "ParetoNorm")
  }
  for(i in 1:length(test_conditions)){
    test.cond_tmp[i] <- str_split(test_conditions[i], "/")
    test.cond_list[i] <- test.cond_tmp[i]
    for(j in 1:length(test.cond_tmp[[i]])){
      if(!any(is.element(c("NULL", "MeanCenter", "AutoNorm", "ParetoNorm", "RangeNorm",
                           "LogNorm", "CrNorm",
                           "QuantileNorm", "CompNorm", "SumNorm", "MedianNorm", "SpecNorm"), test.cond_tmp[[i]][j]))){
        stop(paste0("'", test.cond_tmp[[i]][j], "' is not a valid option for 'rowNorm', 'transNorm', or 'scaleNorm'. Valid options are:\n",
                    "rowNorm: 'NULL', 'QuantileNorm', 'CompNorm', 'SumNorm', 'MedianNorm', 'SpecNorm'\n",
                    "transNorm: 'NULL', 'LogNorm', 'CrNorm'\n",
                    "scaleNorm: 'NULL', 'MeanCenter', 'AutoNorm', 'ParetoNorm', 'RangeNorm'"))
      }
    }
    if(length(test.cond_tmp[[i]])==1 & any(is.element(c("LogNorm", "CrNorm"), test.cond_tmp[[i]][1]))){
      test.cond_list[[i]][2] <- test.cond_tmp[[i]][1]
      test.cond_list[[i]][1] <- "NULL"
      test.cond_list[[i]][3] <- "NULL"
    }
    if(length(test.cond_tmp[[i]])==1 & any(is.element(c("NULL", "MeanCenter", "AutoNorm", "ParetoNorm", "RangeNorm"), test.cond_tmp[[i]][1]))){
      test.cond_list[[i]][3] <- test.cond_tmp[[i]][1]
      test.cond_list[[i]][1] <- "NULL"
      test.cond_list[[i]][2] <- "NULL"
    }
    if(length(test.cond_tmp[[i]])==2 & any(is.element(c("LogNorm", "CrNorm"), test.cond_tmp[[i]][1]))){
      test.cond_list[[i]][2] <- test.cond_tmp[[i]][1]
      test.cond_list[[i]][1] <- "NULL"
      test.cond_list[[i]][3] <- test.cond_tmp[[i]][2]
    }
    if(length(test.cond_tmp[[i]])==2 & any(is.element(c("NULL", "MeanCenter", "AutoNorm", "ParetoNorm", "RangeNorm"), test.cond_tmp[[i]][1]))){
      test.cond_list[[i]][3] <- test.cond_tmp[[i]][1]
      test.cond_list[[i]][1] <- "NULL"
      test.cond_list[[i]][2] <- test.cond_tmp[[i]][2]
    }
    if(length(test.cond_tmp[[i]])==2 & (is.element("NULL", test.cond_tmp[[i]][1]) | is.element("NULL", test.cond_tmp[[i]][2]))){
      if(any(is.element(c("LogNorm", "CrNorm"), test.cond_list[[i]][1]))){
        test.cond_list[[i]][1] <- "NULL"
        test.cond_list[[i]][2] <- test.cond_tmp[[i]][1]
        test.cond_list[[i]][3] <- "NULL"
      } else if(any(is.element(c("NULL", "MeanCenter", "AutoNorm", "ParetoNorm", "RangeNorm"), test.cond_tmp[[i]][1]))){
        test.cond_list[[i]][1] <- "NULL"
        test.cond_list[[i]][2] <- "NULL"
        test.cond_list[[i]][3] <- test.cond_tmp[[i]][1]
      }
    }
    names(test.cond_list)[i] <- paste0(if(test.cond_list[[i]][1]=="NULL" &
                                          test.cond_list[[i]][2]=="NULL" &
                                          test.cond_list[[i]][3]=="NULL") {"untransformed"},
                                       if_else(test.cond_list[[i]][1]!="NULL",test.cond_list[[i]][1],""),
                                       if_else(test.cond_list[[i]][1]!="NULL", "_", ""),
                                       if_else(test.cond_list[[i]][2]!="NULL",test.cond_list[[i]][2],""),
                                       if_else(test.cond_list[[i]][2]!="NULL" & test.cond_list[[i]][3]!="NULL", "_", ""),
                                       if_else(test.cond_list[[i]][3]!="NULL",test.cond_list[[i]][3],"")
    )
  }
  list_names <- names(test.cond_list)
  mSet_list <- list()
  for(i in 1:length(test.cond_list)){
    mSet_list[[i]] <- met.normalize(mSetObj, rowNorm = test.cond_list[[i]][1],
                                    transNorm = test.cond_list[[i]][2],
                                    scaleNorm = test.cond_list[[i]][3], ref = ref)
  }
  names(mSet_list) <- list_names
  mSet_list <- sapply(1:length(mSet_list), function(x) met.ANOVA.Anal(mSet_list[[x]], nonpar = nonpar, thresh = alpha, post.hoc = posthoc_method, all_results = TRUE, silent = TRUE), simplify = FALSE,
                      USE.NAMES = TRUE )
  mSet_list <- sapply(1:length(mSet_list), function(x) PCA.Anal(mSet_list[[x]]), simplify = FALSE, USE.NAMES = TRUE )
  mSet_list <- sapply(1:length(mSet_list), function(x) met.PLSR.Anal(mSet_list[[x]], reg = class_order, data = pls.data), simplify = FALSE, USE.NAMES = TRUE )
  options(warn=1)
  for (i in 1:length(mSet_list)){
    cat(paste0("Performing cross-validation test for PLS-DA model with condition: ", names(test.cond_list)[i], "\n"))
    mSet_list[[i]] <-
      met.PLSDA.CV(
        mSet_list[[i]],
        methodName = "L",
        compNum = MetaboAnalystR:::GetDefaultPLSCVComp(mSet_list[[i]]),
        data = pls.data
      )

  }
  options(warn=0)
  names(mSet_list) <- list_names

  anova_significant <- list()
  for (i in 1:length(mSet_list)){
    mSet_list[[i]][["analSet"]][["plsda"]][["vip_sig_comp1"]] <- mSet_list[[i]][["analSet"]][["plsda"]][["vip.mat"]][mSet_list[[i]][["analSet"]][["plsda"]][["vip.mat"]][,1]>1, 1]
    mSet_list[[i]][["analSet"]][["plsda"]][["class_order"]] <- class_order
    mSet_list[[i]][["analSet"]][["plsda"]][["vip_sig_comp2"]] <- mSet_list[[i]][["analSet"]][["plsda"]][["vip.mat"]][mSet_list[[i]][["analSet"]][["plsda"]][["vip.mat"]][,2]>1, 2]
    mSet_list[[i]][["analSet"]][["aov"]][["significant"]] <- mSet_list[[i]][["analSet"]][["aov"]][["inx.imp"]][mSet_list[[i]][["analSet"]][["aov"]][["inx.imp"]]==TRUE]
    mSet_list[[i]][["analSet"]][["contrasts"]] <- cntrst

  }
  names(mSet_list) <- list_names

  mSet_list$vip_sig_comp1 <- lapply(1:length(list_names), function(x) names(mSet_list[[x]][["analSet"]][["plsda"]][["vip_sig_comp1"]]) )
  mSet_list$vip_sig_comp2 <- lapply(1:length(list_names), function(x) names(mSet_list[[x]][["analSet"]][["plsda"]][["vip_sig_comp2"]]) )
  mSet_list$anova_significant <- lapply(1:length(list_names), function(x) names(mSet_list[[x]][["analSet"]][["aov"]][["significant"]]) )
  names(mSet_list$vip_sig_comp1) <- names(mSet_list$vip_sig_comp2) <- names(mSet_list$anova_significant) <- list_names
  mSet_list$vip_sig_comp1 <- mSet_list$vip_sig_comp1[!sapply(mSet_list$vip_sig_comp1,is.null)]
  mSet_list$vip_sig_comp2 <- mSet_list$vip_sig_comp2[!sapply(mSet_list$vip_sig_comp2,is.null)]
  mSet_list$anova_significant <- mSet_list$anova_significant[!sapply(mSet_list$anova_significant,is.null)]

  f_venn <- function(list, title, export = T, format = "pdf", imgName = "Venn") {
    venn <- ggVennDiagram(
      list,
      color = "black",
      lwd = 0.8,
      lty = 1,
      label_alpha = 0,
      label_size = 4.5,
      set_size = 5.5) +
      scale_fill_gradient(low = "#fcfeff", high = "#4981BF") +
      scale_x_continuous(expand = expansion(mult = 0.3)) + ggtitle(title) +
      theme(plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5),
        legend.position = c(0.9, 0.5),
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))
    if (export == TRUE) {
      w = 10
      h = 8
      if (format == "pdf") {
        pdf(
          file = imgName,
          width = w,
          height = h,
          bg = "white",
          onefile = FALSE)
      }
      else {
        Cairo::Cairo(
          file = imgName,
          unit = "in",
          dpi = dpi,
          width = w,
          height = h,
          type = format,
          bg = "white")
      }
      print(venn)
      dev.off()
    }
    return(venn)
  }
  dir.create("Plots/Venn Diagrams", showWarnings = F)
  # Venn Diagram of features with VIP1 > 1
  if(length(mSet_list$vip_sig_comp1)<5){
    mSet_list$plot.venn_VIP1 <- f_venn(mSet_list$vip_sig_comp1, "PLS-DA Component 1 \u2013 VIP score > 1", export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_PLSDA_vip_comp1.", export.format))
  } else if(length(mSet_list$vip_sig_comp1)<7){
    mSet_list$plot.venn_VIP1 <- list()
    mSet_list$plot.venn_VIP1[[1]] <-   f_venn(mSet_list$vip_sig_comp1[1:4], paste0("PLS-DA Component 1 \u2013 VIP score > " ,vip.thresh), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_PLSDA_vip_comp1_A.", export.format))
    mSet_list$plot.venn_VIP1[[2]] <-   f_venn(mSet_list$vip_sig_comp1[c(1,5:length(mSet_list$vip_sig_comp1))], paste0("PLS-DA Component 1 \u2013 VIP score > " ,vip.thresh), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_PLSDA_vip_comp1_B.", export.format))
  } else if(length(mSet_list$vip_sig_comp1)<10){
    mSet_list$plot.venn_VIP1 <- list()
    mSet_list$plot.venn_VIP1[[1]] <- f_venn(mSet_list$vip_sig_comp1[1:5], paste0("PLS-DA Component 1 \u2013 VIP score > " ,vip.thresh), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_PLSDA_vip_comp1_A.", export.format))
    mSet_list$plot.venn_VIP1[[2]] <- f_venn(mSet_list$vip_sig_comp1[c(1,6:length(mSet_list$vip_sig_comp1))], paste0("PLS-DA Component 1 \u2013 VIP score > " ,vip.thresh), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_PLSDA_vip_comp1_B.", export.format))
  } else {
    mSet_list$plot.venn_VIP1 <- list()
    mSet_list$plot.venn_VIP1[[1]] <- f_venn(mSet_list$vip_sig_comp1[1:5], paste0("PLS-DA Component 1 \u2013 VIP score > " ,vip.thresh), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_PLSDA_vip_comp1_A.", export.format))
    mSet_list$plot.venn_VIP1[[2]] <- f_venn(mSet_list$vip_sig_comp1[c(1,6:10)], paste0("PLS-DA Component 1 \u2013 VIP score > " ,vip.thresh), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_PLSDA_vip_comp1_B.", export.format))
    mSet_list$plot.venn_VIP1[[3]] <- f_venn(mSet_list$vip_sig_comp1[c(1,11:length(mSet_list$vip_sig_comp1))], paste0("PLS-DA Component 1 \u2013 VIP score > " ,vip.thresh), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_PLSDA_vip_comp1_C.", export.format))
  }

  # Venn Diagram of features with VIP2 > 1
  if(length(mSet_list$vip_sig_comp2)<5){
    mSet_list$plot.venn_VIP2 <- f_venn(mSet_list$vip_sig_comp2, paste0("PLS-DA Component 2 \u2013 VIP score > " ,vip.thresh), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_PLSDA_vip_comp2.", export.format))

  } else if(length(mSet_list$vip_sig_comp2)<7){
    mSet_list$plot.venn_VIP2 <- list()
    mSet_list$plot.venn_VIP2[[1]] <- f_venn(mSet_list$vip_sig_comp2[1:4], paste0("PLS-DA Component 2 \u2013 VIP score > " ,vip.thresh), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_PLSDA_vip_comp2_A.", export.format))
    mSet_list$plot.venn_VIP2[[2]] <- f_venn(mSet_list$vip_sig_comp2[c(1,5:length(mSet_list$vip_sig_comp2))], paste0("PLS-DA Component 2 \u2013 VIP score > " ,vip.thresh), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_PLSDA_vip_comp2_B.", export.format))

  } else if(length(mSet_list$vip_sig_comp2)<10){
    mSet_list$plot.venn_VIP2 <- list()
    mSet_list$plot.venn_VIP2[[1]] <- f_venn(mSet_list$vip_sig_comp2[1:5], paste0("PLS-DA Component 2 \u2013 VIP score > " ,vip.thresh), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_PLSDA_vip_comp2_A.", export.format))
    mSet_list$plot.venn_VIP2[[2]] <- f_venn(mSet_list$vip_sig_comp2[c(1,6:length(mSet_list$vip_sig_comp2))], paste0("PLS-DA Component 2 \u2013 VIP score > " ,vip.thresh), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_PLSDA_vip_comp2_B.", export.format))
  } else {
    mSet_list$plot.venn_VIP2 <- list()
    mSet_list$plot.venn_VIP2[[1]] <- f_venn(mSet_list$vip_sig_comp2[1:5], paste0("PLS-DA Component 2 \u2013 VIP score > " ,vip.thresh), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_PLSDA_vip_comp2_A.", export.format))
    mSet_list$plot.venn_VIP2[[2]] <- f_venn(mSet_list$vip_sig_comp2[c(1,6:10)], paste0("PLS-DA Component 2 \u2013 VIP score > " ,vip.thresh), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_PLSDA_vip_comp2_B.", export.format))
    mSet_list$plot.venn_VIP2[[3]] <- f_venn(mSet_list$vip_sig_comp2[c(1,11:length(mSet_list$vip_sig_comp2))], paste0("PLS-DA Component 2 \u2013 VIP score > " ,vip.thresh), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_PLSDA_vip_comp2_C.", export.format))
  }

  # Venn Diagram of significant features based on ANOVA test
  if(length(mSet_list$anova_significant)<5){
    mSet_list$plot.venn_anova_signif <- f_venn(mSet_list$anova_significant, paste0("Features with FDR ", "\u2264 ", alpha), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_ANOVA_signif.", export.format))
  } else if(length(mSet_list$anova_significant)<7){
    mSet_list$plot.venn_anova_signif <- list()
    mSet_list$plot.venn_anova_signif[[1]] <- f_venn(mSet_list$anova_significant[1:4], paste0("Features with FDR ", "\u2264 ", alpha), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_ANOVA_signif_A.", export.format))
    mSet_list$plot.venn_anova_signif[[2]] <- f_venn(mSet_list$anova_significant[c(1,5:length(mSet_list$anova_significant))], paste0("Features with FDR ", "\u2264 ", alpha), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_ANOVA_signif_B.", export.format))

  } else if(length(mSet_list$anova_significant)<10){
    mSet_list$plot.venn_anova_signif <- list()
    mSet_list$plot.venn_anova_signif[[1]] <- f_venn(mSet_list$anova_significant[1:5], paste0("Features with FDR ", "\u2264 ", alpha), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_ANOVA_signif_A.", export.format))
    mSet_list$plot.venn_anova_signif[[2]] <- f_venn(mSet_list$anova_significant[c(1,6:length(mSet_list$anova_significant))], paste0("Features with FDR ", "\u2264 ", alpha), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_ANOVA_signif_B.", export.format))
  } else {
    mSet_list$plot.venn_anova_signif <- list()
    mSet_list$plot.venn_anova_signif[[1]] <- f_venn(mSet_list$anova_significant[1:5], paste0("Features with FDR ", "\u2264 ", alpha), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_ANOVA_signif_A.", export.format))
    mSet_list$plot.venn_anova_signif[[2]] <- f_venn(mSet_list$anova_significant[c(1,6:10)], paste0("Features with FDR ", "\u2264 ", alpha), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_ANOVA_signif_B.", export.format))
    mSet_list$plot.venn_anova_signif[[3]] <- f_venn(mSet_list$anova_significant[c(1,11:length(mSet_list$anova_significant))], paste0("Features with FDR ", "\u2264 ", alpha), export = export, format = export.format, imgName= paste0("Plots/Venn Diagrams/Venn_ANOVA_signif_C.", export.format))
  }

  # Create Venn Diagram of significant features and VIP>vip.thresh components
  anova_vip1_list <- list()
  names_common <- intersect(names(mSet_list$anova_significant), names(mSet_list$vip_sig_comp1))
  for (i in 1:length(names_common)){
    anova_vip1_list[[i]] <- list(mSet_list$anova_significant[[match(names_common[i], names(mSet_list$anova_significant))]], mSet_list$vip_sig_comp1[[match(names_common[i], names(mSet_list$vip_sig_comp1))]])
    names(anova_vip1_list[[i]]) <- c("ANOVA_signif", paste0("PLS-DA_VIP1 > ", vip.thresh))
  }
  names(anova_vip1_list) <- names_common

  mSet_list$plot_venn_anova_vs_vip1 <- list()
  only_vip1_list <- list()

  for(i in 1:length(names(anova_vip1_list))){
    if(!is.null(anova_vip1_list[[i]][[1]]) && !is.null(anova_vip1_list[[i]][[2]])){
      mSet_list$plot_venn_anova_vs_vip1[[i]] <- ggVennDiagram(anova_vip1_list[[i]], color = "black", lwd = 0.8, lty = 1,label_alpha = 0, label_size = 4.5, set_size = 5.5) +
        scale_fill_gradient(low = "#fcfeff", high = "#4981BF") + ggtitle(names(anova_vip1_list[i])) +
        theme(plot.title = element_text(hjust=0.5, vjust = 5, size = 20, face = "bold"),
              legend.key.size = unit(1, "cm"),
              legend.title = element_text(size=16),
              legend.text = element_text(size=14))
    } else {
      mSet_list$plot_venn_anova_vs_vip1[[i]] <- ggplot()+geom_blank()+theme_bw()
    }
    df_i = if (nrow(data.frame(
      "VIP1 - ANOVA" = Setdiff(anova_vip1_list[[i]][2], anova_vip1_list[[i]][1]),
      check.names = F
    )) != 0) {
      data.frame(
        "VIP1 - ANOVA" = Setdiff(anova_vip1_list[[i]][2], anova_vip1_list[[i]][1]),
        check.names = F
      )
    } else {
      df_i = ""
    }
    only_vip1_list[[i]] <- ggtexttable(df_i, rows = NULL, theme = ttheme(base_style = "light", base_size = 12)) %>%
      {if(!is.null(nrow(df_i))) tab_add_hline(tab = ., at.row = 1:2, row.side = "top", linewidth = 2) else .}

    mSet_list$plot_venn_anova_vs_vip1[[i]] <- ggarrange(mSet_list$plot_venn_anova_vs_vip1[[i]], only_vip1_list[[i]],
                                                        ncol = 2, nrow = 1,
                                                        widths = c(10,3))
  }
  names(mSet_list$plot_venn_anova_vs_vip1) <- names_common

  # Extract elements of Venn Diagram VIP1
  mSet_list$vip_sig_comp1[["combs"]] <-
    unlist(lapply(1:length(mSet_list$vip_sig_comp1), function(j)
      combn(names(mSet_list$vip_sig_comp1), j, simplify = FALSE)),
      recursive = FALSE)
  names(mSet_list$vip_sig_comp1[["combs"]]) <-
    sapply(mSet_list$vip_sig_comp1[["combs"]], function(i)
      paste0(i, collapse = " & "))
  mSet_list$vip_sig_comp1[["combs>0"]] <-
    lapply(mSet_list$vip_sig_comp1[["combs"]], function(i)
      Setdiff(mSet_list$vip_sig_comp1[i], mSet_list$vip_sig_comp1[setdiff(names(mSet_list$vip_sig_comp1), i)])) %>%
    Filter(function(x)
      any(length(x) != 0), .)

  # Extract elements of Venn Diagram VIP2
  mSet_list$vip_sig_comp2[["combs"]] <-
    unlist(lapply(1:length(mSet_list$vip_sig_comp2), function(j)
      combn(names(mSet_list$vip_sig_comp2), j, simplify = FALSE)),
      recursive = FALSE)
  names(mSet_list$vip_sig_comp2[["combs"]]) <-
    sapply(mSet_list$vip_sig_comp2[["combs"]], function(i)
      paste0(i, collapse = " & "))
  mSet_list$vip_sig_comp2[["combs>0"]] <-
    lapply(mSet_list$vip_sig_comp2[["combs"]], function(i)
      Setdiff(mSet_list$vip_sig_comp2[i], mSet_list$vip_sig_comp2[setdiff(names(mSet_list$vip_sig_comp2), i)])) %>%
    Filter(function(x)
      any(length(x) != 0), .)

  # Extract elements of Venn ANOVA significant
  mSet_list$anova_significant[["combs"]] <-
    unlist(lapply(1:length(mSet_list$anova_significant), function(j)
      combn(names(mSet_list$anova_significant), j, simplify = FALSE)),
      recursive = FALSE)
  names(mSet_list$anova_significant[["combs"]]) <-
    sapply(mSet_list$anova_significant[["combs"]], function(i)
      paste0(i, collapse = " & "))
  mSet_list$anova_significant[["combs>0"]] <-
    lapply(mSet_list$anova_significant[["combs"]], function(i)
      Setdiff(mSet_list$anova_significant[i], mSet_list$anova_significant[setdiff(names(mSet_list$anova_significant), i)])) %>%
    Filter(function(x)
      any(length(x) != 0), .)
  if (export == TRUE) {
    w = 10
    h = 8
    for(i in 1:length(mSet_list$plot_venn_anova_vs_vip1)){
      imgName = paste0("Plots/Venn Diagrams/Venn_ANOVA-signif_vs_PLSDA-VIP1_", names(mSet_list$plot_venn_anova_vs_vip1)[i], ".", export.format)
      if (export.format == "pdf") {
        pdf(
          file = imgName,
          width = w,
          height = h,
          bg = "white",
          onefile = FALSE
        )
      }
      else {
        Cairo::Cairo(
          file = imgName,
          unit = "in",
          dpi = dpi,
          width = w,
          height = h,
          type = export.format,
          bg = "white"
        )
        print(mSet_list$plot_venn_anova_vs_vip1[[i]])
        dev.off()
      }
    }
  }
  names(mSet_list) <- c(list_names, names(mSet_list[(length(list_names)+1):length(mSet_list)]))

  mSet_list <- c(lapply(1:length(list_names), function(x) met.plot_FeatureNormSummary(mSetObj = mSet_list[[x]], imgName = paste0("Plots/FeatureNormSummary_", names(mSet_list)[x]),
                                                                                      plot = F, show_prenorm = F, export = export, format = export.format)),
                 mSet_list[(length(list_names)+1):length(mSet_list)])

  names(mSet_list) <- c(list_names, names(mSet_list[(length(list_names)+1):length(mSet_list)]))

  mSet_list <- c(lapply(1:length(list_names), function(x) met.plot_SampleNormSummary(mSetObj = mSet_list[[x]], imgName = paste0("Plots/SampleNormSummary_", names(mSet_list)[x]),
                                                                                     plot = F, show_prenorm = F, export = export, format = export.format)),
                 mSet_list[(length(list_names)+1):length(mSet_list)])


  names(mSet_list) <- c(list_names, names(mSet_list[(length(list_names)+1):length(mSet_list)]))

  dir.create("Plots/ANOVA", showWarnings = F)

  mSet_list <- c(lapply(1:length(list_names), function(x) met.plot_ANOVA(mSetObj = mSet_list[[x]], imgName = paste0("Plots/ANOVA/ANOVA-Plot_", names(mSet_list)[x]),
                                                                         format = export.format, subtitle = TRUE,
                                                                         plot = F, export = export, dpi = if_else(export.format=="pdf", 72, 300))),
                 mSet_list[(length(list_names)+1):length(mSet_list)])
  names(mSet_list) <- c(list_names, names(mSet_list[(length(list_names)+1):length(mSet_list)]))

  dir.create("Plots/PCA", showWarnings = F)

  mSet_list <- c(lapply(1:length(list_names), function(x) met.plot_PCA2DScore(mSetObj = mSet_list[[x]], imgName = paste0("Plots/PCA/PCA-2DScore_", names(mSet_list)[x]),
                                                                              pcx = 1, pcy = 2, reg = 0.95, show = 0, format = export.format,
                                                                              plot = F, export = export, dpi = if_else(export.format=="pdf", 72, 300), subtitle = TRUE)),
                 mSet_list[(length(list_names)+1):length(mSet_list)])
  names(mSet_list) <- c(list_names, names(mSet_list[(length(list_names)+1):length(mSet_list)]))

  mSet_list <- c(lapply(1:length(list_names), function(x) met.plot_PCA2DLoading(mSetObj = mSet_list[[x]], imgName = paste0("Plots/PLSDA/PCA-Loadings_", names(mSet_list)[x]),
                                                                              inx1 = 1, inx2 = 2, plot = F, export = export, format = export.format, subtitle = TRUE)),
                 mSet_list[(length(list_names)+1):length(mSet_list)])

  names(mSet_list) <- c(list_names, names(mSet_list[(length(list_names)+1):length(mSet_list)]))

  dir.create("Plots/PLSDA", showWarnings = F)

  mSet_list <- c(lapply(1:length(list_names), function(x) met.plot_PLS2DScore(mSetObj = mSet_list[[x]], imgName = paste0("Plots/PLSDA/PLS-2DScore_", names(mSet_list)[x]),
                                                                              inx1 = 1, inx2 = 2, reg = 0.95, show = 0, format = export.format, subtitle = TRUE,
                                                                              plot = F, export = export, dpi = if_else(export.format=="pdf", 72, 300))),
                 mSet_list[(length(list_names)+1):length(mSet_list)])
  names(mSet_list) <- c(list_names, names(mSet_list[(length(list_names)+1):length(mSet_list)]))

  mSet_list <- c(lapply(1:length(list_names), function(x) met.plot_PLS2DLoading(mSetObj = mSet_list[[x]], imgName = paste0("Plots/PLSDA/PLS-Loadings_", names(mSet_list)[x]),
                                                                              inx1 = 1, inx2 = 2, plot = F, export = export, format = export.format)),
                 mSet_list[(length(list_names)+1):length(mSet_list)])

  names(mSet_list) <- c(list_names, names(mSet_list[(length(list_names)+1):length(mSet_list)]))

  mSet_list <- c(lapply(1:length(list_names), function(x) suppressWarnings(met.plot_PLS_Imp(mSetObj = mSet_list[[x]], imgName = paste0("Plots/PLSDA/PLSDA-coef_coef.mean_", names(mSet_list)[x]),
                                                                                            type = "coef", feat.nm = "coef.mean", feat.num = 15,
                                                                                            color.BW = FALSE, plot = F, export = export, format = export.format))),
                 mSet_list[(length(list_names)+1):length(mSet_list)])

  names(mSet_list) <- c(list_names, names(mSet_list[(length(list_names)+1):length(mSet_list)]))

  mSet_list <- c(lapply(1:length(list_names), function(x) suppressWarnings(met.plot_PLS_Imp(mSetObj = mSet_list[[x]], imgName = paste0("Plots/PLSDA/PLSDA-VIP1_", names(mSet_list)[x]),
                                                                                            type = "vip", feat.nm = "Comp. 1", feat.num = 15,
                                                                                            color.BW = FALSE, plot = F, export = export, format = export.format))),
                 mSet_list[(length(list_names)+1):length(mSet_list)])

  names(mSet_list) <- c(list_names, names(mSet_list[(length(list_names)+1):length(mSet_list)]))

  for(i in levels(mSet_list[[1]][["dataSet"]][["cls"]])){
    names(mSet_list) <- c(list_names, names(mSet_list[(length(list_names)+1):length(mSet_list)]))
    mSet_list <- c(lapply(1:length(list_names), function(x) suppressWarnings(met.plot_PLS_Imp(mSetObj = mSet_list[[x]], imgName = paste0("Plots/PLSDA/PLSDA-coef_", i, "_", names(mSet_list)[x]),
                                                                                              type = "coef", feat.nm = i, feat.num = 15,
                                                                                              color.BW = FALSE, plot = F, export = export, format = export.format))),
                   mSet_list[(length(list_names)+1):length(mSet_list)])
  }

  mSet_list <- c(lapply(1:length(list_names), function(x) PLSDA.Permut(mSetObj = mSet_list[[x]], permut.num, "bw")),
                 mSet_list[(length(list_names)+1):length(mSet_list)])

  names(mSet_list) <- c(list_names, names(mSet_list[(length(list_names)+1):length(mSet_list)]))

  mSet_list <- c(lapply(1:length(list_names), function(x) suppressWarnings(met.PLS_ImpScatter(mSetObj = mSet_list[[x]], imgName = paste0("Plots/PLSDA/PLS_Scatter_coef.mean_vs_VIP1_", names(mSet_list)[x]),
                                                                                              feat.nm = "coef.mean",  vip.nm=c("Comp. 1", "Comp. 2"),  dpi=72,
                                                                                              plot = F, export = export, format = export.format, show.title = TRUE,
                                                                                              title = paste0("(",
                                                                                                             str_replace(mSet_list[[x]][["dataSet"]][["trans.method"]],"N/A","none"),
                                                                                                             "/",
                                                                                                             str_replace(mSet_list[[x]][["dataSet"]][["scale.method"]],"N/A","none"),
                                                                                                             ")")))),
                 mSet_list[(length(list_names)+1):length(mSet_list)])

  for(i in levels(mSet_list[[1]][["dataSet"]][["cls"]])){
    names(mSet_list) <- c(list_names, names(mSet_list[(length(list_names)+1):length(mSet_list)]))
    mSet_list <- c(lapply(1:length(list_names), function(x) suppressWarnings(met.PLS_ImpScatter(mSetObj = mSet_list[[x]], imgName = paste0("Plots/PLSDA/PLS_Scatter_coef-", i, "_vs_VIP1_", names(mSet_list)[x]),
                                                                                                feat.nm = i,  vip.nm=c("Comp. 1", "Comp. 2"), dpi=72, vip.thresh = vip.thresh,
                                                                                                plot = F, export = export, format = export.format, show.title = TRUE,
                                                                                                title = paste0("(",
                                                                                                               str_replace(mSet_list[[x]][["dataSet"]][["trans.method"]],"N/A","none"),
                                                                                                               "/",
                                                                                                               str_replace(mSet_list[[x]][["dataSet"]][["scale.method"]],"N/A","none"),
                                                                                                               ")")))),
                   mSet_list[(length(list_names)+1):length(mSet_list)])
  }

  names(mSet_list) <- c(list_names, names(mSet_list[(length(list_names)+1):length(mSet_list)]))

  try(mSet_list <- c(lapply(1:length(list_names), function(x) met.plot_PLS.Classification(mSetObj = mSet_list[[x]],
                                                                                          imgName = paste0("Plots/PLSDA/PLSDA-CrossValidation_", names(mSet_list)[x]),
                                                                                          export = export, format = export.format, title = TRUE)),
                     mSet_list[(length(list_names)+1):length(mSet_list)]))

  mSet_list <- c(lapply(1:length(list_names), function(x) met.plot_PLS.Permutation(mSetObj = mSet_list[[x]],
                                                                                   imgName = paste0("Plots/PLSDA/PLSDA-Permutation_", names(mSet_list)[x]),
                                                                                   plot=F, export = export, format = export.format, title = TRUE)),
                 mSet_list[(length(list_names)+1):length(mSet_list)])

  names(mSet_list) <- c(list_names, names(mSet_list[(length(list_names)+1):length(mSet_list)]))

  if(!is.null(dev.list())){
    dev.off()
  }
  if(report==TRUE) {
    stack_size <- getOption("pandoc.stack.size", default = "512m")
    met.report_test_normalization(
      mSet_list,
      report.nm = report.nm,
      list_names = list_names,
      volcano = volcano,
      volcano.contrast = volcano.contrast,
      contrasts = cntrst,
      alpha = alpha,
      lfc = lfc,
      vip.thresh = vip.thresh,
      names_common = names_common
    )
  }
  setwd(wd)
  return(mSet_list)
}



$met.UpdateData
function (mSetObj = NA,
          filt.metab = c(""),
          filt.smpl = c(""),
          filt.grp = c(""))
{

  mSetObj$dataSet$edit <- NULL
  if (is.null(mSetObj$dataSet$filt) && !is.null(mSetObj$dataSet$data_proc)){
    data <- mSetObj$dataSet$data_proc
    cls <- mSetObj$dataSet$proc.cls
    if (substring(mSetObj$dataSet$format, 4, 5) == "ts") {
      facA <- mSetObj$dataSet$proc.facA
      facB <- mSetObj$dataSet$proc.facB
    }
  } else if (is.null(mSetObj$dataSet$filt) && !is.null(mSetObj$dataSet$preproc)) {
    data <- mSetObj$dataSet$preproc
    cls <- mSetObj$dataSet$proc.cls
    if (substring(mSetObj$dataSet$format, 4, 5) == "ts") {
      facA <- mSetObj$dataSet$proc.facA
      facB <- mSetObj$dataSet$proc.facB
    }
  } else {
    data <- mSetObj$dataSet$filt
    cls <- mSetObj$dataSet$filt.cls
    if (substring(mSetObj$dataSet$format, 4, 5) == "ts") {
      facA <- mSetObj$dataSet$filt.facA
      facB <- mSetObj$dataSet$filt.facB
    }
  }
  feat.hit.inx <- colnames(data) %in% filt.metab
  data <- CleanDataMatrix(data[, !feat.hit.inx, drop = FALSE])
  smpl.hit.inx <- rownames(data) %in% filt.smpl
  data <- CleanDataMatrix(data[!smpl.hit.inx, , drop = FALSE])
  cls <- as.factor(as.character(cls[!smpl.hit.inx]))
  if (substring(mSetObj$dataSet$format, 4, 5) == "ts") {
    facA <- as.factor(as.character(facA[!smpl.hit.inx]))
    facB <- as.factor(as.character(facB[!smpl.hit.inx]))
  }
  grp.hit.inx <- cls %in% filt.grp
  data <- CleanDataMatrix(data[!grp.hit.inx, , drop = FALSE])
  cls <- droplevels(factor(cls[!grp.hit.inx]))
  if (substring(mSetObj$dataSet$format, 4, 5) == "ts") {
    facA <- droplevels(factor(facA[!grp.hit.inx]))
    facB <- droplevels(factor(facB[!grp.hit.inx]))
  }
  print("Successfully updated the data!")
  mSetObj$dataSet$edit <- data
  mSetObj$dataSet$edit.cls <- cls
  if (substring(mSetObj$dataSet$format, 4, 5) == "ts") {
    mSetObj$dataSet$edit.facA <- facA
    mSetObj$dataSet$edit.facB <- facB
  }
  mSetObj <- met.PreparePrenormData(mSetObj)
  return(mSetObj)
}

$met.workflow
function(mSetObj,
         rowNorm = "NULL", # Option for row-wise normalization
         transNorm = "NULL", # Option to transform the data, "LogNorm" for Log Normalization, and "CrNorm" for Cubic Root Transformation.
         scaleNorm = "NULL", # Option for scaling the data, "MeanCenter" for Mean Centering, "AutoNorm" for Autoscaling, "ParetoNorm" for Pareto Scaling, amd "RangeNorm" for Range Scaling.
         ref = NULL, # Input the name of the reference sample or the reference feature, use " " around the name (for rowNorm="CompNorm" or "SpecNorm").
         contrast.type = "all",
         control = NULL,
         contrast = NULL,
         class_order = FALSE, # Class order matters (i.e. implying time points, disease severity, etc.)
         alpha = 0.05, # Significance threshold for adj. p values.
         lfc = 2, # Relevance threshold for log2(fold change) values.
         posthoc_method = "tukey",
         permut.num = 500,
         vip.thresh = 1,
         volcano.contrast = NULL,
         volcano.test = "anova",
         volcano.add_names = TRUE, # Display names next to symbols in volcano plot.
         volcano.label_size = 3, # Size of labels in volcano plot if
         volcano.adjusted = FALSE, # Display adjusted p-values on y axis of volcano plot? (only for volcano.test = "ttest"),
         pls.data = "all",
         plot = FALSE, # Shall plots be returned?
         export = FALSE, # Shall plots be exported as PDF and PNG files?
         export.format = "pdf",
         report = FALSE, # Shall a report (HTML and PDF) be created?
         report.nm = NULL, # Folder name for created report (if report = TRUE)
         export.dir = "met.ProcessedData"
)
{
  # Show error if inputs are not the required classes
  if(is.integer(alpha)) alpha <- as.numeric(alpha)
  if(is.integer(lfc)) lfc <- as.numeric(lfc)
  assertthat::assert_that(is.numeric(alpha),
                          length(alpha) == 1,
                          is.numeric(lfc),
                          length(lfc) == 1)
  if(!(rowNorm %in% c('NULL', 'QuantileNorm', 'CompNorm', 'SumNorm', 'MedianNorm', 'SpecNorm'))){
    stop(paste0(rowNorm, " is not a valid option for 'rowNorm'. Valid options are:\n",
                "'NULL', 'QuantileNorm', 'CompNorm', 'SumNorm', 'MedianNorm', 'SpecNorm'"))
  }
  if(!(transNorm %in% c('NULL', 'LogNorm', 'CrNorm'))){
    stop(paste0(transNorm, " is not a valid option for 'transNorm'. Valid options are:\n",
                "'NULL', 'LogNorm', 'CrNorm'"))
  }
  if(!(scaleNorm %in% c('NULL', 'MeanCenter', 'AutoNorm', 'ParetoNorm', 'RangeNorm'))){
    stop(paste0(scaleNorm, " is not a valid option for 'scaleNorm'. Valid options are:\n",
                "'NULL', 'MeanCenter', 'AutoNorm', 'ParetoNorm', 'RangeNorm'"))
  }
  dir.create(paste0(getwd(), "/", export.dir), showWarnings = F)
  wd <- getwd()
  file.copy("data_orig.qs", paste0(wd, "/", export.dir), overwrite = T)
  file.copy("preproc.qs", paste0(wd, "/", export.dir), overwrite = T)
  file.copy("prenorm.qs", paste0(wd, "/", export.dir), overwrite = T)
  file.copy("data_proc.qs", paste0(wd, "/", export.dir), overwrite = T)

  setwd(paste0(wd, "/", export.dir))

  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
    dir.create(paste0(getwd(), "/Plots/PLSDA"), showWarnings = F)
  }
  conditions <- levels(mSetObj$dataSet$cls)
  if (contrast.type == "all") {
    cntrst <- rev(apply(utils::combn(rev(conditions), 2), 2, paste,
                        collapse = "-"))
    if (!is.null(control)) {
      flip <- grep(paste("^", control, sep = ""),
                   cntrst)
      if (length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>% gsub(paste(control,
                                                    "- ", sep = " "), "", .) %>%
          paste(" - ", control, sep = "")
      }
    }
  }
  if (contrast.type == "control") {
    if (is.null(control))
      stop(
        paste0(
          "run test_diff(type = 'control') with a 'control' argument.\n",
          "Possible control conditions are: '",
          paste(conditions,
                collapse = "', '"),
          "'"
        )
      )
    if ((!control %in% conditions)) {
      stop(
        "Run met.workflow(contrast.type = 'control') with a valid condition in 'control'",
        ".\nValid controls are: '",
        paste0(conditions, collapse = "', '"),
        "'.",
        call. = FALSE
      )
    }
    cntrst <- paste(conditions[!conditions %in% control],
                    control, sep = " - ")
  }
  if (contrast.type == "manual") {
    if (is.null(contrast)) {
      stop(
        paste0(
          "Run met.workflow(contrast.type = 'manual') with a 'contrast' argument.\n",
          "Valid contrasts should contain combinations of: '",
          paste(conditions,
                collapse = "', '"
          ),
          "', for example '",
          paste0(conditions[1],
                 "_vs_", conditions[2]),
          "'."
        ) ,
        call. = FALSE)
    }
    assertthat::assert_that(is.character(contrast))
    if (any(!unlist(strsplit(contrast, "_vs_")) %in% conditions)) {
      stop("Run met.workflow() with valid contrasts in 'contrast'",
           ".\nValid contrasts should contain combinations of: '",
           paste0(conditions, collapse = "', '"),
           "', for example '", paste0(conditions[1],
                                      "_vs_", conditions[2]), "'.", call. = FALSE)
    }
    cntrst <- gsub("_vs_", " - ", contrast)
  }
  message("Tested contrasts: ", paste(gsub(" - ",
                                           "_vs_", cntrst), collapse = ", "))

  mSetObj[["analSet"]][["contrasts"]] <- cntrst

  mSetObj <- met.normalize(mSetObj, rowNorm = rowNorm, transNorm = transNorm, scaleNorm = scaleNorm, ref = ref)

  mSetObj <- met.ANOVA.Anal(mSetObj, nonpar = FALSE, thresh = alpha, post.hoc = posthoc_method, all_results = TRUE)
  mSetObj <- PCA.Anal(mSetObj)
  mSetObj <- met.PLSR.Anal(mSetObj, reg = class_order, data = pls.data)
  options(warn=1)
  mSetObj <- met.PLSDA.CV(mSetObj, methodName = "L", data = pls.data)
  options(warn=0)

  mSetObj[["analSet"]][["plsda"]][["class_order"]] <- class_order

  mSetObj <- met.plot_FeatureNormSummary(mSetObj = mSetObj, imgName = "Plots/norm_summary_features_", format = export.format,
                                         dpi = if_else(export.format=="pdf", 72, 300), width=NA, show_prenorm = TRUE,
                                         plot = plot, export = export)

  mSetObj <- met.plot_SampleNormSummary(mSetObj = mSetObj, imgName = "Plots/norm_summary_samples_", format = export.format,
                                        dpi = if_else(export.format=="pdf", 72, 300), width=NA, show_prenorm = TRUE,
                                        plot = plot, export = export)

  try(mSetObj <- met.plot_CorrHeatMap_Samples(mSetObj, imgName = "Plots/correlation_samples_", format = export.format,
                                              dpi = if_else(export.format=="pdf", 72, 300), plot = plot, export = export))

  mSetObj <- met.plot_CorrHeatMap_Features(mSetObj, imgName = "Plots/correlation_features_", format = export.format,
                                           dpi = if_else(export.format=="pdf", 72, 300), cor.method = "pearson",
                                           plot = plot, export = export)

  mSetObj <- met.plot_ANOVA(mSetObj, imgName = paste0("Plots/ANOVA-Plot"), format = export.format, subtitle = FALSE,
                            plot = F, export = export, dpi = if_else(export.format=="pdf", 72, 300))

  mSetObj <- met.plot_PCAScree(mSetObj, imgName = "Plots/PCA_ScreePlot_", format = export.format,
                               dpi = if_else(export.format=="pdf", 72, 300), width=NA, scree.num = 6,
                               plot = plot, export = export)

  mSetObj <- met.plot_PCA2DScore(mSetObj, imgName = "Plots/PCA_2DScores_PC1-PC2", format = export.format,
                                 dpi = if_else(export.format=="pdf", 72, 300), width=NA, pcx = 1, pcy = 2, reg = 0.95,
                                 show = 0, grey.scale = 0, plot = plot, export = export)

  mSetObj <- met.plot_PCA2DScore(mSetObj, imgName = "Plots/PCA_2DScores_PC2-PC3", format = export.format,
                                 dpi = if_else(export.format=="pdf", 72, 300), width=NA, pcx = 1, pcy = 3, reg = 0.95,
                                 show = 0, grey.scale = 0, plot = plot, export = export)

  mSetObj <- met.plot_PCA2DLoading(mSetObj, imgName = "Plots/PCA_LoadingPlot_", format = export.format,
                                 dpi = if_else(export.format=="pdf", 72, 300), width=NA, inx1 = 1, inx2 = 2,
                                 plot = plot, export = export)

  mSetObj <- met.plot_PCA2DLoading(mSetObj, imgName = "Plots/PCA_LoadingPlot_", format = export.format,
                                 dpi = if_else(export.format=="pdf", 72, 300), width=NA, inx1 = 1, inx2 = 3,
                                 plot = plot, export = export)

  mSetObj <- met.plot_PLS2DScore(mSetObj, "Plots/PLSDA/PLSDA_2DScores_C1-C2", format = export.format, dpi = if_else(export.format=="pdf", 72, 300),
                                 width = NA, inx1 = 1, inx2 = 2, reg = 0.95, show = 0, grey.scale = 0,
                                 plot = plot, export = export)

  mSetObj <- met.plot_PLS2DScore(mSetObj, "Plots/PLSDA/PLSDA_2DScores_C2-C3", format = export.format, dpi = if_else(export.format=="pdf", 72, 300),
                                 width = NA, inx1 = 1, inx2 = 3, reg = 0.95, show = 0, grey.scale = 0,
                                 plot = plot, export = export)

  mSetObj <- met.plot_PLS2DLoading(mSetObj, "Plots/PLSDA/PLSDA_LoadingPlot_", format = export.format, dpi = if_else(export.format=="pdf", 72, 300),
                                 width=NA, inx1 = 1, inx2 = 2, plot = plot, export = export)

  mSetObj <- met.plot_PLS2DLoading(mSetObj, "Plots/PLSDA/PLSDA_LoadingPlot_", format = export.format, dpi = if_else(export.format=="pdf", 72, 300),
                                 width=NA, inx1 = 1, inx2 = 3, plot = plot, export = export)

  mSetObj <- PLSDA.Permut(mSetObj, permut.num, "bw")

  mSetObj <- met.plot_PLS.Permutation(mSetObj, "Plots/PLSDA/PLS-DA_permutation_", format = export.format, dpi = if_else(export.format=="pdf", 72, 300), width=NA)

  mSetObj <- met.plot_PLS.Classification(mSetObj, imgName = paste0("Plots/PLSDA/PLSDA-CrossValidation"),
                                         export = export, format = export.format, title = TRUE)

  mSetObj <- met.plot_PLS_Imp(mSetObj, paste0("Plots/PLSDA/PLSDA_imp_Coeff_", "coef.mean"), format = export.format, dpi = if_else(export.format=="pdf", 72, 300), width=NA,
                              type = "coef", feat.nm = "coef.mean", feat.num = 20, color.BW = FALSE, export = export)

  for(i in unique(mSetObj$dataSet$cls)){
    mSetObj <- met.plot_PLS_Imp(mSetObj, paste0("Plots/PLSDA/PLSDA_imp_Coeff_", i), format = export.format, dpi = if_else(export.format=="pdf", 72, 300), width=NA,
                                type = "coef", feat.nm = i, feat.num = 20, color.BW = FALSE, export = export)
  }

  mSetObj <- met.plot_PLS_Imp(mSetObj, "Plots/PLSDA/PLSDA_imp_VIP_Comp1", format = export.format, dpi = if_else(export.format=="pdf", 72, 300), width=NA,
                              type = "vip", feat.nm = "Comp. 1", feat.num = 20, color.BW = FALSE, export = export)

  mSetObj <- met.plot_PLS_Imp(mSetObj, "Plots/PLSDA/PLSDA_imp_VIP_Comp2", format = export.format, dpi = if_else(export.format=="pdf", 72, 300), width=NA,
                              type = "vip", feat.nm = "Comp. 2", feat.num = 20, color.BW = FALSE, export = export)

  mSetObj <- met.PLS_ImpScatter(mSetObj, "Plots/PLSDA/PLSDA_Scatter_coef.mean_vs_VIP1",
                                feat.nm = "coef.mean",  vip.nm = c("Comp. 1", "Comp. 2"),  dpi=72,
                                plot = F, export = export, format = export.format)

  for (i in levels(mSetObj[["dataSet"]][["cls"]])) {
    mSetObj <-
      met.PLS_ImpScatter(mSetObj,
                         paste0("Plots/PLSDA/PLSDA_Scatter_coef-", i, "_vs_VIP1_"), feat.nm = i,  vip.nm = c("Comp. 1", "Comp. 2"),
                         dpi = 72, vip.thresh = vip.thresh, plot = F, export = export, format = export.format)
  }

  mSetObj <- met.plot_PCA3DScore(mSetObj, "Plots/PLSDA/PLSDA_3DScores", format = "json", export=export)
  mSetObj <- met.plot_PCA3DLoading(mSetObj, "Plots/PLSDA/PLSDA_3DScores", format = "json", export=export)
  mSetObj <- met.plot_PLS3DScore(mSetObj, "Plots/PLSDA/PLSDA_3DScores", format = "json", export=export)
  mSetObj <- met.plot_PLS3DLoading(mSetObj, "Plots/PLSDA/PLSDA_3DScores", format = "json", export=export)



  for(i in 1:length(cntrst)){
    grp1 <- gsub("-.+", "", cntrst[i])
    grp2 <- gsub(".+-", "", cntrst[i])
    mSetObj <- met.FC.Anal(mSetObj, log2fc.thresh = lfc, grp1 = grp1, grp2 = grp2, paired= FALSE)
    mSetObj <- met.plot_volcano(mSetObj, grp1 = grp2, grp2 = grp1,  label_size = volcano.label_size, test = volcano.test,
                                plot = plot, export = export, threshp = alpha, dpi = if_else(export.format=="pdf", 72, 300),
                                pval.type = if_else(volcano.adjusted==TRUE, "fdr", "raw"), format = export.format )
  }

  setwd(wd)
  if(report==TRUE){
    stack_size <- getOption("pandoc.stack.size", default = "512m")
    met.report(mSetObj, report.nm=report.nm, volcano=volcano, volcano.contrast = volcano.contrast, contrasts=cntrst, alpha=alpha, lfc=lfc, vip.thresh=vip.thresh)
  }
  return(mSetObj)


}
