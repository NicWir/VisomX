## Internal utilities for correlating gene/pathway expression with sample level covariates
## These are intentionally lightweight to avoid new dependencies.

#' Correlate each row (feature) of a matrix with a numeric vector.
#'
#' @param mat Numeric matrix (features x samples).
#' @param vec Numeric vector (length == ncol(mat)) containing a covariate.
#' @param feature_ids Character vector of feature identifiers (length == nrow(mat)).
#' @param method Correlation method: "spearman" (default) or "pearson".
#' @param min_n Minimum number of non-missing pairs required to compute correlation (default 3).
#' @return A data.frame with columns: feature, cor, pvalue.
correlate_features_with_vector <- function(mat, vec, feature_ids = rownames(mat),
                                           method = c("spearman","pearson"),
                                           min_n = 6) {
  method <- match.arg(method)
  keep <- which(!is.na(vec))
  vec2 <- vec[keep]
  mat2 <- if (length(keep) < ncol(mat)) mat[, keep, drop = FALSE] else mat
  n <- ncol(mat2)

  if (identical(method, "pearson")) {
    if (n < min_n) {
      return(data.frame(feature = feature_ids, cor = NA_real_, pvalue = NA_real_, n = NA_real_,
                        row.names = NULL, check.names = FALSE))
    }
    v <- scale(vec2, center = TRUE, scale = TRUE)
    mat2_z <- t(scale(t(mat2), center = TRUE, scale = TRUE))
    r <- drop(mat2_z %*% v) / (n - 1)
    r[is.nan(r)] <- NA_real_
    n_per_row <- rowSums(!is.na(mat2))
    denom <- pmax(1e-12, 1 - r^2)
    df <- pmax(n_per_row - 2, 1)
    tval <- r * sqrt(pmax(0, n_per_row - 2) / denom)
    pval <- rep(NA_real_, length(r))
    ok <- n_per_row >= min_n & is.finite(r)
    pval[ok] <- 2 * stats::pt(abs(tval[ok]), df = df[ok], lower.tail = FALSE)
    r[!ok] <- NA_real_
    return(data.frame(feature = feature_ids, cor = as.numeric(r), pvalue = as.numeric(pval),
                      n = as.numeric(n_per_row), row.names = NULL, check.names = FALSE))
  }

  res_list <- apply(mat2, 1, function(x){
    ok <- which(!is.na(x) & !is.na(vec2)); n_ok <- length(ok)
    if (n_ok < min_n) return(c(cor = NA_real_, pvalue = NA_real_, n = NA_real_))
    ct <- suppressWarnings(cor.test(x[ok], vec2[ok], method = method))
    c(cor = unname(ct$estimate), pvalue = ct$p.value, n = n_ok)
  })
  data.frame(feature = feature_ids, t(res_list), row.names = NULL, check.names = FALSE)
}

#' High-level wrapper to correlate all features with multiple covariates.
#'
#' @param mat Numeric matrix (features x samples).
#' @param covariates Data.frame with rownames matching sample IDs (columns numeric or coercible).
#' @param method Correlation method.
#' @param adjust_method P-value adjustment method (BH by default) applied within each covariate.
#' @param feature_type Character label ("gene" or "pathway") used in output column naming.
#' @param var_filter_quantile Optional numeric quantile (0-1) to filter out low-variance features.
#'   E.g. 0.5 to keep only features with variance in the top 50%. Default NULL (no filtering).
#' @param min_sd Optional numeric minimum standard deviation threshold to filter out low-variance features.
#'   Default NULL (no filtering).
#' @param min_n Minimum number of non-missing pairs required; features with fewer pairs are set to NA (default 6).
#' @param low_n_threshold Minimum number of overlapping samples required to compute correlation (default 6).
#' @return Long-format data.frame: feature, variable, cor, pvalue, padj.
correlate_matrix_with_covariates <- function(mat, covariates,
                                             method = c("spearman","pearson"),
                                             adjust_method = "BH",
                                             feature_type = c("gene","pathway"),
                                             var_filter_quantile = NULL,
                                             min_sd = NULL,
                                             min_n = 6,
                                             low_n_threshold = 6) {
  method <- match.arg(method)
  feature_type <- match.arg(feature_type)

  # Keep only numeric columns, track categorical
  all_cols <- colnames(covariates)
  num_cols <- vapply(covariates, is.numeric, logical(1))
  categorical_cov <- all_cols[!num_cols]
  if (!any(num_cols)) stop("No numeric covariate columns provided.")
  cov_df <- covariates[, num_cols, drop = FALSE]

  # Align samples
  common_samples <- intersect(colnames(mat), rownames(cov_df))
  if (length(common_samples) < 3) stop("Fewer than 3 overlapping samples between expression matrix and covariates.")
  mat <- mat[, common_samples, drop = FALSE]
  cov_df <- cov_df[common_samples, , drop = FALSE]

  # Sanity: skip numeric covariates with <4 unique values; flag integer-like low-level as categorical-like
  is_integer_like <- function(x) {
    y <- x[!is.na(x)]
    if (!length(y)) return(FALSE)
    all(is.finite(y)) && all(abs(y - round(y)) < .Machine$double.eps^0.5)
  }
  num_cat_like <- names(which(vapply(cov_df, function(x) is_integer_like(x) && length(unique(stats::na.omit(x))) <= 8, logical(1))))
  low_unique <- names(which(vapply(cov_df, function(x) length(unique(stats::na.omit(x))) < 4, logical(1))))
  if (length(low_unique) > 0) {
    cov_df <- cov_df[, setdiff(colnames(cov_df), low_unique), drop = FALSE]
  }
  if (ncol(cov_df) == 0) {
    out <- data.frame()
    attr(out, "skipped_covariates_low_unique") <- low_unique
    attr(out, "categorical_covariates") <- categorical_cov
    attr(out, "numeric_categorical_covariates") <- num_cat_like
    return(out)
  }

  # Optional feature variance filters
  if (!is.null(var_filter_quantile) || !is.null(min_sd)) {
    vars <- MatrixGenerics::rowVars(mat, na.rm = TRUE)
    keep_var <- is.finite(vars)
    if (!is.null(min_sd)) keep_var <- keep_var & sqrt(vars) >= min_sd
    if (!is.null(var_filter_quantile)) {
      qcut <- stats::quantile(vars[keep_var], probs = var_filter_quantile, na.rm = TRUE)
      keep_var <- keep_var & vars >= qcut
    }
    mat <- mat[keep_var, , drop = FALSE]
  }

  out_list <- vector("list", ncol(cov_df))
  for (j in seq_len(ncol(cov_df))) {
    v <- cov_df[[j]]
    cor_res <- correlate_features_with_vector(mat, v, feature_ids = rownames(mat),
                                              method = method, min_n = min_n)
    cor_res$variable <- colnames(cov_df)[j]
    cor_res$padj <- p.adjust(cor_res$pvalue, method = adjust_method)
    out_list[[j]] <- cor_res
  }
  out <- do.call(rbind, out_list)
  if (!is.null(low_n_threshold) && "n" %in% colnames(out)) {
    out$low_n <- !is.na(out$n) & out$n < low_n_threshold
  }
  colnames(out)[1] <- feature_type
  keep_cols <- c(feature_type, "variable", "cor", "pvalue", "padj", "n", "low_n")
  out <- out[, intersect(keep_cols, colnames(out))]

  attr(out, "skipped_covariates_low_unique") <- low_unique
  attr(out, "categorical_covariates") <- categorical_cov
  attr(out, "numeric_categorical_covariates") <- num_cat_like
  out
}

#' Compute simple pathway scores (mean expression) and correlate with covariates.
#'
#' @param expr_mat Gene expression matrix (genes x samples) from rlog/vst.
#' @param pathway_list Named list of character vectors of gene IDs.
#' @param covariates Data.frame of sample covariates (rownames = sample IDs).
#' @param min_genes Minimum genes required per pathway to compute a score (default 5).
#' @param method Correlation method.
#' @param adjust_method P-value adjustment method.
#' @return data.frame pathway correlations (pathway, variable, cor, pvalue, padj, n_genes_used).
correlate_pathways <- function(expr_mat, pathway_list, covariates, min_genes = 5, method = c("spearman","pearson"), adjust_method = "BH") {
  method <- match.arg(method)
  # Build pathway score matrix
  pw_scores <- lapply(names(pathway_list), function(pw){
    genes <- intersect(pathway_list[[pw]], rownames(expr_mat))
    if (length(genes) < min_genes) return(NULL)
    sc <- colMeans(expr_mat[genes, , drop = FALSE])
    attr(sc, "n_genes") <- length(genes)
    sc
  })
  names(pw_scores) <- names(pathway_list)
  pw_scores <- pw_scores[!vapply(pw_scores, is.null, logical(1))]
  if (length(pw_scores) == 0) return(data.frame())
  score_mat <- do.call(rbind, lapply(pw_scores, function(x) x))
  pathway_sizes <- vapply(pw_scores, function(x) attr(x, "n_genes"), numeric(1))
  # Correlate pathway scores with covariates
  cor_df <- correlate_matrix_with_covariates(score_mat, covariates,
                                             method = method,
                                             adjust_method = adjust_method,
                                             feature_type = "pathway",
                                             min_n = 6,
                                             low_n_threshold = 6)
  cor_df$n_genes_used <- pathway_sizes[match(cor_df$pathway, names(pathway_sizes))]
  cor_df
}

#' Run GSEA on gene–covariate correlation ranks
run_gsea_on_correlation <- function(gene_cor_df, id_map, gmt_file,
                                    pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                    minGSSize = 5, maxGSSize = 500) {
  if (is.null(gene_cor_df) || nrow(gene_cor_df) == 0) return(list())
  if (is.null(gmt_file) || !file.exists(gmt_file)) stop("GMT file not found: ", gmt_file)
  gmt_df <- clusterProfiler::read.gmt(gmt_file)
  term2gene <- data.frame(term = gmt_df$term, gene = gmt_df$gene, stringsAsFactors = FALSE)
  out <- list()
  for (v in unique(gene_cor_df$variable)) {
    sub <- gene_cor_df[gene_cor_df$variable == v, c("gene","cor")]
    sub$id <- unname(id_map[match(sub$gene, names(id_map))])
    gl <- sub$cor; names(gl) <- sub$id
    gl <- gl[!is.na(names(gl)) & !is.na(gl)]
    if (length(gl) < minGSSize) { out[[v]] <- NULL; next }
    gl <- sort(gl, decreasing = TRUE)
    out[[v]] <- tryCatch(
      clusterProfiler::GSEA(geneList = gl, TERM2GENE = term2gene,
                            pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff,
                            minGSSize = minGSSize, maxGSSize = maxGSSize,
                            verbose = FALSE, seed = TRUE, by = "fgsea"),
      error = function(e) NULL
    )
  }
  out
}
