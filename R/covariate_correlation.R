## Internal utilities for correlating gene/pathway expression with sample level covariates
## These are intentionally lightweight to avoid new dependencies.

#' Correlate each row (feature) of a matrix with a numeric vector.
#'
#' @param mat Numeric matrix (features x samples).
#' @param vec Numeric vector (length == ncol(mat)) containing a covariate.
#' @param feature_ids Character vector of feature identifiers (length == nrow(mat)).
#' @param method Correlation method: "spearman" (default) or "pearson".
#' @return A data.frame with columns: feature, cor, pvalue.
correlate_features_with_vector <- function(mat, vec, feature_ids = rownames(mat), method = c("spearman", "pearson")) {
  method <- match.arg(method)
  stopifnot(length(vec) == ncol(mat))
  # Remove samples with NA covariate
  keep <- which(!is.na(vec))
  vec2 <- vec[keep]
  mat2 <- if (length(keep) < ncol(mat)) mat[, keep, drop = FALSE] else mat
  if (length(unique(vec2)) < 3) {
    warning("Covariate has <3 unique non-NA values; correlations may be unreliable.")
  }
  res_list <- apply(mat2, 1, function(x){
    ct <- suppressWarnings(cor.test(x, vec2, method = method))
    c(cor = unname(ct$estimate), pvalue = ct$p.value)
  })
  res_df <- data.frame(feature = feature_ids, t(res_list), row.names = NULL, check.names = FALSE)
  res_df
}

#' High-level wrapper to correlate all features with multiple covariates.
#'
#' @param mat Numeric matrix (features x samples).
#' @param covariates Data.frame with rownames matching sample IDs (columns numeric or coercible).
#' @param method Correlation method.
#' @param adjust_method P-value adjustment method (BH by default) applied within each covariate.
#' @param feature_type Character label ("gene" or "pathway") used in output column naming.
#' @return Long-format data.frame: feature, variable, cor, pvalue, padj.
correlate_matrix_with_covariates <- function(mat, covariates, method = c("spearman","pearson"), adjust_method = "BH", feature_type = c("gene","pathway")) {
  method <- match.arg(method)
  feature_type <- match.arg(feature_type)
  # Keep only numeric columns
  num_cols <- vapply(covariates, is.numeric, logical(1))
  if (!any(num_cols)) stop("No numeric covariate columns provided.")
  cov_df <- covariates[, num_cols, drop = FALSE]
  # Reorder mat columns to match covariate rows if possible
  common_samples <- intersect(colnames(mat), rownames(cov_df))
  if (length(common_samples) < 3) stop("Fewer than 3 overlapping samples between expression matrix and covariates.")
  mat <- mat[, common_samples, drop = FALSE]
  cov_df <- cov_df[common_samples, , drop = FALSE]
  out_list <- vector("list", ncol(cov_df))
  for (j in seq_len(ncol(cov_df))) {
    v <- cov_df[[j]]
    cor_res <- correlate_features_with_vector(mat, v, feature_ids = rownames(mat), method = method)
    cor_res$variable <- colnames(cov_df)[j]
    # Adjust p-values within this variable
    cor_res$padj <- p.adjust(cor_res$pvalue, method = adjust_method)
    out_list[[j]] <- cor_res
  }
  out <- do.call(rbind, out_list)
  colnames(out)[1] <- feature_type
  out <- out[, c(feature_type, "variable", "cor", "pvalue", "padj")] 
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
  cor_df <- correlate_matrix_with_covariates(score_mat, covariates, method = method, adjust_method = adjust_method, feature_type = "pathway")
  cor_df$n_genes_used <- pathway_sizes[match(cor_df$pathway, names(pathway_sizes))]
  cor_df
}
