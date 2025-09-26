#' Multi-omics RNA/Proteomics workflow for VisomX
#'
#' @description
#' `multiomics.workflow()` coordinates the existing RNA sequencing (`rna.workflow`) and
#' proteomics (`prot.workflow`) pipelines, optionally aligning samples, sharing sample-level
#' covariates, and augmenting the proteomics results with the covariate correlation analyses
#' that were previously only available for RNA data. The helper makes it easy to run both
#' assays side-by-side and collect shared metadata such as overlapping samples or contrasts.
#'
#' @param rna_se A `SummarizedExperiment` containing RNA-seq counts prepared for
#'   `rna.workflow()`.
#' @param prot_se A `SummarizedExperiment` containing proteomics intensities prepared for
#'   `prot.workflow()`.
#' @param sample_covariates Optional `data.frame` of sample-level covariates (rownames must
#'   correspond to sample identifiers). Numeric columns are used for covariate correlation
#'   analyses.
#' @param match_samples Either "intersection" (default) to restrict both assays and
#'   covariates to the overlapping sample names, or "none" to run the workflows on the
#'   supplied objects without enforcing overlap.
#' @param rna_args Named list of additional arguments forwarded to `rna.workflow()`.
#'   Entries in this list override defaults set by `multiomics.workflow()`.
#' @param prot_args Named list of additional arguments forwarded to `prot.workflow()`.
#' @param prot_covariate Logical (default `TRUE`); if `TRUE` and covariates are supplied,
#'   the function forwards them to `prot.workflow()` to enable protein/pathway
#'   correlation analyses.
#' @param prot_covariate_options Named list whose entries override the corresponding
#'   arguments in `prot.workflow()` (e.g. `correlate_proteins`, `correlate_pathways`,
#'   `pathway_min_genes`, `correlation_method`, `correlation_p_adjust`,
#'   `correlation_var_filter`, `correlation_min_sd`, `correlation_low_n`,
#'   `gsea_covariates`).
#' @param quiet Logical; passed through to `rna.workflow()` when not already supplied via
#'   `rna_args`.
#'
#' @return A list with RNA and proteomics workflow outputs, enriched with shared metadata.
#'   The object contains
#'   
#'   * `rna`: Output from `rna.workflow()`.
#'   * `proteomics`: Output from `prot.workflow()` (with optional covariate correlations).
#'   * `samples`: Sample name bookkeeping (`shared`, `rna`, `proteomics`, `rna_only`,
#'     `proteomics_only`).
#'   * `sample_covariates`: The (potentially subset) covariate table used for analyses.
#'   * `shared_contrasts`: Intersection of contrast identifiers between the RNA and
#'     proteomics analyses when available.
#'   * `metadata`: Additional run metadata (matching mode, whether covariate analysis was
#'     executed for proteomics).
#'
#' @export
#'
multiomics.workflow <- function(
  rna_se,
  prot_se,
  sample_covariates = NULL,
  match_samples = c("intersection", "none"),
  rna_args = list(),
  prot_args = list(),
  prot_covariate = TRUE,
  prot_covariate_options = list(),
  quiet = FALSE
) {
  match_mode <- match.arg(match_samples)

  if (!inherits(rna_se, "SummarizedExperiment")) {
    stop("'rna_se' must be a SummarizedExperiment.", call. = FALSE)
  }
  if (!inherits(prot_se, "SummarizedExperiment")) {
    stop("'prot_se' must be a SummarizedExperiment.", call. = FALSE)
  }
  if (!is.null(sample_covariates) && !is.data.frame(sample_covariates)) {
    stop("'sample_covariates' must be a data.frame or NULL.", call. = FALSE)
  }
  if (!is.list(rna_args)) stop("'rna_args' must be a list.", call. = FALSE)
  if (!is.list(prot_args)) stop("'prot_args' must be a list.", call. = FALSE)

  rna_samples <- colnames(rna_se)
  prot_samples <- colnames(prot_se)
  shared_samples <- intersect(rna_samples, prot_samples)

  if (match_mode == "intersection" && length(shared_samples) == 0) {
    stop("No overlapping samples between RNA and proteomics inputs when match_samples='intersection'.", call. = FALSE)
  }

  rna_use <- rna_se
  prot_use <- prot_se
  cov_use <- sample_covariates

  if (match_mode == "intersection") {
    rna_use <- rna_use[, rna_samples %in% shared_samples, drop = FALSE]
    prot_use <- prot_use[, prot_samples %in% shared_samples, drop = FALSE]
    if (!is.null(cov_use)) {
      cov_use <- cov_use[intersect(rownames(cov_use), shared_samples), , drop = FALSE]
    }
  } else {
    if (!is.null(cov_use)) {
      valid_rows <- intersect(rownames(cov_use), union(rna_samples, prot_samples))
      cov_use <- cov_use[valid_rows, , drop = FALSE]
    }
  }

  # Prepare argument lists for underlying workflows -----------------------------------
  rna_args$se <- rna_use
  if (!"quiet" %in% names(rna_args)) {
    rna_args$quiet <- quiet
  }
  if (!is.null(cov_use)) {
    if (!"sample_covariates" %in% names(rna_args)) {
      rna_args$sample_covariates <- cov_use[intersect(rownames(cov_use), colnames(rna_use)), , drop = FALSE]
    } else {
      rna_args$sample_covariates <- rna_args$sample_covariates[intersect(rownames(rna_args$sample_covariates), colnames(rna_use)), , drop = FALSE]
    }
  }

  rna_res <- do.call(rna.workflow, rna_args)

  prot_args$se <- prot_use
  if (!"quiet" %in% names(prot_args)) {
    prot_args$quiet <- quiet
  }
  if (prot_covariate && !is.null(cov_use)) {
    overlap <- intersect(rownames(cov_use), colnames(prot_use))
    if (length(overlap) < 3) {
      warning("<3 overlapping samples between proteomics data and 'sample_covariates'; skipping proteomics covariate analysis.")
    } else {
      cov_subset <- cov_use[overlap, , drop = FALSE]
      if (!"sample_covariates" %in% names(prot_args)) {
        prot_args$sample_covariates <- cov_subset
      } else {
        prot_args$sample_covariates <- prot_args$sample_covariates[intersect(rownames(prot_args$sample_covariates), overlap), , drop = FALSE]
      }
      # defaults & overrides from prot_covariate_options
      opts <- prot_covariate_options
      if (!"correlate_proteins" %in% names(prot_args)) {
        prot_args$correlate_proteins <- if (!is.null(opts$correlate_proteins)) opts$correlate_proteins else TRUE
      }
      override_keys <- c("correlate_pathways", "pathway_min_genes", "correlation_method",
                         "correlation_p_adjust", "correlation_var_filter", "correlation_min_sd",
                         "correlation_low_n", "gsea_covariates")
      for (nm in override_keys) {
        if (!(nm %in% names(prot_args)) && !is.null(opts[[nm]])) {
          prot_args[[nm]] <- opts[[nm]]
        }
      }
    }
  }

  prot_res <- do.call(prot.workflow, prot_args)

  shared_contrasts <- intersect(names(rna_res$results), infer_proteomics_contrasts(prot_res))

  out <- list(
    rna = rna_res,
    proteomics = prot_res,
    samples = list(
      rna = colnames(rna_use),
      proteomics = colnames(prot_use),
      shared = intersect(colnames(rna_use), colnames(prot_use)),
      rna_only = setdiff(colnames(rna_se), colnames(prot_use)),
      proteomics_only = setdiff(colnames(prot_se), colnames(rna_use))
    ),
    sample_covariates = cov_use,
    shared_contrasts = shared_contrasts,
    metadata = list(
      match_samples = match_mode,
      prot_covariate_analysis = prot_covariate && !is.null(cov_use),
      prot_covariate_overlap = prot_res$covariate_overlap
    )
  )

  if (!is.null(prot_res$protein_covariate_correlations)) {
    out$proteomics_covariate_correlations <- prot_res$protein_covariate_correlations
  }
  if (!is.null(prot_res$pathway_covariate_correlations)) {
    out$proteomics_pathway_covariate_correlations <- prot_res$pathway_covariate_correlations
  }
  if (!is.null(prot_res$gsea_covariate_results)) {
    out$proteomics_covariate_gsea <- prot_res$gsea_covariate_results
  }

  return(out)
}

#' Infer proteomics contrast identifiers from workflow output
#'
#' @param prot_res Output list from `prot.workflow()`.
#'
#' @return Character vector of contrast names inferred from result column names.
#'
#' @keywords internal
#'
infer_proteomics_contrasts <- function(prot_res) {
  if (is.null(prot_res$results) || !is.data.frame(prot_res$results)) {
    return(character())
  }
  fc_cols <- grep("_log2fc$", colnames(prot_res$results), value = TRUE)
  gsub("_log2fc$", "", fc_cols)
}
