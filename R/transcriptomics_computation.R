#' @title Read transcriptomics data
#' @description Imports count tables (CSV/TSV/XLS/XLSX) or data frames and returns a `SummarizedExperiment`, with optional RSD and missing-value filtering.
#' @param data Data frame or file path for a count table.
#' @param files.ind Character; prefix for matching multiple files in the working directory.
#' @param expdesign Experimental design as data frame or file path.
#' @param csvsep Field separator for CSV files (default `";"`).
#' @param dec Decimal separator for CSV/TSV/TXT files (default `"."`).
#' @param sheet Sheet number or name for Excel files.
#' @param name Column name for primary gene IDs.
#' @param id Column name for alternative gene IDs.
#' @param values Column name for count values when using `files.ind`.
#' @param id2name.table Optional ID→name mapping table or file path.
#' @param id2name.id Column name of alternative IDs in mapping table.
#' @param id2name.name Column name of primary IDs in mapping table.
#' @param pfx.counts Prefix for abundance column names.
#' @param rsd_thresh Numeric; RSD threshold (%) for filtering genes.
#' @param filt_type Character; one of `"complete"`, `"condition"`, or `"fraction"` for missing-value filtering.
#' @param filt_thr Numeric; threshold for missing-value count per condition.
#' @param filt_min Numeric; minimum fraction of valid values required.
#' @return A `SummarizedExperiment` object with filtered and processed counts.
#' @export
#' @importFrom assertthat assert_that
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom dplyr select filter rename
#' @importFrom magrittr %>%
rna.read_data <- function (data = NULL,
                           files.ind = NULL,
                           expdesign = NULL,
                           csvsep = ";",
                           dec = ".",
                           sheet = 1,
                           name = "SymbolID",
                           id = 'gene_id',
                           values = "FPKM",
                           id2name.table = NULL,
                           id2name.id = NULL,
                           id2name.name = NULL,
                           pfx.counts = "counts.",
                           rsd_thresh = NULL,
                           filt_type =  c("condition", "complete", "fraction", NULL),
                           filt_thr = 3,
                           filt_min = NULL
) {
  assertthat::assert_that(is.character(name),
                          length(name) == 1,
                          is.character(id),
                          length(id) == 1)
  ## --- (1) Read data -----------------------------------------------------------
  if(!is.null(data)){
    if (is.character(data)) {
      dat <- read_file(data, csvsep = csvsep, dec = dec, sheet = sheet)
    } else if (is.data.frame(data) || is.matrix(data)) {
      dat <- data
    } else {
      stop("`data` must be a file path, data.frame, or matrix.")
    }
  } else {
    # Combine separate data files into a single dataframe
    if(!is.null(files.ind)){
      if(!is.character(files.ind)){
        stop("files.ind needs to be a string or vector of several strings.")
      } else {
        pfx.counts <- values
        dfs <- list()
        files <- list.files()[grep(files.ind, list.files())]
        for(i in seq_along(files)){
          dfs[[i]] <- read_file(files[i], csvsep = csvsep)
          names(dfs)[i] <- gsub(files.ind, "", files[i])
        }
        # Join dataframes in list
        df_joined <- data.frame(id = NA,
                                name = NA,
                                "Values" = NA)
        names(df_joined)[1] <- paste0(id)
        names(df_joined)[2] <- paste0(name)
        names(df_joined)[3] <- paste0(values)
        for(i in seq_along(dfs)){
          df_joined <- merge(df_joined, dfs[[i]],
                             by = c(names(df_joined)[1], names(df_joined)[2]),
                             all.x = TRUE, all.y = TRUE)
          colnames(df_joined) <- gsub(pattern = paste0(values, ".y"),
                                      replacement = paste0(pfx.counts, names(dfs)[i]),
                                      colnames(df_joined))
          colnames(df_joined) <- gsub(pattern = paste0(values, "$"),
                                      replacement = paste0(pfx.counts, names(dfs)[i]),
                                      colnames(df_joined))
          df_joined <- df_joined %>% select(grep(paste(values, id, name, sep = "|"), colnames(df_joined)))
        }
        dat <- df_joined[, -grep(paste0(values, ".x"), colnames(df_joined))]
      }
    } else {
      stop("Please run 'rna.read_data' with an R dataframe object or a table file containing read data for all samples in the 'data' argument, or prefixes of several files in the 'files.ind' argument.")
    }
  }

  ## --- (2) Fix column headers if needed --------------------------------------
  # if (all(grepl("^V\\d", colnames(dat)))){
  #   new_colnames <- gsub('["\']', '', dat[1,])
  #   colnames(dat) <- new_colnames
  #   dat <- dat[-1,]
  # }
  # dat <- dat[, !is.na(colnames(dat))]

  # ---(2) Ensure first row is used as column headers if needed)------

  if (all(grepl("^V\\d", colnames(dat)))){
    new_colnames <- gsub('["\']', '', dat[1,])
    colnames(dat) <- new_colnames
    dat <- dat[-1,]
  }

  ## --- (3) Optionally read id2name table --------------------------------------
  if (is.character(id2name.table)) {
    id2name <- read_file(id2name.table, csvsep = csvsep)
  }

  ## --- (4) Check that the abundance prefix exists -----------------------------
  if (!(any(grepl(pfx.counts, colnames(dat))))) {
    stop(paste0("The prefix '", pfx.counts, "' does not exist in any column of '",
                data, "'. Please provide a valid prefix to identify columns with gene abundances."), call. = F)
  }

  ## --- (5) Rename duplicate abundance columns -------------------------------
  suffix_all_duplicates <- function(x, sep = "_") {
    counts <- ave(seq_along(x), x, FUN = seq_along)
    paste0(x, sep, counts)
  }

  abundance_idx <- grep(paste0("^", pfx.counts), colnames(dat))
  if (length(abundance_idx) > 0) {
    ab_names <- colnames(dat)[abundance_idx]
    ab_names_unique <- suffix_all_duplicates(ab_names)
    colnames(dat)[abundance_idx] <- ab_names_unique
  }

  ## --- (6) Make row-wise gene names unique ------------------------------------
  if (any(colnames(dat) %in% name) && any(colnames(dat) %in% id)) {
    row.rm <- which(is.na(dat[name]) & is.na(dat[id]), arr.ind = TRUE)
    if (length(row.rm) > 1) {
      dat.rm <- dat[-row.rm, ]
      message(paste0("Removing genes with no defined name in column \"", name,
                     "\" or \"", id, "\". (Removed: ",
                     (nrow(dat) - nrow(dat.rm)), " genes)"))
    }

    if (length(row.rm) > 1) {
      rna_unique <- prot.make_unique(proteins = dat.rm, names = name, ids = id, delim = ";")
    } else {
      rna_unique <- prot.make_unique(proteins = dat, names = name, ids = id, delim = ";")
    }

  } else if (!any(colnames(dat) %in% name) && !any(colnames(dat) %in% id)) {
    stop("\"", name, "\" and \"", id, "\" are not columns in \"", data, "\". Please provide valid column names with gene identifiers as \"name=\" and \"id=\".", call. = F)
  } else if (!any(colnames(dat) %in% name)) {
    stop("\"", name, "\" is not a column in \"", data, "\". Please provide a valid column name with gene names as \"name=\".", call. = F)
  } else if (!any(colnames(dat) %in% id)) {
    stop("\"", id, "\" is not a column in \"", data, "\". Please provide a valid column name with gene identifiers as \"id=\".", call. = F)
  }
  ## --- (7) Build experimental design if none is provided -----------------------
  if (is.null(expdesign)) {
    label <- rna_unique %>%
      dplyr::select(contains(pfx.counts)) %>%
      colnames() %>%
      gsub(pfx.counts, "", .)

    condition <- rna_unique %>%
      dplyr::select(contains(pfx.counts)) %>%
      colnames() %>%
      gsub(pfx.counts, "", .) %>%
      gsub(".[[:digit:]]+$", "", .)

    replicate <- rna_unique %>%
      dplyr::select(contains(pfx.counts)) %>%
      colnames() %>%
      stringr::str_extract("[:digit:]{1,}$")

    experimental_design <- data.frame(label, condition, replicate, check.names = FALSE)
    message(paste0("Writing experimental design file to: ", getwd(), "/experimental_design.txt"))
    utils::write.table(experimental_design,
                       file = "experimental_design.txt",
                       sep = "\t",
                       row.names = FALSE)
  } else if (is.list(expdesign)) {
    experimental_design <- expdesign
  } else if (is.character(expdesign) && file.exists(expdesign)) {
    ext <- str_replace_all(expdesign, ".{1,}\\.", "")
    if (ext == "csv") {
      experimental_design <- utils::read.csv(expdesign, sep = csvsep, header = TRUE,
                                             stringsAsFactors = FALSE, fill = TRUE, na.strings = "",
                                             quote = "", comment.char = "", check.names = FALSE)
    } else if (ext %in% c("xls", "xlsx")) {
      experimental_design <- readxl::read_excel(expdesign)
    } else if (ext == "tsv") {
      experimental_design <- utils::read.csv(expdesign, sep = "\t", header = TRUE,
                                             stringsAsFactors = FALSE, fill = TRUE, na.strings = "",
                                             quote = "", comment.char = "", check.names = FALSE)
    } else if (ext == "txt") {
      experimental_design <- utils::read.table(expdesign, sep = "\t", header = TRUE,
                                               stringsAsFactors = FALSE, fill = TRUE, na.strings = "", comment.char = "", check.names = FALSE)
    } else {
      stop("No compatible file format for experimental design provided.
           Supported formats are: .txt (tab delimited), .csv, .tsv, .xls, and .xlsx")
    }
  } else {
    stop(paste0("File \"", expdesign, "\" does not exist."), call. = FALSE)
  }

  ## --- (8) Create SummarizedExperiment ----------------------------------------
  abundance_columns <- grep(pfx.counts, colnames(rna_unique))

  # Convert abundance columns to numeric with proper NA handling
  if (any(sapply(rna_unique[, abundance_columns], function(x) is.character(x)))) {
    rna_unique[, abundance_columns] <- apply(rna_unique[, abundance_columns], 2,
                                             function(x) as.numeric(dplyr::na_if(x, "NA")))
  }

  message("Generating SummarizedExperiment.")
  rna_se <- rna.make_se(rna_unique, abundance_columns, experimental_design)

  ## --- (9) Apply RSD threshold, if any ----------------------------------------
  if(!is.null(rsd_thresh)){
    rsd_thresh <- rsd_thresh/100
    int.mat <- SummarizedExperiment::assay(rna_se)
    int.mat[is.na(int.mat)] <- 0

    ls.sd <- lapply(seq_along(unique(rna_se$condition)), function (x)
      apply(int.mat[, unique(rna_se$condition)[x] == gsub("_[[:digit:]]+$", "", colnames(int.mat))],
            1, stats::sd, na.rm = FALSE))
    ls.mns <- lapply(seq_along(unique(rna_se$condition)), function (x)
      apply(int.mat[, unique(rna_se$condition)[x] == gsub("_[[:digit:]]+$", "", colnames(int.mat))],
            1, stats::mean, na.rm = FALSE))
    ls.rsd <- mapply("/", ls.sd, ls.mns, SIMPLIFY = FALSE)
    mat.rsd <- do.call(rbind, ls.rsd)
    filter.val <- apply(mat.rsd, 2, max, na.rm = TRUE)
    keep <- filter.val <= rsd_thresh
    filtered_rsd <- rna_se[keep, ]
  } else {
    filtered_rsd <- rna_se
  }

  if ((nrow(SummarizedExperiment::assay(rna_se)) - nrow(SummarizedExperiment::assay(filtered_rsd))) != 0) {
    number_removed_rsd <- nrow(SummarizedExperiment::assay(rna_se)) - nrow(SummarizedExperiment::assay(filtered_rsd))
    cat(paste0(number_removed_rsd, " out of ",
               nrow(SummarizedExperiment::assay(rna_se)), " genes were removed from the dataset due to too high RSD.\n\n"))
    filtered_rsd@metadata$n.filtered.rsd <- number_removed_rsd
  }
  filtered_rsd@metadata$rsd.thresh <- rsd_thresh
  # record key parameters in metadata for reproducibility
  meta <- filtered_rsd@metadata
  if (is.null(meta) || !is.list(meta)) meta <- list()
  meta$filt.type   <- if (is.null(filt_type)) NA else filt_type[1]
  meta$filt.thr    <- filt_thr
  meta$filt.min    <- filt_min
  meta$pfx.counts  <- pfx.counts
  meta$csvsep      <- csvsep
  meta$dec         <- dec
  meta$sheet       <- sheet
  meta$name.col    <- name
  meta$id.col      <- id
  filtered_rsd@metadata <- meta

  ## --- (10) Remove samples (columns) with all NA values -----------------------
  assay_data <- SummarizedExperiment::assay(filtered_rsd)
  na_cols <- which(colSums(!is.na(assay_data)) == 0)
  na_sample_names <- colnames(assay_data)[na_cols]
  if (length(na_sample_names) > 0) {
    cat("Samples with all genes as NA:\n")
    print(na_sample_names)
  } else {
    cat("No samples have all genes as NA.\n")
  }

  if (length(na_cols) > 0) {
    filtered_rsd <- filtered_rsd[, -na_cols]
    cat("Removed", length(na_cols), "samples with all NA values.\n")
  } else {
    cat("No samples to remove.\n")
  }

  ## --- (11) Filter missing values based on specified criteria ------------------
  if(!is.null(filt_type)) {
    rna_se <- rna.filter_missing(filtered_rsd, type = filt_type, thr = filt_thr, min = filt_min)
  } else {
    rna_se <- filtered_rsd
  }

  # carry over metadata to final SummarizedExperiment
  if (is.list(filtered_rsd@metadata)) {
    rna_se@metadata <- filtered_rsd@metadata
  }

  cat(paste0("Identified conditions:\n ", paste(unique(rna_se$condition), collapse = ", "), "\n"))

  return(rna_se)
}

#' Filter genes based on missing values
#'
#' \code{rna.filter_missing} filters a transcriptomics dataset based on missing values. Different types of filtering can be applied, which range from only keeping genes without missing values to keeping genes with a certain percent valid values in all samples or keeping genes that are complete in at least one condition.
#'
#' @param se \code{SummarizedExperiment} object, transcriptomics data parsed with \code{\link{rna.read_data}}.
#' @param type (Character string) "complete", "condition" or "fraction", Sets the type of filtering applied. "complete" will only keep genes with valid values in all samples. "condition" will keep genes that have a maximum of \code{thr} missing values in at least one condition. "fraction" will keep genes that have a \code{min} fraction of valid values in all samples.
#' @param thr (Integer) Sets the threshold for the allowed number of missing values in at least one condition if \code{type = "condition"}. In other words: "keep genes that have a maximum of 'thr' missing values in at least one condition."
#' @param min (Numeric) Sets the threshold for the minimum fraction of valid values allowed for any gene if \code{type = "fraction"}.
#'
#' @return A filtered SummarizedExperiment object.
#' @export
rna.filter_missing <- function (se, type = c("complete", "condition", "fraction", NULL),
                                 thr = NULL, min = NULL)
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  type <- match.arg(type)
  if (any(!c("name", "ID") %in% colnames(SummarizedExperiment::rowData(se,
                                                                       use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)), "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in%
          colnames(SummarizedExperiment::colData(se)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)), "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  if (!is.null(type)){
    if (type == "complete") {
      keep <- !apply(SummarizedExperiment::assay(se), 1, function(x) any(is.na(x)))
      filtered <- se[keep, ]
    }
    if (type == "condition") {
      assertthat::assert_that(is.numeric(thr), length(thr) ==
                                1)
      max_repl <- max(SummarizedExperiment::colData(se)$replicate)
      if (thr < 0 | thr > max_repl) {
        stop("invalid filter threshold 'thr' applied",
             "\nRun filter() with a threshold ranging from 0 to the maximum number of replicates: ",
             max_repl)
      }
      filtered <- filter_missval(se, thr = thr)
    }
    if (type == "fraction") {
      assertthat::assert_that(is.numeric(min), length(min) ==
                                1)
      if (min < 0 | min > 1) {
        stop("invalid filter threshold 'min' applied",
             "\nRun filter() with a percent ranging from 0 to 1")
      }
      # bin_data <- SummarizedExperiment::assay(se)
      # idx <- is.na(SummarizedExperiment::assay(se))
      # bin_data[!idx] <- 1
      # bin_data[idx] <- 0
      # keep <- bin_data %>% as.data.frame() %>% tibble::rownames_to_column() %>%
      #   gather(ID, value, -rowname) %>% group_by(rowname) %>%
      #   dplyr::summarize(n = n(), valid = sum(value), frac = valid/n) %>%
      #   filter(frac >= min)
      # filtered <- se[keep$rowname, ]
      mat <- SummarizedExperiment::assay(se)
      keep <- rowMeans(!is.na(mat)) >= min
      filtered <- se[keep, ]
    }
  } else {
    filtered <- se
  }
  if( (nrow(SummarizedExperiment::assay(se)) - nrow(SummarizedExperiment::assay(filtered))) != 0 ){
    number_removed <- nrow(SummarizedExperiment::assay(se)) - nrow(SummarizedExperiment::assay(filtered))
    cat(paste0(number_removed, " out of ",
               nrow(SummarizedExperiment::assay(se)), " genes were removed from the dataset due to missing values.\n\n"))

    # Get indices of removed rows
    removed_idx <- base::setdiff(rownames(se), rownames(filtered))

    # Randomly pick up to 10
    sample_rows <- sample(removed_idx, min(10, length(removed_idx)))

    # Extract assay values for those rows
    preview <- SummarizedExperiment::assay(se)[sample_rows, , drop = FALSE]

    cat("Random sample of removed genes:\n")
    print(preview)

    filtered@metadata$n.filtered <- number_removed
  }
  filtered@metadata$filt_type <- type
  filtered@metadata$n.pre_filt <- nrow(SummarizedExperiment::assay(se))
  return(filtered)
}



#' Generate a transcriptomics analysis report
#'
#' This function generates a report based on the results of an RNA analysis workflow.
#'
#' @param results A list object containing the results of an RNA analysis workflow.
#' @param report.dir A character string specifying the directory in which to save the report.
#' @param ... Additional parameters to be passed to the report.
#' @return A report in the specified directory.
#' @export
#'
#' @examples
#'
#' rna.report(results, report.dir = NULL)
rna.report <- function(res, ..., report.dir = NULL, param = NULL) {
  results <- get0("results", ifnotfound = get0("res", ifnotfound = NULL))
  assertthat::assert_that(is.list(results))
  if (any(!c("data", "se",
             "imputed", "dds", "pca", "results",
             "param") %in% names(results))) {
    stop("run rna.report() with appropriate input generated by rna.workflow",
         call. = FALSE)
  }
  args <- list(...)
  for(i in 1:length(args)){
    assign(names(args)[i], args[[i]])
  }
  data <- results$se
  norm <- results$norm
  imp <- results$imputed
  pca <- results$pca
  dds <- results$dds
  #param <- results$param
  res <- results

  sample_covariates <- if (exists("sample_covariates")) get("sample_covariates") else NULL
  if (!is.null(sample_covariates)) {
    common <- intersect(rownames(sample_covariates), colnames(dds))
    sample_covariates <- if (length(common) >= 3) sample_covariates[common, , drop = FALSE] else NULL
  }

  # Check if plot_volcano has been defined in the environment
  if (exists("plot_volcano", envir = .GlobalEnv)) {
    # if yes, plot_volcano_report is set to TRUE
    plot_volcano_report <- TRUE
  } else {
    # if no, plot_volcano_report is set to FALSE
    plot_volcano_report <- FALSE
  }

  if("pora_kegg_up" %in% names(results)){
    pora_kegg_up <- results$pora_kegg_up
  }
  if("pora_kegg_dn" %in% names(results)){
    pora_kegg_dn <- results$pora_kegg_dn
  }
  if("pora_custom_up" %in% names(results)){
    pora_custom_up <- results$pora_custom_up
  }
  if("pora_custom_dn" %in% names(results)){
    pora_custom_dn <- results$pora_custom_dn
  }

  if(!is.null(report.dir)){
    wd <- paste0(getwd(), "/", report.dir)
  } else {
    wd <- paste(getwd(), "/Report.rna_", format(Sys.time(),
                                                 "%Y%m%d_%H%M%S"), sep = "")
  }
  dir.create(wd, showWarnings = F)
  message("Save RData object")
  save(results, file = paste(wd, "results.RData", sep = "/"))
  # Prefer local repo template if available; else use installed package copy
  local_rmd <- "/Users/ncw/VisomX/inst/Report_RNA.Rmd"
  if (file.exists(local_rmd)) {
    file <- local_rmd
  } else {
    file <- system.file("Report_RNA.Rmd", package = "VisomX")
    if (!nzchar(file) || !file.exists(file)) {
      stop("Report_RNA.Rmd not found in local repo or installed package.")
    }
  }
  #file <- "/Users/ncw/VisomX/inst/Report_RNA.Rmd"
  message("Render reports...")
  rmarkdown::render(file, output_format = "all", output_dir = wd,
                    intermediates_dir = tempdir(), quiet = TRUE,
                      params = list(param = param),
                    envir = environment())
  message(paste0("Files saved in: '", wd, "'"))
}


#' Impute missing values in a SummarizedExperiment object
#'
#' This function imputes missing values in a SummarizedExperiment object,
#' using the specified imputation method.
#'
#' @param se A SummarizedExperiment object.
#' @param fun A character vector specifying the imputation method.
#' Options include "zero", "bpca", "knn", "QRILC", "MLE", "MinDet",
#' "MinProb", "man", "min", "zero", "mixed", "nbavg", and "SampMin".
#' @param ... Additional arguments passed to the imputation method.
#'
#' @return A SummarizedExperiment object with missing values imputed.
#' @export
#'
#' @details "SampMin" replaces missing values with the minimum value found in each sample. For information about the remaining imputation methods, see \code{help("imputeMethods", "MsCoreUtils")}
#'
#'
#' @importFrom assertthat assert_that
#' @importFrom MSnbase impute
#' @importFrom methods as
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment assay
#'
rna.impute <- function (se, fun = c("zero", "bpca", "knn", "QRILC",
                                     "MLE", "MinDet", "MinProb", "man",
                                     "min", "zero", "mixed", "nbavg", "SampMin"), ...)
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(fun))
  fun <- match.arg(fun)
  if (any(!c("name", "ID") %in% colnames(SummarizedExperiment::rowData(se,
                                                                       use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)), "'\nRun prot.make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if (!any(is.na(SummarizedExperiment::assay(se)))) {
    warning("No missing values in '", deparse(substitute(se)),
            "'. ", "Returning the unchanged object.",
            call. = FALSE)
    return(se)
  }
  SummarizedExperiment::rowData(se)$imputed <- apply(is.na(SummarizedExperiment::assay(se)), 1, any)
  SummarizedExperiment::rowData(se)$num_NAs <- rowSums(is.na(SummarizedExperiment::assay(se)))
  if(fun != "zero" && length(fun) == 1){
      if (fun == "man") {
        se <- manual_impute(se, ...)
      } else if (fun == "SampMin"){
        SummarizedExperiment::assay(se) <- SummarizedExperiment::assay(se) %>% data.frame() %>%
          mutate_if(is.numeric, function(x) ifelse(is.na(x), min(x, na.rm = T), x)) %>%
          as.matrix()
        cat("Missing values were replaced with the minimum value in each sample.")
      } else {
        MSnSet_data <- methods::as(se, "MSnSet")
        MSnSet_imputed <- MSnbase::impute(MSnSet_data, method = fun,
                                          ...)
        SummarizedExperiment::assay(se) <- MSnbase::exprs(MSnSet_imputed)
      }
  } else {
    SummarizedExperiment::assay(se) <- SummarizedExperiment::assay(se) %>% data.frame() %>%
      mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x)) %>%
      as.matrix()
    cat("Missing values were replaced with 0.\n\n")
  }
  se@metadata["imp_fun"] <- fun
  return(se)
}

#' Create a SummarizedExperiment object
#'
#' This function takes gene expression data stored in a data frame with unique gene names
#' and columns of expression data, as well as an experimental design data frame, and
#' creates a SummarizedExperiment object.
#'
#' @param genes_unique A data frame with unique gene names and columns of expression data.
#' @param columns A vector of the column positions or names with expression data in \code{genes_unique}.
#' @param expdesign A data frame of experimental design.
#'
#' @return A SummarizedExperiment object.
#'
#' @export
#' @details used internally by \code{rna.read_data()}.
rna.make_se <- function(genes_unique, columns, expdesign)
{
  assertthat::assert_that(is.data.frame(genes_unique), is.integer(columns),
                          is.data.frame(expdesign))
  if (any(!c("name", "ID") %in% colnames(genes_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)), "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns",
         "are not present in the experimental design", call. = FALSE)
  }
  if (any(!apply(genes_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric", "\nRun make_se_parse() with the appropriate columns as argument",
         call. = FALSE)
  }
  if (tibble::is_tibble(genes_unique))
    genes_unique <- as.data.frame(genes_unique)
  if (tibble::is_tibble(expdesign))
    expdesign <- as.data.frame(expdesign)
  genes_unique$name <- make.unique(genes_unique$name)
  rownames(genes_unique) <- genes_unique$name
  raw <- genes_unique[, columns]
  raw[is.na(raw)] <- 0
  expdesign <- mutate(expdesign, condition = make.names(condition)) %>%
    tidyr::unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID
  matched <- match(make.names(delete_prefix(expdesign$label)),
                   make.names(delete_prefix(colnames(raw))))
  if (any(is.na(matched))) {
    stop("None of the labels in the experimental design match ",
         "with column names in 'genes_unique'", "\nRun make_se() with the correct labels in the experimental design",
         "and/or correct columns specification")
  }
  colnames(raw)[matched] <- expdesign$ID
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
  row_data <- genes_unique[, -columns]
  rownames(row_data) <- row_data$name
  mat <- as.matrix(sapply(raw, as.integer))
  rownames(mat) <- rownames(raw)
  mat[mat==0] <- NA
  se <- SummarizedExperiment::SummarizedExperiment(assays = mat, colData = expdesign,
                                                   rowData = row_data)
  return(se)
}
