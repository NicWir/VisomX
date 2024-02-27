#' Read Transcriptomics Data
#'
#' Reads transcriptomics data from one/several file(s) or dataframe and creates a \code{SummarizedExperiment} object.
#'
#' @param data File or dataframe containing transcriptomics data, if \code{files.ind} is not used.
#' @param files.ind Prefixes of several files in the working directory containing transcriptomics data, if \code{data} is not used.
#' @param expdesign Experimental design as file path or data frame, if made previously.
#' @param csvsep Delimiter if reading CSV file(s).
#' @param dec Decimal separator if reading CSV, TSV, or TXT files.
#' @param sheet Sheet of an Excel file to be read.
#' @param name Header of column containing primary gene IDs (e.g., gene names).
#' @param id Header of column containing alternative gene IDs (e.g., locus IDs).
#' @param values Name of the column containing the values (if \code{files.ind != NULL}).
#' @param id2name.table File containing a table with ID to name mappings.
#' @param id2name.id Header of column containing alternative gene IDs in \code{id2name.table}.
#' @param id2name.name Header of column containing primary gene IDs in \code{id2name.table}.
#' @param pfx.counts Prefix in headers of columns containing gene abundances (if \code{data != NULL}).
#' @param rsd_thresh Provide a relative standard deviation (RSD) threshold **in %** for genes. The RSD is calculated for each condition and if the maximum RSD value determined for a given protein exceeds \code{rsd_thresh}, the gene is discarded. The RSD filter is applied **before** further missing value filters based on the three \code{filt_} arguments.
#' @param filt_type (Character string) "complete", "condition" or "fraction", Sets the type of filtering applied. "complete" will only keep genes with valid values in all samples. "condition" will keep genes that have a maximum of \code{filt_thr} missing values in at least one condition. "fraction" will keep genes that have a \code{filt_min} fraction of valid values in all samples.
#' @param filt_thr (Integer) Sets the threshold for the allowed number of missing values in at least one condition if \code{filt_type = "condition"}. In other words: "keep genes that have a maximum of 'filt_thr' missing values in at least one condition."
#' @param filt_min (Numeric) Sets the threshold for the minimum fraction of valid values allowed for any genes if \code{filt_type = "fraction"}.
#'
#' @return Returns a \code{SummarizedExperiment} object.
#'
#' @export
#'
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
  # Read data
  if(!is.null(data)){
    if (is.character(data)) {
    # Read table file
    dat <- read_file(data, csvsep = csvsep, dec = dec, sheet = sheet)
    } #if (is.character(data))
    else if(exists(paste(quote(data)))){
      dat <- data
    }
  } else {
    # Combine separate data files into a single dataframe
    if(!is.null(files.ind)){
      if(!is.character(files.ind)){
        stop(
          "files.ind needs to be a string or vector of several strings."
        )
      } else {
        pfx.counts <- values
        # Read data files and collect them into a list
        dfs <- list()
        files <- list.files()[grep(files.ind, list.files())]
        for(i in 1:length(files)){
          dfs[[i]] <- read_file(files[i], csvsep = csvsep)
          names(dfs)[i] <- gsub(files.ind, "", files[i])
        }
        #join dataframes in list
        df_joined <- data.frame(id = NA,
                                name = NA,
                                "Values" = NA)
        names(df_joined)[1] <- paste0(id)
        names(df_joined)[2] <- paste0(name)
        names(df_joined)[3] <- paste0(values)
        for(i in 1:length(dfs)){
          df_joined <- merge(df_joined, dfs[[i]], by = c(names(df_joined)[1], names(df_joined)[2]), all.x = TRUE, all.y = TRUE)
          colnames(df_joined) <- gsub(pattern = paste0(values, ".y"), paste0(pfx.counts,names(dfs)[i]), colnames(df_joined))
          colnames(df_joined) <- gsub(pattern = paste0(values, "$"), paste0(pfx.counts,names(dfs)[i]), colnames(df_joined))
          df_joined <- df_joined %>% select(grep(paste(values,id,name, sep = "|"), colnames(df_joined)))
        }
        dat <- df_joined[, -grep(paste0(values,".x"), colnames(df_joined))]
      } # else of if(!is.character(files.ind))
    } else { # of if(!is.null(files.ind))
      stop(
        "Please run 'rna.read_data' with an R dataframe object or a table file containing read data for all samples in the 'data' argument, or prefixes of several files in the 'file.pfx' argument."
      )
    }
  }
  # Check if column headers exist. If not, take first row as column headers
  if (all(grepl("^V\\d", colnames(dat)))){
    new_colnames <- gsub('["\']', '', dat[1,])
    colnames(dat) <- new_colnames
    dat <- dat[-1,]
  }
  if (is.character(id2name.table)) {
    # Read table file
    id2name <- read_file(id2name.table, csvsep=csvsep)
  }

  # Test for occurence of prefix for abundance columns.
  if (!(any(grepl(pfx.counts, colnames(dat))))) {
    stop(paste0("The prefix '", pfx.counts, "' does not exist in any column of '",
                data, "'. Please provide a valid prefix to identify columns with gene abundances."), call. = F)
  }

  # Make unique names using the annotations in the as name and id defined columns as primary and
  # secondary names, respectively.
  if (any(colnames(dat) %in% name) && any(colnames(dat) %in% id)) {
    # Remove rows that have 'NA' in both specified name and id columns
    row.rm <- which(is.na(dat[name]) & is.na(dat[id]), arr.ind = T)
    if (length(row.rm) > 1) {
      dat.rm <- dat[-c(row.rm),]
      message(
        paste0(
          "Removing genes with no defined name in column \"",
          name,
          "\" or \"",
          id,
          "\". (Removed: ",
          (nrow(dat) - nrow(dat.rm)),
          " genes)"
        )
      )
    }


    if (length(row.rm) > 1) {
      rna_unique <- prot.make_unique(proteins = dat.rm, names = name, ids = id, delim = ";")
    } else {
      rna_unique <- prot.make_unique(proteins = dat, names = name, ids = id, delim = ";")
    }

  } else if (!any(colnames(dat) %in% name) && !any(colnames(dat) %in% id)) {
    stop("\"", name, "\" and \"", id, "\" are not columns in \"", data, "\".", "Please provide valid column names with gene identifiers as \"name= \" and \"id= \".",
         call. = F)
  } else if (!any(colnames(dat) %in% name)) {
    stop("\"", name, "\" is not a column in \"", data, "\".", "Please provide a valid column name with gene names as \"name= \".",
         call. = F)
  } else if (!any(colnames(dat) %in% id)) {
    stop("\"", id, "\" is not a column in \"", data, "\".", "Please provide a valid column name with gene identifiers as \"id= \".",
         call. = F)
  }


  if (is.null(expdesign)) {
    # Create experimental design based on column names if neither file nor data frame is provided
    label <- rna_unique %>%
      select(., contains(pfx.counts)) %>%
      colnames() %>%
      gsub(pfx.counts, "", .)

    condition <- rna_unique %>%
      select(., contains(pfx.counts)) %>%
      colnames() %>%
      gsub(pfx.counts, "", .) %>%
      gsub(".[[:digit:]]+$", "", .)  # Remove prefix and replicate number from sample name

    replicate <- rna_unique %>%
      select(., contains(pfx.counts)) %>%
      colnames() %>%
      str_extract(., "[:digit:]{1,}$")  # Remove string before replicate number

    experimental_design <- data.frame(label, condition, replicate, check.names = FALSE)
    message(paste0(
      "Writing experimental design file to: " ,
      getwd(),
      "/experimental_design.txt"
    ))
    utils::write.table(
      experimental_design,
      file = "experimental_design.txt",
      sep = "\t",
      row.names = F
      #col.names = NA
    )

  } else if (is.list(expdesign)) {
    experimental_design <- expdesign
  } else if (is.character(expdesign) && file.exists(expdesign)) {
    if (str_replace_all(expdesign, ".{1,}\\.", "") == "csv") {
      experimental_design <- utils::read.csv( expdesign, sep = csvsep, header = T,
                                              stringsAsFactors = F, fill = T, na.strings = "",
                                              quote = "", comment.char = "", check.names=FALSE )
    } else if (str_replace_all(expdesign, ".{1,}\\.", "") == "xls" |
               str_replace_all(expdesign, ".{1,}\\.", "") == "xlsx") {
      experimental_design <- readxl::read_excel(expdesign)
    } else if (str_replace_all(expdesign, ".{1,}\\.", "") == "tsv") {
      experimental_design <- utils::read.csv(expdesign, sep = "\t", header = T,
                                             stringsAsFactors = F, fill = T, na.strings = "",
                                             quote = "", comment.char = "", check.names=FALSE )
    } else if (str_replace_all(expdesign, ".{1,}\\.", "") == "txt") {
      experimental_design <- utils::read.table(expdesign, sep = "\t", header = T, stringsAsFactors = F,
                                               fill = T, na.strings = "", comment.char = "", check.names=FALSE )
    } else {
      stop(
        "No compatible file format for experimental design provided.
             Supported formats are: \\.txt (tab delimited), \\.csv (delimiters can be specified with the argument \"csvsep = \", \\.tsv, \\.xls, and \\.xlsx"
      )
    }
  } else {
    stop(paste0("File \"", expdesign, "\" does not exist."), call. = F)
  }

  # Generate a SummarizedExperiment object using an experimental design
  abundance_columns <- grep(pfx.counts, colnames(rna_unique))  # get abundance column numbers
  rna_unique[,abundance_columns] <- suppressWarnings(matrix(as.numeric(as.matrix(rna_unique[,abundance_columns])),
                                           ncol = ncol(as.matrix(rna_unique[,abundance_columns])) ) )
  message("Generating SummarizedExperiment.")
  rna_se <- rna.make_se(rna_unique, abundance_columns, experimental_design)

  # Apply RSD threshold
  if(!is.null(rsd_thresh)){
    rsd_thresh <- rsd_thresh/100

    int.mat <- assay(rna_se)
    # replace NA with 0
    int.mat[is.na(int.mat)] <- 0

    # Calculate the standard deviation for each compound in each group. Create a list of n vectors for n sample groups.
    ls.sd <- lapply(1:length(unique(rna_se$condition)), function (x)
      apply(int.mat[,unique(rna_se$condition)[x] == gsub("_[[:digit:]]+$", "", colnames(int.mat))], 1, stats::sd, na.rm = F))
    # Calculate the average for each compound in each group. Create a list of n vectors for n sample groups.
    ls.mns <-
      lapply(1:length(unique(rna_se$condition)), function (x)
        apply(int.mat[,unique(rna_se$condition)[x] == gsub("_[[:digit:]]+$", "", colnames(int.mat))], 1, mean, na.rm = F))
    # Calculate the RSD for each compound in each group. Create a list of n vectors for n sample groups.
    ls.rsd <- mapply("/", ls.sd, ls.mns, SIMPLIFY = FALSE)
    # Combine group vectors in RSD list into matrix
    mat.rsd <- do.call(rbind, ls.rsd)
    # Create vector with the maximum RSD value among all groups for each compound.
    filter.val <- apply(mat.rsd, 2, max, na.rm = T)
    # Check if each compound has an RSD value below the threshold.
    keep <- filter.val <= rsd_thresh

    # Apply RSD filter to data matrix
    filtered_rsd <- rna_se[keep, ]
  } else {
    filtered_rsd <- rna_se
  }
  if( (nrow(SummarizedExperiment::assay(rna_se)) - nrow(SummarizedExperiment::assay(filtered_rsd))) != 0 ){
    number_removed_rsd <- nrow(SummarizedExperiment::assay(rna_se)) - nrow(SummarizedExperiment::assay(filtered_rsd))
    cat(paste0(number_removed_rsd, " out of ",
               nrow(SummarizedExperiment::assay(rna_se)), " genes were removed from the dataset due to too high RSD.\n\n"))
    filtered_rsd@metadata$n.filtered.rsd <- number_removed_rsd
  }
  filtered_rsd@metadata$rsd.thresh <- rsd_thresh
  ## Drop genes with missing values based on defined type filter "filt_thr"
  if(!is.null(filt_type)) rna_se <- rna.filter_missing(filtered_rsd, type = filt_type, thr = filt_thr, min = filt_min)

  cat(paste0("Identified conditions:\n ", paste(str_c(unique(rna_se$condition), collapse = ", ")), "\n"))

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
             "\nRun filter() with a threshold ranging from 0 to ",
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
      bin_data <- SummarizedExperiment::assay(se)
      idx <- is.na(SummarizedExperiment::assay(se))
      bin_data[!idx] <- 1
      bin_data[idx] <- 0
      keep <- bin_data %>% as.data.frame() %>% tibble::rownames_to_column() %>%
        gather(ID, value, -rowname) %>% group_by(rowname) %>%
        dplyr::summarize(n = n(), valid = sum(value), frac = valid/n) %>%
        filter(frac >= min)
      filtered <- se[keep$rowname, ]
    }
  } else {
    filtered <- se
  }
  if( (nrow(SummarizedExperiment::assay(se)) - nrow(SummarizedExperiment::assay(filtered))) != 0 ){
    number_removed <- nrow(SummarizedExperiment::assay(se)) - nrow(SummarizedExperiment::assay(filtered))
    cat(paste0(number_removed, " out of ",
               nrow(SummarizedExperiment::assay(se)), " genes were removed from the dataset due to missing values.\n\n"))
    filtered@metadata$n.filtered <- number_removed
  }
  filtered@metadata$filt_type <- type
  filtered@metadata$n.pre_filt <- nrow(SummarizedExperiment::assay(se))
  return(filtered)
}

#' RNA sequencing workflow
#'
#' This function performs a complete RNA sequencing workflow, including imputation of missing values, normalization,
#' principal component analysis, differential expression analysis, and pathway analysis. The function also provides
#' several options for plotting, exporting plots, and creating a report.
#'
#' @param se A SummarizedExperiment object, generated with read_prot().
#' @param imp_fun (Character string)  Function used for data imputation. "SampMin", "man", "bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb", "min", "zero", "mixed", or "nbavg". See (\code{\link{rna.impute}}) for details.
#' @param q (Numeric) q value for imputing missing values with method \code{imp_fun = 'MinProb'}.
#' @param knn.rowmax (Numeric) The maximum percent missing data allowed in any row for \code{imp_fun = 'knn'}. Default: 0.5.
#' @param type (Character string) Type of differential analysis to perform. "all" (contrast each condition with every other condition), "control" (contrast each condition to a defined control condition), "manual" (manually define selected conditions).
#' @param design Formula for the design matrix.
#' @param size.factors Optional: Manually define size factors for normalization.
#' @param altHypothesis  Specify those genes you are interested in finding. The test provides p values for the null hypothesis, the complement of the set defined by altHypothesis. For further details, see \code{\link[DESeq2]{results}}.
#' @param control Control condition; required if type = "control".
#' @param contrast (String or vector of strings) Defined test(s) for differential analysis in the form "A_vs_B"; required if type = "manual".
#' @param controlGenes Specifying those genes to use for size factor estimation (e.g. housekeeping or spike-in genes).
#' @param pAdjustMethod Method for adjusting p values. Available options are "IHW","BH".
#' @param alpha Significance threshold for adjusted p values.
#' @param alpha.independent Adjusted p value threshold for independent filtering or NULL. If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
#' @param alpha_pathways Significance threshold for pathway analysis.
#' @param lfcShrink Use shrinkage to calculate log2 fold change values.
#' @param shrink.method Method for shrinkage. Available options are "apeglm", "ashr", "normal". See \code{\link[DESeq2]{lfcShrink}} for details.
#' @param lfc Relevance threshold for absolute log2(fold change) values. Used to filter unshrunken lfc values
#' or in shrinkage method "apeglm" or "normal".
#' @param heatmap.show_all Shall all samples be displayed in the heatmap or only the samples contained in the defined
#' "contrast"? (only applicable for type = "manual")
#' @param heatmap.kmeans Shall the proteins be clustered in the heat map?
#' @param k Number of protein clusters in heat map if kmeans = TRUE.
#' @param heatmap.col_limit Define the breaks in the heat map legends.
#' @param heatmap.show_row_names Show protein names in heat map?
#' @param heatmap.row_font_size Font size of protein names if show_row_names = TRUE.
#' @param volcano.add_names Display names next to symbols in volcano plot.
#' @param volcano.label_size Size of labels in volcano plot.
#' @param volcano.adjusted Display adjusted p-values on y axis of volcano plot?
#' @param plot Shall plots be returned in the Plots pane?
#' @param export Shall plots be exported as PDF and PNG files?
#' @param report Shall a report (HTML and PDF) be created?
#' @param report.dir Folder name for created report (if report = TRUE)
#' @param pathway_enrichment Perform pathway over-representation analysis for each tested contrast
#' @param pathway_kegg Perform pathway over-representation analysis with gene sets in the KEGG database
#' @param kegg_organism Name of the organism in the KEGG database (if 'pathway_kegg = TRUE')
#' @param custom_pathways Dataframe providing custom pathway annotations.
#' @param quiet Suppress messages and warnings.
#'
#' @return The function returns a SummarizedExperiment object with added columns for log2 fold change, p-values and adjusted p-values for each comparison.
#' It also includes a column for significant genes for each comparison and a column for significant genes overall.
#' Additionally, the function generates various plots and a report (if specified).

#' @export
rna.workflow <- function(se, # SummarizedExperiment, generated with read_prot().
                          imp_fun = c("zero", "man", "bpca", "knn", "QRILC", "MLE", "MinDet",
                                      "MinProb", "min", "zero", "mixed", "nbavg", "SampMin"),
                          q = 0.01,
                          knn.rowmax = 0.5,
                          type = c("all", "control", "manual"),
                          design = "~ condition",
                          size.factors = NULL,
                          altHypothesis = c("greaterAbs", "lessAbs", "greater", "less"),
                          control = NULL,
                          contrast = NULL,
                          controlGenes = NULL,
                          pAdjustMethod = c("IHW","BH"),
                          alpha = 0.05,
                          alpha.independent = 0.1,
                          alpha_pathways = 0.1,
                          lfcShrink = TRUE,
                          shrink.method = c("apeglm", "ashr", "normal"),
                          lfc = 2,
                          heatmap.show_all = TRUE,
                          heatmap.kmeans = F,
                          k = 6,
                          heatmap.col_limit = NA,
                          heatmap.show_row_names = TRUE,
                          heatmap.row_font_size = 6,
                          volcano.add_names = FALSE,
                          volcano.label_size = 2.5,
                          volcano.adjusted = TRUE,
                          plot = FALSE,
                          export = FALSE,
                          report = TRUE,
                          report.dir = NULL,
                          pathway_enrichment = FALSE,
                          pathway_kegg = FALSE,
                          kegg_organism = NULL,
                          custom_pathways = NULL,
                          quiet = FALSE
)
{
  # Show error if inputs are not the required classes
  imp_fun <- match.arg(imp_fun)
  shrink.method <- match.arg(shrink.method)
  pAdjustMethod <- match.arg(pAdjustMethod)
  if(is.integer(alpha)) alpha <- as.numeric(alpha)
  if(is.integer(lfc)) lfc <- as.numeric(lfc)
  assertthat::assert_that(is.character(imp_fun),
                          is.character(type),
                          is.numeric(alpha),
                          length(alpha) == 1,
                          is.numeric(lfc),
                          length(lfc) == 1
                          )

  # Show error if inputs are not valid
  if (!type %in% c("all", "control", "manual")) {
    stop("run workflow_proteomics() with a valid type",
         "\nValid types are: 'all', 'control' and 'manual'.",
         call. = FALSE)
  }
  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
  }
  # Impute missing values
    if (imp_fun == "MinProb"){
        se_imp <- rna.impute(se, fun = imp_fun, q = q)
    } else if (imp_fun == "knn"){
        se_imp <- rna.impute(se, fun = imp_fun, rowmax = knn.rowmax)
    } else {
        se_imp <- rna.impute(se, fun = imp_fun)
    }

  # Create DESeqDataSet and define contrasts
  if(!quiet) message(">> Creating DESeqDataSet object <<")
  conditions <- as.character(unique(se$condition))
  se_imp$condition <- as.factor(se_imp$condition)
  if(!quiet){
    ddsSE <- DESeq2::DESeqDataSet(se_imp, design = as.formula(design))
  } else {
    ddsSE <- suppressMessages( DESeq2::DESeqDataSet(se_imp, design = as.formula(design)) )
  }
  if (type == "control") {
    if (is.null(control))
      stop("run test_diff(type = 'control') with a 'control' argument")
    cntrst <- paste(conditions[!conditions %in% control],
                    control, sep = "_vs_")
    ddsSE$condition <- relevel(ddsSE$condition, control)
  }
  if (type == "manual") {
    if (is.null(contrast)) {
      stop("run test_diff(type = 'manual') with a 'test' argument")
    }
    assertthat::assert_that(is.character(contrast))
    if (any(!unlist(strsplit(contrast, "_vs_")) %in% conditions)) {
      stop("run rna.workflow() with valid contrasts in 'contrast'",
           ".\nValid contrasts should contain combinations of: '",
           paste0(conditions, collapse = "', '"),
           "', for example '", paste0(conditions[1],
                                      "_vs_", conditions[2]), "'.", call. = FALSE)
    }
    ddsSE$condition <- relevel(ddsSE$condition, unlist(strsplit(contrast, "_vs_"))[2])
    cntrst <- contrast
  }
  if (type == "all") {
    cntrst <- apply(utils::combn(conditions, 2), 2, paste,
                    collapse = "_vs_")
    if (!is.null(control)) {
      flip <- grep(paste("^", control, "_vs_", sep = ""),
                   cntrst)
      if (length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>% gsub(paste(control,
                                                    "_vs_", sep = ""), "", .) %>%
          paste("_vs_", control, sep = "")
      }
    }
  }


  # Pre-filtering:  keep only rows that have at least 10 reads total
  keep <- BiocGenerics::rowSums(BiocGenerics::counts(ddsSE) >= 5) >= 3
  ddsSE <- ddsSE[keep,]

  # Differential expression analysis; Performs 'replaceOutliers' (replaces outlier counts flagged by extreme Cook's distances)
  dds <- DESeq2::DESeq(ddsSE, sfType = "ratio", quiet = quiet)

  # Perform PCA Analysis
  if(!quiet) message("performing PCA analysis")
  norm.counts <- BiocGenerics::counts(dds, normalized = TRUE)
  rlog.counts <- tryCatch(DESeq2::rlog(dds, fitType = 'mean'), error = function(e) { rlog(dds, fitType = 'mean') })
  rna.pca <- prot.pca(SummarizedExperiment::assay(rlog.counts))

  # Create list with test results for defined contrasts; Perform independent hypothesis weighting if pAdjustMethod == "IHW"
  if(!quiet) message("assembling results for defined contrasts")
  results <- list()
  for(i in 1:length(cntrst)){
    trmt <- str_split(cntrst[i], "_vs_")[[1]][1]
    ctrl <- str_split(cntrst[i], "_vs_")[[1]][2]
    if(!is.null(altHypothesis)){
      results[[i]] <- DESeq2::results(dds, contrast=c("condition", trmt, ctrl), independentFiltering = ifelse(!is.null(alpha.independent), TRUE, FALSE),
                                      altHypothesis = altHypothesis, lfcThreshold = lfc,
                                      alpha = alpha.independent, filterFun = ifelse(pAdjustMethod == "IHW", IHW::ihw, NULL))
    } else {
      results[[i]] <- DESeq2::results(dds, contrast=c("condition", trmt, ctrl), independentFiltering = ifelse(!is.null(alpha.independent), TRUE, FALSE),
                                      alpha = alpha.independent, filterFun = ifelse(pAdjustMethod == "IHW", IHW::ihw, NULL))
    }

  }
  names(results) <- cntrst

  # Store condition levels in object
  levels <- levels(ddsSE$condition)

  # Creation of shrunken results list after re-running nbinomWaldTest with re-ordered conditions
  if(lfcShrink == TRUE){
    if(!quiet) message(paste0("performing lfc shrinkage with method '", shrink.method, "'"))
    results_shrink <- list()
    if(shrink.method == "apeglm"){
      cntrst_ndx <- lapply(2:length(DESeq2::resultsNames(dds)), function(x) match(gsub("condition_", "", DESeq2::resultsNames(dds)[x]), cntrst))
      cntrst_ndx <- cntrst_ndx[!sapply(cntrst_ndx,is.na)]
      names(cntrst_ndx) <- cntrst[unlist(cntrst_ndx)]
      if(length(cntrst_ndx) > 0){
        for(i in 1:length(cntrst_ndx)){
          results_shrink[[i]] <- DESeq2::lfcShrink(dds, coef = match(names(cntrst_ndx)[i], gsub("condition_", "", DESeq2::resultsNames(dds))),
                                                   type = "apeglm", returnList = F, quiet = T, )
        }
      }
      names(results_shrink) <- names(cntrst_ndx)
      if (type == "all") {
        for(i in 2:length(levels)){
          dds$condition <- relevel(dds$condition, levels[i])
          dds <- DESeq2::nbinomWaldTest(object = dds, quiet = T)
          cntrst_ndx <- lapply(2:length(DESeq2::resultsNames(dds)), function(x) match(gsub("condition_", "", DESeq2::resultsNames(dds)[x]), cntrst))
          cntrst_ndx <- cntrst_ndx[!sapply(cntrst_ndx,is.na)]
          names(cntrst_ndx) <- cntrst[unlist(cntrst_ndx)]
          if(length(cntrst_ndx) > 0){
            for(j in 1:length(cntrst_ndx)){
              results_shrink[[length(results_shrink)+1]] <- DESeq2::lfcShrink(dds, coef = match(names(cntrst_ndx)[j], gsub("condition_", "", DESeq2::resultsNames(dds))),
                                                                              type = "apeglm", returnList = F, quiet = T)
              names(results_shrink)[length(results_shrink)] <- names(cntrst_ndx)[j]
            }
          }
        }
      }
    }
    if(shrink.method == "ashr" || shrink.method == "normal"){
      for (i in 1:length(results)) {
        trmt <- str_split(cntrst[i], "_vs_")[[1]][1]
        ctrl <- str_split(cntrst[i], "_vs_")[[1]][2]
        results_shrink[[i]] <- DESeq2::lfcShrink(dds, contrast = c("condition", trmt, ctrl), type = shrink.method,
                                                 lfcThreshold = ifelse(shrink.method == "normal", lfc, 0), quiet = T)
      }
      names(results_shrink) <- cntrst
    }
    # add shrunken lfc values to results object
    for(i in 1:length(cntrst)){
      results[[match(cntrst[i], names(results))]]$lfc.shrink <- as.vector(results_shrink[[match(cntrst[i], names(results_shrink))]]$log2FoldChange)
      results[[match(cntrst[i], names(results))]]$lfcSE.shrink <- as.vector(results_shrink[[match(cntrst[i], names(results_shrink))]]$lfcSE)
    }
  } # if(lfcShrink == TRUE)

  # add significant column for each contrast
  for(i in 1:length(cntrst)){
    if(lfcShrink == TRUE){
      if(!is.null(altHypothesis)){
        results[[match(cntrst[i], names(results))]]$significant <-
          results[[match(cntrst[i], names(results))]]$padj <= alpha
      } else {
        results[[match(cntrst[i], names(results))]]$significant <-
          abs(results[[match(cntrst[i], names(results))]]$lfc.shrink) >= lfc &
          results[[match(cntrst[i], names(results))]]$padj <= alpha
      }
    }
    else{
      if(!is.null(altHypothesis)){
        results[[match(cntrst[i], names(results))]]$significant <-
          results[[match(cntrst[i], names(results))]]$padj <= alpha
      } else {
        results[[match(cntrst[i], names(results))]]$significant <-
          abs(results[[match(cntrst[i], names(results))]]$log2FoldChange) >= lfc &
          results[[match(cntrst[i], names(results))]]$padj <= alpha
      }
    }
  }

  if(!quiet) message("gathering results into list")
  # store results values for each contrast in dds object
  for(i in 1:length(results)){
    SummarizedExperiment::rowData(dds)[,paste0("lfc.", names(results)[i])] <- results[[i]]$log2FoldChange
    SummarizedExperiment::rowData(dds)[,paste0("lfcSE.", names(results)[i])] <- results[[i]]$lfcSE
    if(lfcShrink == TRUE){
      SummarizedExperiment::rowData(dds)[,paste0("lfc_shrink.", names(results)[i])] <- results[[i]]$lfc.shrink
      SummarizedExperiment::rowData(dds)[,paste0("lfcSE_shrink.", names(results)[i])] <- results[[i]]$lfcSE.shrink
    }
    SummarizedExperiment::rowData(dds)[,paste0("pvalue.", names(results)[i])] <- results[[i]]$pvalue
    SummarizedExperiment::rowData(dds)[,paste0("padj.", names(results)[i])] <- results[[i]]$padj
    SummarizedExperiment::rowData(dds)[,paste0("significant.", names(results)[i])] <- results[[i]]$significant
  }
  # Add 'significant' column if gene was significant in any contrast
  if(length(grep("significant.", colnames(SummarizedExperiment::rowData(dds)))) > 1){
    SummarizedExperiment::rowData(dds)[,"significant"] <- apply(SummarizedExperiment::rowData(dds)[,grep("significant.", colnames(SummarizedExperiment::rowData(dds)))], 1, any)
  }
  else{
    SummarizedExperiment::rowData(dds)[,"significant"] <- SummarizedExperiment::rowData(dds)[,grep("significant.", colnames(SummarizedExperiment::rowData(dds)))]
  }
  if(!quiet) message("Significantly different genes identified:")
  if(!quiet){
    for (i in 1:length(cntrst)) {
      message(paste0("* ",
                     nrow(dds[SummarizedExperiment::rowData(dds)[[paste0("significant.", cntrst[i])]][!is.na(SummarizedExperiment::rowData(dds)[[paste0("significant.", cntrst[i])]])],]),
                     " genes for contrast ", cntrst[i], "\n"))
    }
  }

  # Perform pathway enrichment analysis
  if (pathway_enrichment) {
    if(!quiet) message("performing pathway enrichment analysis")
    res.pathway <- enrich_pathways(dds, contrasts = cntrst, alpha_pathways = alpha_pathways, pathway_kegg=pathway_kegg, kegg_organism, custom_pathways)
  }
  else {
    res.pathway <- NA
  }

  # Write output
  df <- data.frame(SummarizedExperiment::rowData(dds))
  res.table <- data.frame("gene_symbol" = df$ID,
                          "gene_name" = df$name,
                          "imputed" = df$imputed,
                          "baseMean" = df$baseMean,
                          "baseVar" = df$baseVar,
                          df[,grep("lfc\\.", colnames(df))[1]:match("significant", colnames(df))]
                          )
  if(!quiet) message(paste0("Save results as tab-delimited table to: ", getwd(), "/results.txt"))
  utils::write.table(res.table, paste(getwd(), "results.txt",
                                  sep = "/"), row.names = FALSE, sep = "\t")

  param <- data.frame(type, design, altHypothesis = ifelse(!is.null(altHypothesis), altHypothesis, NA), controlGenes = ifelse(is.null(controlGenes), "none", controlGenes), pAdjustMethod,
                      lfcShrink, shrink.method, alpha, lfc, check.names = FALSE)

  res <- list(data = SummarizedExperiment::rowData(se), se = se, imputed = se_imp,
              dds = dds, pca = rna.pca,
              results = results, param = param)

  if (pathway_enrichment == T && pathway_kegg) {
    res <- c(res, pora_kegg_up = list(res.pathway$ls.pora_kegg_up), pora_kegg_dn = list(res.pathway$ls.pora_kegg_dn))
    if(!quiet) message(paste0("Writing results of KEGG pathway enrichment analysis to: ", getwd(), "/pora_kegg_contrast...txt"))
    for(i in 1:length(res.pathway$ls.pora_kegg_up)){
      if(!is.null(res.pathway$ls.pora_kegg_up[[i]])){
        utils::write.table(res.pathway$ls.pora_kegg_up[[i]]@result, paste(getwd(), paste0("pora_kegg_", names(res.pathway$ls.pora_kegg_up)[i], "_up.txt"),
                                                                          sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
    for(i in 1:length(res.pathway$ls.pora_kegg_dn)){
      if(!is.null(res.pathway$ls.pora_kegg_dn[[i]])){
        utils::write.table(res.pathway$ls.pora_kegg_dn[[i]]@result, paste(getwd(), paste0("pora_kegg_", names(res.pathway$ls.pora_kegg_dn)[i], "_down.txt"),
                                                                          sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
  }
  if (pathway_enrichment == T && !is.null(custom_pathways)) {
    res <- c(res, pora_custom_up = list(res.pathway$ls.pora_custom_up), pora_custom_dn = list(res.pathway$ls.pora_custom_dn))
    if(!quiet) message(paste0("Writing results of custom pathway enrichment analysis to: ", getwd(), "/pora_custom_contrast...txt"))
    for(i in 1:length(res.pathway$ls.pora_custom_up)){
      if(!is.null(res.pathway$ls.pora_custom_up[[i]])){
        utils::write.table(res.pathway$ls.pora_custom_up[[i]]@result, paste(getwd(), paste0("pora_custom_", names(res.pathway$ls.pora_custom_up)[i], "_up.txt"),
                                                                            sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
    for(i in 1:length(res.pathway$ls.pora_custom_dn)){
      if(!is.null(res.pathway$ls.pora_custom_dn[[i]])){
        utils::write.table(res.pathway$ls.pora_custom_dn[[i]]@result, paste(getwd(), paste0("pora_custom_", names(res.pathway$ls.pora_custom_dn)[i], "_down.txt"),
                                                                            sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
  }

  # Create Plots
  if(export == TRUE || plot == TRUE){
    if (export == TRUE){
      dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
      if(!quiet) message(paste0("Rendering and exporting figures as PDF and PNG files to: ", getwd(), "/Plots"))
    }
    suppress_ggrepel <- function(w) {
      if (any(grepl("ggrepel", w)))
        invokeRestart("muffleWarning")
    }
    #### Box plot of Cook's distance for each sample (outlier detection)
    if(plot==TRUE){
      par(mar=c(8,5,2,2))
      boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
    }
    if (export == TRUE){
      grDevices::pdf(paste0("Plots/BoxPlot_CooksDistance", ".pdf"),
                     width = 6*length(dds$condition)/20, height = 6)
      par(mar=c(8,5,2,2))
      boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
      grDevices::dev.off()

      grDevices::png(paste0("Plots/BoxPlot_CooksDistance", ".png"),
                     width = 6*length(dds$condition)/20, height = 6, units = 'in', res = 300)
      par(mar=c(8,5,2,2))
      boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
      grDevices::dev.off()
    }
    try(suppressMessages(
      prot.plot_missval(se, plot = plot, export = export)
    ))

    se_log2 <- se
    assay(se_log2) <- log2(assay(se_log2))
    try(suppressMessages(
      prot.plot_detect(se_log2, basesize = 10, plot = plot, export = export)
    ))
    suppressMessages(
      rna.plot_imputation(log2(SummarizedExperiment::assay(se)+1),
                          "Imputed" = log2(SummarizedExperiment::assay(se_imp)+1),
                          "Normalized" = log2(norm.counts+1), colData = SummarizedExperiment::colData(se),
                          plot = plot, export = export, basesize = 12)
    )
    suppressMessages(
      suppressWarnings(
        prot.plot_screeplot(rna.pca, axisLabSize = 18, titleLabSize = 22, plot = plot, export = export)
      ) )
    withCallingHandlers(suppressMessages(
      prot.plot_loadings(rna.pca,labSize = 3, plot = plot, export = export)
    ) , warning = suppress_ggrepel)
    suppressMessages(
      rna.plot_pca(dds,x = 1, y = 2, point_size = 4, basesize = 14, title = "PC Scores - PC2 vs. PC1",
                    plot = plot,export = export)
    )
    suppressMessages(
      rna.plot_pca(dds,x = 1, y = 3, point_size = 4, basesize = 14, title = "PC Scores - PC3 vs. PC1",
                    plot = plot,export = export)
    )
    suppressMessages(
      rna.plot_heatmap(dds, type = "centered", kmeans = heatmap.kmeans, show_all = heatmap.show_all, contrast = cntrst,
                        k = k, col_limit = heatmap.col_limit,show_row_names = heatmap.show_row_names,
                        row_font_size = heatmap.row_font_size, indicate = c("condition"),
                        plot = plot, export = export)
    )
    suppressMessages(
      rna.plot_heatmap(dds, type = "contrast", contrast = cntrst,  kmeans = heatmap.kmeans, k = k, col_limit = heatmap.col_limit,
                        show_row_names = heatmap.show_row_names, row_font_size = heatmap.row_font_size,
                        plot = plot, export = export)
    )
    for (i in 1:length(cntrst)){
      suppressMessages(
        suppressWarnings(
          volcano.tmp <- rna.plot_volcano(dds, contrast = cntrst[i],
                            add_names = volcano.add_names, label_size = volcano.label_size, adjusted =  volcano.adjusted,
                            plot = plot, export = export, lfc = lfc, alpha = alpha)
        ) )
    }
    if (pathway_kegg) {
      for (i in 1:length(cntrst)) {
        if(!(nrow(as.data.frame(res.pathway$ls.pora_kegg_up[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(res.pathway$ls.pora_kegg_up[[i]], title = paste0("Upregulated pathways", " - KEGG"),
                                   subtitle =  str_replace(cntrst[i], "_vs_", " vs. "), plot = plot, export = export, kegg = TRUE)
            ) )
        }
        if(!(nrow(as.data.frame(res.pathway$ls.pora_kegg_dn[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(res.pathway$ls.pora_kegg_dn[[i]], title = paste0("Downregulated pathways", " - KEGG"),
                                   subtitle = str_replace(cntrst[i], "_vs_", " vs. "), plot = plot, export = export, kegg = TRUE)
            ) )
        }
      }
    }
    if(!is.null(custom_pathways)){
      for (i in 1:length(cntrst)) {
        if(!(nrow(as.data.frame(res.pathway$ls.pora_custom_up[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(res.pathway$ls.pora_custom_up[[i]], title = "Upregulated pathways",
                                   subtitle =  str_replace(cntrst[i], "_vs_", " vs. "), plot = plot, export = export, kegg = FALSE)
            ) )
        }
        if(!(nrow(as.data.frame(res.pathway$ls.pora_custom_dn[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(res.pathway$ls.pora_custom_dn[[i]], title = "Downregulated pathways",
                                   subtitle = str_replace(cntrst[i], "_vs_", " vs. "), plot = plot, export = export, kegg = FALSE)
            ) )
        }
      }
    } # if(!is.null(custom_pathways))
  } #if(export == TRUE || plot == TRUE)

  if(report == TRUE){
    rna.report(res,
               volcano.adjusted = volcano.adjusted,
               volcano.add_names = volcano.add_names,
               pathway_enrichment = pathway_enrichment,
               heatmap.show_all = heatmap.show_all,
               heatmap.kmeans = heatmap.kmeans,
               heatmap.row_font_size = heatmap.row_font_size,
               export = export,
               k = k,
               report.dir = report.dir)
  }
  return(res)
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
rna.report <- function(results, report.dir = NULL, ...){
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
  param <- results$param
  res <- results
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
    wd <- paste(getwd(), "/Report.prot_", format(Sys.time(),
                                                 "%Y%m%d_%H%M%S"), sep = "")
  }
  dir.create(wd, showWarnings = F)
  message("Save RData object")
  save(results, file = paste(wd, "results.RData", sep = "/"))
  file <- paste("C:/Users/nicwir/Documents/DTU_Biosustain/Scripts_and_Modelling/fluctuator/220111/R_package/VisomX/inst", "/Report_RNA.Rmd",
                sep = "")
  message("Render reports...")
  rmarkdown::render(file, output_format = "all", output_dir = wd,
                    quiet = TRUE)
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
