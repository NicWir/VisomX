#' Filter proteins based on missing values
#'
#' \code{prot.filter_missing} filters a proteomics dataset based on missing values. Different types of filtering can be applied, which range from only keeping proteins without missing values to keeping proteins with a certain percent valid values in all samples or keeping proteins that are complete in at least one condition.
#'
#' @param se \code{SummarizedExperiment} object, proteomics data parsed with \code{\link{prot.read_data}}.
#' @param type (Character string) "complete", "condition" or "fraction", Sets the type of filtering applied. "complete" will only keep proteins with valid values in all samples. "condition" will keep proteins that have a maximum of \code{thr} missing values in at least one condition. "fraction" will keep proteins that have a \code{min} fraction of valid values in all samples.
#' @param thr (Integer) Sets the threshold for the allowed number of missing values in at least one condition if \code{type = "condition"}. In other words: "keep proteins that have a maximum of 'thr' missing values in at least one condition."
#' @param min (Numeric) Sets the threshold for the minimum fraction of valid values allowed for any protein if \code{type = "fraction"}.
#'
#' @return A filtered SummarizedExperiment object.
#' @export
prot.filter_missing <- function (se, type = c("complete", "condition", "fraction", NULL),
                                 thr = NULL, min = NULL)
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  if(is.character(type)) type <- match.arg(type)
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
               nrow(SummarizedExperiment::assay(se)), " proteins were removed from the dataset due to missing values.\n\n"))
    filtered@metadata$n.filtered <- number_removed
  }
  filtered@metadata$filt_type <- type
  filtered@metadata$n.pre_filt <- nrow(SummarizedExperiment::assay(se))
  return(filtered)
}

#' Read proteomics data in table format and create SummarizedExperiment
#'
#' \code{prot.read_data} takes a table file containing proteomics data and filters proteins based on missing values (see \code{\link{prot.filter_missing}}) or a given relative standard deviation threshold, and creates a \code{SummarizedExperiment} object.
#'
#' @param data An R dataframe object or a table file with extension '.xlsx', '.xls', '.csv', '.tsv', or '.txt' containing proteomics data.
#' The table must contain:
#' \enumerate{
#'    \item A column with protein IDs (e.g., accession numbers). The header of this column is provided as argument \code{id}.
#'    \item A column with protein names. The header of this column is provided as argument \code{name}.
#'    \item X columns containing abundance values for X samples. The column headers must have a prefix (e.g., "abundances.") that is provided as argument \code{pfx}. Replicates are identified by identical column headers followed by an underscore and the replicate number (e.g., "abundances.ConditionA_1", "abundances.ConditionA_2", "abundances.ConditionA_3", ...).
#' }
#' @param expdesign (_optional, if made previously_) An R dataframe object or a table file containing the columns 'label', 'condition', and 'replicate' with label = "condition_replicate". If \code{NULL}, an experimental design table will be created automatically.
#' @param csvsep (Character string) separator used in CSV files (ignored for other file types). Default: \code{";"}
#' @param dec (Character string) decimal separator used in CSV, TSV or TXT files (ignored for other file types). Default: \code{"."}
#' @param na.strings A character vector of strings which are to be interpreted as NA values.
#' @param sheet (Integer or Character string) Number or name of the sheet with proteomics data in XLS or XLSX files (_optional_).
#' @param filter (Character string or vector of strings) Provide the header of a column containing "+" or "-" to indicate if proteins should be discarded or kept, respectively.
#' @param rsd_thresh (Numeric, optional) Provide a relative standard deviation (RSD) threshold **in %** for proteins. The RSD is calculated for each condition and if the maximum RSD value determined for a given protein exceeds \code{rsd_thresh}, the protein is discarded. The RSD filter is applied **before** further missing value filters based on the three \code{filt_} arguments.
#' @param name (Character string) Provide the header of the column containing protein names.
#' @param id (Character string) Provide the header of the column containing protein IDs
#' @param pfx (Character string) Provide the common prefix for headers containing abundance values (e.g., "abundances.").
#' @param filt_type (Character string) "complete", "condition" or "fraction", Sets the type of filtering applied. "complete" will only keep proteins with valid values in all samples. "condition" will keep proteins that have a maximum of \code{filt_thr} missing values in at least one condition. "fraction" will keep proteins that have a \code{filt_min} fraction of valid values in all samples.
#' @param filt_thr (Integer) Sets the threshold for the allowed number of missing values in at least one condition if \code{filt_type = "condition"}. In other words: "keep proteins that have a maximum of 'filt_thr' missing values in at least one condition."
#' @param filt_min (Numeric) Sets the threshold for the minimum fraction of valid values allowed for any protein if \code{filt_type = "fraction"}.
#'
#' @return A filtered SummarizedExperiment object.
#' @export
#'
#' @importFrom SummarizedExperiment assay rowData colData
prot.read_data <- function (data = "dat_prot.csv",
                            expdesign = NULL,
                            csvsep = ";",
                            dec = ".",
                            na.strings = "",
                            sheet = 1,
                            filter = c("Reverse", "Potential contaminant"),
                            rsd_thresh = NULL,
                            name = 'Gene Symbol',
                            id = 'Ensembl Gene ID',
                            pfx = "abundances.",
                            filt_type = c("condition", "complete", "fraction", NULL),
                            filt_thr = 3,
                            filt_min = NULL
) {
  assertthat::assert_that(is.character(name),
                          length(name) == 1,
                          is.character(id),
                          length(id) == 1)
  if(is.character(filt_type)) filt_type <- match.arg(filt_type)
  data.object <- data
  # Read data file
  if(!is.null(data)){
    if (is.character(data)) {
      # Read table file
      prot <- read_file(data, csvsep = csvsep, dec = dec, sheet = sheet)
    } #if (is.character(data))
    else if(exists(paste(quote(data)))){
      prot <- data
    }
  } else {
    stop("Please provide proteomics data as either an R dataframe object or a CSV/TXT/XLS/XLSX/TSV table file.")
  }

  # Test for occurence of prefix for abundance columns.
  if (!(any(grepl(pfx, colnames(prot))))) {
    stop(paste0("The prefix '", pfx, "' does not exist in any column of '",
                paste(quote(data.object)), "'. Please provide a valid prefix to identify columns with protein abundances."), call. = F)
  }

  # Filter out the positive proteins (indicated by '+')
  # in the pre-defined "filter" columns
  cols_filt <- grep(paste("^", filter, "$", sep = "", collapse = "|"),  # The columns to filter on
                    colnames(prot))
  if (!is.null(cols_filt)) {
    if (length(cols_filt) == 1) {
      rows <- which(prot[, cols_filt] == "+")
      if(length(rows) > 0) prot <- prot[-rows,]
    } else {
      rows <- which(apply(prot[, cols_filt] == "+", 1, any))
      if(length(rows) > 0) prot <- prot[-rows,]
    }
  }

  # Make unique names using the annotations in the as name and id defined columns as primary and
  # secondary names, respectively.
  if (any(colnames(prot) %in% name) && any(colnames(prot) %in% id)) {
    # Remove rows that have 'NA' in both specified name and id columns
    row.rm <- which(is.na(prot[name]) & is.na(prot[id]), arr.ind = T)
    if (length(row.rm) > 1) {
      prot.rm <- prot[-c(row.rm),]
      message(
        paste0(
          "Removing proteins with no defined name in column \"",
          name,
          "\" or \"",
          id,
          "\". (Removed: ",
          (nrow(prot) - nrow(prot.rm)),
          " proteins)"
        )
      )
    }


    if (length(row.rm) > 1) {
      prot_unique <- prot.make_unique(proteins = prot.rm, names = name, ids = id, delim = ";")
    } else {
      prot_unique <- prot.make_unique(proteins = prot, names = name, ids = id, delim = ";")
    }

  } else if (!any(colnames(prot) %in% name) && !any(colnames(prot) %in% id)) {
    stop("\"", name, "\" and \"", id, "\" are not columns in \"", data, "\".", "Please provide valid column names with protein identifiers as \"name= \" and \"id= \".",
         call. = F)
  } else if (!any(colnames(prot) %in% name)) {
    stop("\"", name, "\" is not a column in \"", data, "\".", "Please provide a valid column name with protein names as \"name= \".",
         call. = F)
  } else if (!any(colnames(prot) %in% id)) {
    stop("\"", id, "\" is not a column in \"", data, "\".", "Please provide a valid column name with protein identifiers as \"id= \".",
         call. = F)
  }


  if (is.null(expdesign)) {
    # Create experimental design based on column names if neither file nor data frame is provided
    label <- prot_unique %>%
      select(., contains(pfx, ignore.case = F)) %>%
      colnames() %>%
      gsub(pfx, "", .)

    condition <- prot_unique %>%
      select(., contains(pfx, ignore.case = F)) %>%
      colnames() %>%
      gsub(pfx, "", .) %>%
      gsub(".[[:digit:]]+$", "", .)  # Remove prefix and replicate number from sample name

    replicate <- prot_unique %>%
      select(., contains(pfx, ignore.case = F)) %>%
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
  abundance_columns <- grep(pfx, colnames(prot_unique))  # get abundance column numbers
  if(any(sapply(prot_unique[,abundance_columns], function(x) any(is.character(x))))){
    prot_unique[,abundance_columns] <- apply(prot_unique[,abundance_columns], 2, function(x) as.numeric(na_if(x, "NA"))) # replace "NA" with NA
    prot_unique[, abundance_columns] <- # convert to numeric
      sapply(1:length(abundance_columns), function (x)
        as.numeric(prot_unique[, abundance_columns[x]]))
  }
  message("Generating SummarizedExperiment.")
  prot_se <- make_se(prot_unique, abundance_columns, experimental_design)

  # Apply RSD threshold
  if(!is.null(rsd_thresh)){
    rsd_thresh <- rsd_thresh/100

    int.mat <- 2^assay(prot_se)
    # replace NA with 0
    int.mat[is.na(int.mat)] <- 0

    # Calculate the standard deviation for each compound in each group. Create a list of n vectors for n sample groups.
    ls.sd <- lapply(1:length(unique(se$condition)), function (x)
      apply(int.mat[,unique(se$condition)[x] == gsub("_[[:digit:]]+$", "", colnames(int.mat))], 1, stats::sd, na.rm = F))
    # Calculate the average for each compound in each group. Create a list of n vectors for n sample groups.
    ls.mns <-
      lapply(1:length(unique(se$condition)), function (x)
        apply(int.mat[,unique(se$condition)[x] == gsub("_[[:digit:]]+$", "", colnames(int.mat))], 1, mean, na.rm = F))
    # Calculate the RSD for each compound in each group. Create a list of n vectors for n sample groups.
    ls.rsd <- mapply("/", ls.sd, ls.mns, SIMPLIFY = FALSE)
    # Combine group vectors in RSD list into matrix
    mat.rsd <- do.call(rbind, ls.rsd)
    # Create vector with the maximum RSD value among all groups for each compound.
    filter.val <- apply(mat.rsd, 2, max, na.rm = T)
    # Check if each compound has an RSD value below the threshold.
    keep <- filter.val <= rsd_thresh

    # Apply RSD filter to data matrix
    filtered_rsd <- prot_se[keep, ]
  } else {
    filtered_rsd <- prot_se
  }
  if( (nrow(SummarizedExperiment::assay(prot_se)) - nrow(SummarizedExperiment::assay(filtered_rsd))) != 0 ){
    number_removed_rsd <- nrow(SummarizedExperiment::assay(prot_se)) - nrow(SummarizedExperiment::assay(filtered_rsd))
    cat(paste0(number_removed_rsd, " out of ",
               nrow(SummarizedExperiment::assay(se)), " proteins were removed from the dataset due to too high RSD.\n\n"))
    filtered_rsd@metadata$n.filtered.rsd <- number_removed_rsd
  }
  filtered_rsd@metadata$rsd.thresh <- rsd_thresh
  ## Drop proteins with missing values based on defined type filter "filt_thr"
  prot_se <- prot.filter_missing(filtered_rsd, type = filt_type, thr = filt_thr, min = filt_min)

  cat(paste0("Identified conditions:\n ", paste(str_c(unique(prot_se$condition), collapse = ", ")), "\n"))

  return(prot_se)
}

#' Run a complete proteomics analysis workflow.
#'
#' \code{prot.workflow} performs variance stabilization normalization (\code{\link{prot.normalize_vsn}}), missing value imputation (\code{\link{prot.impute}}), principal component analysis (\code{\link{prot.pca}}), differential enrichment test (\code{\link{prot.test_diff}}), and pathway enrichment analysis. If desired, standardized plots and a report are generated and exported as separate files.
#'
#' @param se \code{SummarizedExperiment} object, proteomics data parsed with \code{\link{prot.read_data}}.
#' @param normalize (Logical) Should the data be normalized via variance stabilization normalization?
#' @param imp_fun (Character string)  Function used for data imputation. "SampMin", "man", "bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb", "min", "zero", "mixed", or "nbavg". See (\code{\link{prot.impute}}) for details.
#' @param q (Numeric) q value for imputing missing values with method \code{imp_fun = 'MinProb'}.
#' @param knn.rowmax (Numeric) The maximum percent missing data allowed in any row for \code{imp_fun = 'knn'}. Default: 0.5.
#' @param type (Character string) Type of differential analysis to perform. "all" (contrast each condition with every other condition), "control" (contrast each condition to a defined control condition), "manual" (manually define selected conditions).
#' @param control (Character string) The name of the control condition if \code{type = control}.
#' @param contrast (Character string or vector of strings) Define the contrasts to be tested if \code{type = manual} in the form: \code{"ConditionA_vs_ConditionB"}, or \code{c("ConditionA_vs_ConditionC", "ConditionB_vs_ConditionC")}.
#' @param alpha (Numeric) Significance threshold for adjusted p values.
#' @param alpha_pathways (Numeric) Significance threshold for adjusted p values in pathway enrichment analysis.
#' @param lfc (Numeric) Relevance threshold for log2(fold change) values. Only proteins with a |log2(fold change)| value above \code{lfc} for a given contrast are considered "significant" (if they additionally fullfil the \code{alpha} criterion).
#' @param heatmap.show_all (Logical) Shall all samples be displayed in the heat map (\code{TRUE}) or only the samples contained in the defined \code{contrast} (\code{FALSE})?
#' @param heatmap.kmeans (Logical) Shall the proteins be clustered in the heat map with the k-nearest neighbour method (\code{TRUE}) or not \code{FALSE})?
#' @param k (Integer) Number of protein clusters in heat map if \code{heatmap.kmeans = TRUE}.
#' @param heatmap.col_limit (Integer) Define the outer breaks in the heat map legend. Example: if \code{heatmap.col_limit = 3}, the color scale will span from -3 to 3. Alls values below -3 will have the same color as -3, and all values above 3 will have the same color as 3.
#' @param heatmap.show_row_names (Logical) Show protein names in heat map (\code{TRUE}) or not \code{FALSE}).
#' @param heatmap.row_font_size (Numeric) Font size of protein names in heat maps if \code{heatmap.show_row_names = TRUE}.
#' @param volcano.add_names (Logical) Show protein names in volcano plots (\code{TRUE}) or not \code{FALSE}).
#' @param volcano.label_size (Numeric) Font size of protein names in volcano plots if \code{volcano.add_names = TRUE}.
#' @param volcano.adjusted (Logical) Shall adjusted p values be shown on the y axis of volcano plots (\code{TRUE}) or raw p values (\code{FALSE})?.
#' @param plot (Logical) Show the generated plots in the \code{Plots} pane of RStudio (\code{TRUE}) or not \code{FALSE}).
#' @param export (Logical) Exported the generated plots as PNG and PDF files (\code{TRUE}) or not \code{FALSE}).
#' @param report (Logical) Render and export a report in PDF and HTML format that summarizes the results (\code{TRUE}) or not \code{FALSE}).
#' @param report.dir (Character string) Provide the name of or path to a folder into which the report will be saved.
#' @param pathway_enrichment (Logical) Perform pathway over-representation analysis for each tested contrast (\code{TRUE}) or not \code{FALSE}).
#' @param pathway_kegg (Logical) Perform pathway over-representation analysis with gene sets in the KEGG database (\code{TRUE}) or not \code{FALSE}).
#' @param kegg_organism (Character string) Identifier of the organism in the KEGG database (if \code{pathway_kegg = TRUE})
#' @param custom_pathways (a R dataframe object) Data frame providing custom pathway annotations. The table must contain a "Pathway" column listing identified pathway in the studies organism, and an "Accession" column listing the proteins (or genes) each pathway is composed of. The **Accession** entries must match with protein **IDs**.
#' @param out.dir (Character string) absolute path to the location where result TXT files should be exported to.
#'
#' @return A list containing `SummarizedExperiment` object for every computation step of the workflow, a \code{pca} object, and lists of up- or down regulated pathways for each tested contrast and method (KEGG and/or custom).
#' @export
#'
prot.workflow <- function(se, # SummarizedExperiment, generated with read_prot().
                          normalize = TRUE,
                          imp_fun = c("SampMin", "man", "bpca", "knn", "QRILC", "MLE", "MinDet", # Method for imputing of missing values
                                      "MinProb", "min", "zero", "mixed", "nbavg"),
                          q = 0.01, # q value for imputing missing values with method "fun = 'MinProb'".
                          knn.rowmax = 0.5, # The maximum percent missing data allowed in any row (default 50%).
                          # For any rows with more than rowmax% missing are imputed using the overall mean per sample.
                          type = c("all", "control", "manual"), # Type of differential analysis to perform.
                          control = NULL, # Control condition; required if type = "control".
                          contrast = NULL, # Defined test for differential analysis "A_vs_B"; required if type = "manual".
                          alpha = 0.05, # Significance threshold for adj. p values.
                          alpha_pathways = 0.1,
                          lfc = 1, # Relevance threshold for log2(fold change) values.
                          heatmap.show_all = TRUE, # Shall all samples be displayed in the heatmap or only the samples contained in the defined "contrast"?
                          # (only applicable for type = "manual")
                          heatmap.kmeans = F, # Shall the proteins be clustered in the heat map?
                          k = 6, # Number of protein clusters in heat map if kmeans = TRUE.
                          heatmap.col_limit = NA, # Define the breaks in the heat map legends.
                          heatmap.show_row_names = TRUE, # Show protein names in heat map?
                          heatmap.row_font_size = 6, # Font size of protein names if show_row_names = TRUE.
                          volcano.add_names = FALSE, # Display names next to symbols in volcano plot.
                          volcano.label_size = 2.5, # Size of labels in volcano plot if
                          volcano.adjusted = TRUE, # Display adjusted p-values on y axis of volcano plot?
                          plot = FALSE, # Shall plots be returned in the Plots pane?
                          export = FALSE, # Shall plots be exported as PDF and PNG files?
                          report = TRUE, # Shall a report (HTML and PDF) be created?
                          report.dir = NULL, # Folder name for created report (if report = TRUE)
                          pathway_enrichment = FALSE, # Perform pathway over-representation analysis for each tested contrast
                          pathway_kegg = FALSE, # Perform pathway over-representation analysis with gene sets in the KEGG database
                          kegg_organism = NULL, # Name of the organism in the KEGG database (if 'pathway_kegg = TRUE')
                          custom_pathways = NULL, # Dataframe providing custom pathway annotations
                          out.dir = NULL)
{
  # Show error if inputs are not the required classes
  if(is.integer(alpha)) alpha <- as.numeric(alpha)
  if(is.integer(lfc)) lfc <- as.numeric(lfc)
  assertthat::assert_that(is.character(imp_fun),
                          is.character(type),
                          is.numeric(alpha),
                          length(alpha) == 1,
                          is.numeric(lfc),
                          length(lfc) == 1)

  imp_fun <- match.arg(imp_fun)
  # Show error if inputs are not valid
  if (!type %in% c("all", "control", "manual")) {
    stop("run workflow_proteomics() with a valid type",
         "\nValid types are: 'all', 'control' and 'manual'.",
         call. = FALSE)
  }
  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
  }
  # Variance stabilization
  if(normalize == TRUE){
    prot_norm <- suppressMessages(prot.normalize_vsn(se, plot = plot, export = export))
  } else {
    prot_norm <- se
  }
  # Impute missing values
  if (imp_fun == "MinProb"){
    prot_imp <- prot.impute(prot_norm, fun = imp_fun, q = q)
  } else if (imp_fun == "knn"){
    prot_imp <- prot.impute(prot_norm, fun = imp_fun, rowmax = knn.rowmax)
  } else {
    prot_imp <- prot.impute(prot_norm, fun = imp_fun)
  }
  # Perform PCA Analysis
  prot_pca <- prot.pca(SummarizedExperiment::assay(prot_imp), in_workflow = TRUE)
  # Test for differential expression by empirical Bayes moderation
  # of a linear model and defined contrasts
  prot_diff <- prot.test_diff(prot_imp, type = type, control = control, test = contrast)
  # Denote significantly differentially expressed proteins
  prot_dep <- prot.add_rejections(prot_diff, alpha = alpha, lfc = lfc)
  contrasts <- SummarizedExperiment::rowData(prot_dep) %>%
    data.frame(check.names = FALSE) %>%
    select(ends_with("_diff")) %>%
    colnames() %>% str_replace_all("_diff", "")
  # Generate a results table
  results <- prot.get_results(prot_dep)
  n_significant <- results %>% filter(significant) %>% nrow()

  message(paste0(n_significant,
                 " proteins were found to be differentially expressed with ",
                 expression(alpha),
                 " = ", alpha,
                 " and |log2(fold change)| > ",
                 lfc, "."))
  # Perform pathway enrichment analysis
  if (pathway_enrichment) {
    res.pathway <- enrich_pathways(prot_dep, contrasts, alpha_pathways = alpha_pathways,
                                   pathway_kegg = pathway_kegg, kegg_organism = kegg_organism,
                                   custom_pathways = custom_pathways)
  }
  suppress_ggrepel <- function(w) {
    if (any(grepl("ggrepel", w)))
      invokeRestart("muffleWarning")
  }

  if (export == TRUE |
      plot == TRUE) {
    if (export == TRUE){
      message(paste0("Rendering and exporting figures to:\n",
                     getwd(), "/Plots."))
    }
    suppressMessages(
      prot.boxplot_intensity(se, prot_norm, prot_imp, plot = plot, export = export)
    )
    try(suppressMessages(
      prot.plot_missval(prot_norm, plot = plot, export = export)
    ))
    try(suppressMessages(
      prot.plot_detect(prot_norm, basesize = 10, plot = plot, export = export)
    ))
    suppressMessages(
      prot.plot_imputation(se, prot_norm, prot_imp, plot = plot, export = export, basesize = 12)
    )
    suppressMessages(
      suppressWarnings(
        prot.plot_screeplot(prot_pca, axisLabSize = 18, titleLabSize = 22, plot = plot, export = export)
      ) )
    withCallingHandlers(suppressMessages(
      prot.plot_loadings(prot_pca,labSize = 3, plot = plot, export = export)
    ) , warning = suppress_ggrepel)
    suppressMessages(
      prot.plot_pca(prot_imp,x = 1, y = 2, point_size = 4, basesize = 14, title = "PC Scores - PC2 vs. PC1",
                    plot = plot,export = export)
    )
    suppressMessages(
      prot.plot_pca(prot_imp,x = 1, y = 3, point_size = 4, basesize = 14, title = "PC Scores - PC3 vs. PC1",
                    plot = plot,export = export)
    )
    suppressMessages(
      prot.plot_heatmap(prot_dep, type = "centered", kmeans = heatmap.kmeans, show_all = heatmap.show_all, contrast = contrast,
                        k = k, col_limit = heatmap.col_limit,show_row_names = heatmap.show_row_names,
                        row_font_size = heatmap.row_font_size, indicate = c("condition"),
                        plot = plot, export = export)
    )
    suppressMessages(
      prot.plot_heatmap(prot_dep, type = "contrast", contrast = contrast,  kmeans = heatmap.kmeans, k = k, col_limit = heatmap.col_limit,
                        show_row_names = heatmap.show_row_names, row_font_size = heatmap.row_font_size,
                        plot = plot, export = export)
    )

    for (i in 1:length(contrasts)){
      suppressMessages(
        suppressWarnings(
          prot.plot_volcano(prot_dep, contrast = contrasts[i],
                            add_names = volcano.add_names, label_size = volcano.label_size, adjusted =  volcano.adjusted,
                            plot = plot, export = export, lfc = lfc, alpha = alpha)
        ) )
    }
    if (pathway_kegg) {
      for (i in 1:length(contrasts)) {
        if(!(nrow(as.data.frame(res.pathway$ls.pora_kegg_up[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(res.pathway$ls.pora_kegg_up[[i]], title = paste0("Upregulated pathways", " - KEGG"),
                                   subtitle =  str_replace(contrasts[i], "_vs_", " vs. "), plot = plot, export = export, kegg = TRUE)
            ) )
        }
        if(!(nrow(as.data.frame(res.pathway$ls.pora_kegg_dn[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(res.pathway$ls.pora_kegg_dn[[i]], title = paste0("Downregulated pathways", " - KEGG"),
                                   subtitle = str_replace(contrasts[i], "_vs_", " vs. "), plot = plot, export = export, kegg = TRUE)
            ) )
        }
      }
    }
    if(!is.null(custom_pathways)){
      for (i in 1:length(contrasts)) {
        if(!(nrow(as.data.frame(res.pathway$ls.pora_custom_up[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(res.pathway$ls.pora_custom_up[[i]], title = "Upregulated pathways",
                                   subtitle =  str_replace(contrasts[i], "_vs_", " vs. "), plot = plot, export = export, kegg = FALSE)
            ) )
        }
        if(!(nrow(as.data.frame(res.pathway$ls.pora_custom_dn[[i]])) == 0)){
          suppressMessages(
            suppressWarnings(
              prot.plot_enrichment(res.pathway$ls.pora_custom_dn[[i]], title = "Downregulated pathways",
                                   subtitle = str_replace(contrasts[i], "_vs_", " vs. "), plot = plot, export = export, kegg = FALSE)
            ) )
        }
      }
    }
  }

if(!is.null(out.dir)){
  out_dir <- out.dir
  dir.create("out_dir", showWarnings = F)
} else {
  out_dir <- getwd()
}



  param <- data.frame(type, alpha, lfc, check.names = FALSE)
  results <- list(data = SummarizedExperiment::rowData(se), se = se, norm = prot_norm,
                  imputed = prot_imp, pca = prot_pca, diff = prot_diff, dep = prot_dep,
                  results = results, param = param)

  message(paste0("Save results as tab-delimited table to: ", getwd(), "/results.txt"))
  utils::write.table(results$results, paste(getwd(), "results.txt",
                                  sep = "/"), row.names = FALSE, sep = "\t")

  message("Save RData object")
  save(results, file = paste(out_dir, "results.RData", sep = "/"))

  if (pathway_enrichment == T && pathway_kegg) {
    results <- c(results, pora_kegg_up = list(res.pathway$ls.pora_kegg_up), pora_kegg_dn = list(res.pathway$ls.pora_kegg_dn))
    message(paste0("Writing results of KEGG pathway enrichment analysis to: ", out_dir, "'pora_kegg_contrast...txt'"))
    for(i in 1:length(res.pathway$ls.pora_kegg_up)){
      if(!is.null(res.pathway$ls.pora_kegg_up[[i]])){
        utils::write.table(res.pathway$ls.pora_kegg_up[[i]]@result, paste(out_dir, paste0("pora_kegg_", names(res.pathway$ls.pora_kegg_up)[i], "_up.txt"),
                                                sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
    for(i in 1:length(res.pathway$ls.pora_kegg_dn)){
      if(!is.null(res.pathway$ls.pora_kegg_dn[[i]])){
        utils::write.table(res.pathway$ls.pora_kegg_dn[[i]]@result, paste(out_dir, paste0("pora_kegg_", names(res.pathway$ls.pora_kegg_dn)[i], "_down.txt"),
                                                sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
  }
  if (pathway_enrichment == T && !is.null(custom_pathways)) {
    results <- c(results, pora_custom_up = list(res.pathway$ls.pora_custom_up), pora_custom_dn = list(res.pathway$ls.pora_custom_dn))
    message(paste0("Writing results of custom pathway enrichment analysis to: ", out_dir, "'pora_custom_contrast...txt'"))
    for(i in 1:length(res.pathway$ls.pora_custom_up)){
      if(!is.null(res.pathway$ls.pora_custom_up[[i]])){
        utils::write.table(res.pathway$ls.pora_custom_up[[i]]@result, paste(out_dir, paste0("pora_custom_", names(res.pathway$ls.pora_custom_up)[i], "_up.txt"),
                                                sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
    for(i in 1:length(res.pathway$ls.pora_custom_dn)){
      if(!is.null(res.pathway$ls.pora_custom_dn[[i]])){
        utils::write.table(res.pathway$ls.pora_custom_dn[[i]]@result, paste(out_dir, paste0("pora_custom_", names(res.pathway$ls.pora_custom_dn)[i], "_down.txt"),
                                                sep = "/"), row.names = FALSE, sep = "\t")
      }
    }
  }

  if(report == TRUE){
    prot.report(results,
                volcano.adjusted = volcano.adjusted,
                pathway_enrichment = pathway_enrichment,
                heatmap.show_all = heatmap.show_all,
                heatmap.kmeans = heatmap.kmeans,
                volcano.add_names = volcano.add_names,
                k = k,
                report.dir = report.dir)
  }
  return(results)
}

#' Impute missing values in a SummarizedExperiment object
#'
#' Imputes missing values in a SummarizedExperiment object.
#'
#' @param se A SummarizedExperiment object
#' @param fun A character string specifying the imputation method to use.
#'   The available methods are "bpca", "knn", "QRILC", "MLE", "MinDet",
#'   "MinProb", "man", "min", "zero", "mixed", "nbavg", and "SampMin".
#' @param ... Additional arguments passed to the imputation function.
#'
#' @return A SummarizedExperiment object with imputed values
#'
#' @export
#'
#' @details "SampMin" replaces missing values with the minimum value found in each sample. For information about the remaining imputation methods, see \code{help("imputeMethods", "MsCoreUtils")}
#'
#' @examples
#' data(mset)
#' se <- make_se(mset)
#' se <- prot.impute(se, fun = "knn")
#'
#' @importFrom assertthat assert_that
#' @importFrom SummarizedExperiment rowData assay
#' @importFrom MSnbase exprs impute
#'
#' @seealso \code{\link{prot.make_unique}}, \code{\link{prot.make_se}}
#'
#' @keywords imputation, SummarizedExperiment
prot.impute <- function (se, fun = c("bpca", "knn", "QRILC",
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
  if (fun == "man") {
    se <- manual_impute(se, ...)
  } else if (fun == "SampMin"){
    SummarizedExperiment::assay(se) <- SummarizedExperiment::assay(se) %>% data.frame() %>%
      mutate_if(is.numeric, function(x) ifelse(is.na(x), min(x, na.rm = T), x)) %>%
      as.matrix()
  } else {
    MSnSet_data <- methods::as(se, "MSnSet")
    MSnSet_imputed <- MSnbase::impute(MSnSet_data, method = fun,
                                      ...)
    SummarizedExperiment::assay(se) <- MSnbase::exprs(MSnSet_imputed)
  }
  se@metadata["imp_fun"] <- fun
  return(se)
}

#' Normalize the data via variance stabilization normalization
#'
#' @param se A SummarizedExperiment object
#' @param plot Logical. If TRUE, plots the meanSdPlot
#' @param export Logical. If TRUE, exports the meanSdPlot in pdf and png
#' @return A SummarizedExperiment object
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom vsn vsnMatrix predict
#' @importFrom SummarizedExperiment assay
#' @importFrom grDevices pdf png
#'
prot.normalize_vsn <- function (se, plot = TRUE, export = TRUE)
  {
  # Normalize the data (including log2 transformation)
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  se_vsn <- se
  vsn.fit <- suppressMessages(vsn::vsnMatrix(2^SummarizedExperiment::assay(se_vsn)), classes = "message") # Fit the vsn model
  SummarizedExperiment::assay(se_vsn) <- vsn::predict(vsn.fit, 2^SummarizedExperiment::assay(se_vsn)) #Apply the vsn transformation to data

  # Verify the variance stabilisation.
  # "The aim of these plots is to see whether there is a systematic trend in the standard
  # deviation of the data as a function of overall expression. The assumption that
  # underlies the usefulness of these plots is that most genes are not differentially
  # expressed, so that the running median is a reasonable estimator of the standard
  # deviation of feature level data conditional on the mean. After variance stabilisation,
  # this should be approximately a horizontal line. It may have some random fluctuations,
  # but should not show an overall trend. If this is not the case, that usually indicates a
  # data quality problem, or is a consequence of inadequate prior data preprocessing."
  # (source: https://bioconductor.org/packages/release/bioc/vignettes/vsn/inst/doc/A-vsn.html )
  if (export == TRUE){
    dir.create(paste0(getwd(), "/Plots"), showWarnings = F)
    grDevices::pdf("Plots/meanSDPlot.pdf")
    plot_meanSdPlot <- suppressWarnings(meanSdPlot(se_vsn, plot = T, xlab = "Rank(mean)", ylab = "SD"))
    grDevices::dev.off()

    grDevices::png("Plots/meanSDPlot.png",
        width = 6, height = 6, units = 'in', res = 300)
    plot_meanSdPlot <- suppressWarnings(meanSdPlot(se_vsn, plot = T, xlab = "Rank(mean)", ylab = "SD"))
    grDevices::dev.off()
    message(paste0("Exporting meanSdPlot to:\n\"", getwd(), "\"/Plots/meanSdPlot.pdf\" and \".../meanSdPlot.png\""))
  }

  if (plot == TRUE){
    meanSdPlot(se_vsn, plot = T, xlab = "Rank(mean)", ylab = "SD")
  }


  return(se_vsn)
}

#' Constructor for ExactParam objects
#'
#' @param deferred Logical indicating whether to use deferred evaluations
#' @param fold Numeric indicating the fold value
#' @return An ExactParam object
#' @export
#'
#' @importFrom utils head
#'
#' @noRd
#'
#' @keywords internal
ExactParam <- function (deferred = FALSE, fold = Inf)
{
  new("ExactParam", deferred = as.logical(deferred),
      fold = as.numeric(fold))
}

#' @title PCA Analysis
#'
#' @description
#' This function performs principal component analysis (PCA) on a given matrix.
#'
#' @param mat A numeric matrix.
#' @param metadata An optional data frame with rownames matching 'colnames(mat)'.
#' @param center A logical indicating whether to center the data before PCA.
#' @param scale A logical indicating whether to scale the data before PCA.
#' @param rank An integer indicating the number of principal components to calculate.
#' @param removeVar A numeric value between 0 and 1 indicating the fraction of variables to remove based on variance.
#' @param transposed A logical indicating whether the matrix is transposed.
#' @param BSPARAM An object of class \code{ExactParam} or \code{ApproxParam}.
#' @param in_workflow Indicated whether this function is called from within `prot.workflow`.
#'
#' @return A list with components \code{rotated}, \code{loadings}, \code{variance}, \code{sdev}, \code{metadata}, \code{xvars}, \code{yvars}, and \code{components}.
#'
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' prot.pca(mat)
#'
#' @references
#' BiocSingular (https://bioconductor.org/packages/release/bioc/html/BiocSingular.html)
#'
prot.pca <- function (mat, metadata = NULL, center = TRUE, scale = FALSE,
                      rank = NULL, removeVar = NULL, transposed = FALSE, BSPARAM = ExactParam(), in_workflow = FALSE)
{
  if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  }
  if (!transposed) {
    mat <- t(mat)
  }
  if (!is.null(metadata)) {
    if (!identical(rownames(mat), rownames(metadata))) {
      stop("'colnames(mat)' is not identical to 'rownames(metadata)'")
    }
  }
  .center <- if (center)
    NULL
  else 0
  vars <- matrixStats::colVars(as.matrix(DelayedArray::DelayedArray(mat)), center = .center)
  if (!is.null(removeVar)) {
    message("-- removing the lower ", removeVar * 100,
            "% of variables based on variance")
    varorder <- order(vars, decreasing = TRUE)
    keep <- utils::head(varorder, max(1, ncol(mat) * (1 - removeVar)))
    mat <- mat[, keep, drop = FALSE]
    vars <- vars[keep]
  }
  if (is.null(rank)) {
    if (in_workflow == TRUE || is(BSPARAM, "ExactParam")) {
      rank <- min(dim(mat))
    }
    else {
      stop("'rank' must be specified for approximate PCA methods")
    }
  }
  pcaobj <- BiocSingular::runPCA(mat, center = center, scale = scale, rank = rank,
                                 BSPARAM = BSPARAM)
  if (scale) {
    total.var <- length(vars)
  }
  else {
    total.var <- sum(vars)
  }
  proportionvar <- (pcaobj$sdev^2)/total.var * 100
  pcaobj <- list(rotated = data.frame(pcaobj$x), loadings = data.frame(pcaobj$rotation),
                 variance = proportionvar, sdev = pcaobj$sdev, metadata = metadata,
                 xvars = colnames(mat), yvars = rownames(mat), components = colnames(pcaobj$x))
  rownames(pcaobj$rotated) <- pcaobj$yvars
  rownames(pcaobj$loadings) <- pcaobj$xvars
  names(pcaobj$variance) <- pcaobj$components
  class(pcaobj) <- "pca"
  return(pcaobj)
}


#' Test Differential Expression
#'
#' This function tests for differential expression between conditions
#'
#' @param se A SummarizedExperiment object
#' @param type Character vector indicating which type of contrast to use. Options are: "control", "all" or "manual".
#' @param control Character vector indicating the control condition. Required when type is "control".
#' @param test Character vector indicating the contrasts to be tested. Required when type is "manual". Needs to be in the form "Subject_vs_Control".
#' @param design_formula Formula object indicating the design of the experiment.
#' @return A SummarizedExperiment object with the results of the differential expression test that can be accessed with \code{ SummarizedExperiment::rowData(se)}
#' @export
#'
#' @importFrom stats terms.formula
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom fdrtool fdrtool
#' @importFrom purrr map_df
#' @importFrom tidyr spread
#' @importFrom dplyr mutate select
#' @importFrom utils combn
#'
#' @importFrom SummarizedExperiment assay
prot.test_diff <- function (se, type = c("control", "all", "manual"),
                            control = NULL, test = NULL, design_formula = stats::formula(~0 +
                                                                                    condition))
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(type), class(design_formula) == "formula")
  type <- match.arg(type)
  col_data <- SummarizedExperiment::colData(se)
  raw <- SummarizedExperiment::assay(se)
  if (any(!c("name", "ID") %in% colnames(SummarizedExperiment::rowData(se,
                                                 use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)), "'\nRun prot.make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in%
          colnames(col_data))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)), "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  if (any(is.na(raw))) {
    warning("Missing values in '", deparse(substitute(se)),
            "'")
  }
  if (!is.null(control)) {
    assertthat::assert_that(is.character(control), length(control) ==
                              1)
    if (!control %in% unique(col_data$condition)) {
      stop("run test_diff() with a valid control.\nValid controls are: '",
           paste0(unique(col_data$condition), collapse = "', '"),
           "'", call. = FALSE)
    }
  }
  variables <- stats::terms.formula(design_formula) %>% attr(., "variables") %>%
    as.character() %>% .[-1]
  if (any(!variables %in% colnames(col_data))) {
    stop("run make_diff() with an appropriate 'design_formula'")
  }
  if (variables[1] != "condition") {
    stop("first factor of 'design_formula' should be 'condition'")
  }
  for (var in variables) {
    temp <- factor(col_data[[var]])
    assign(var, temp)
  }
  design <- model.matrix(design_formula, data = environment())
  colnames(design) <- gsub("condition", "", colnames(design))
  conditions <- as.character(unique(condition))
  if (type == "all") {
    cntrst <- apply(utils::combn(conditions, 2), 2, paste,
                    collapse = " - ")
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
  if (type == "control") {
    if (is.null(control))
      stop("run test_diff(type = 'control') with a 'control' argument")
    cntrst <- paste(conditions[!conditions %in% control],
                    control, sep = " - ")
  }
  if (type == "manual") {
    if (is.null(test)) {
      stop("run test_diff(type = 'manual') with a 'test' argument")
    }
    assertthat::assert_that(is.character(test))
    if (any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
      stop("run test_diff() with valid contrasts in 'test'",
           ".\nValid contrasts should contain combinations of: '",
           paste0(conditions, collapse = "', '"),
           "', for example '", paste0(conditions[1],
                                      "_vs_", conditions[2]), "'.", call. = FALSE)
    }
    cntrst <- gsub("_vs_", " - ", test)
  }
  cat("Tested contrasts: ", paste(gsub(" - ",
                                       "_vs_", cntrst), collapse = ", "), "\n")
  fit <- limma::lmFit(raw, design = design)
  made_contrasts <- limma::makeContrasts(contrasts = cntrst, levels = design)
  contrast_fit <- limma::contrasts.fit(fit, made_contrasts)
  if (any(is.na(raw))) {
    for (i in cntrst) {
      covariates <- strsplit(i, " - ") %>% unlist
      single_contrast <- limma::makeContrasts(contrasts = i, levels = design[,
                                                                      covariates])
      single_contrast_fit <- limma::contrasts.fit(fit[, covariates],
                                           single_contrast)
      contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[,
                                                                         1]
      contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[,
                                                                             1]
    }
  }
  eB_fit <- limma::eBayes(contrast_fit)
  retrieve_fun <- function(comp, fit = eB_fit) {
    res <- limma::topTable(fit, sort.by = "t", coef = comp,
                           number = Inf, confint = TRUE)
    res <- res[!is.na(res$t), ]
    fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
    res$qval <- fdr_res$qval
    res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- tibble::rownames_to_column(res)
    return(res)
  }
  limma_res <- purrr::map_df(cntrst, retrieve_fun)
  table <- limma_res %>% select(rowname, logFC, CI.L, CI.R,
                                P.Value, qval, comparison) %>% dplyr::mutate(comparison = gsub(" - ",
                                                                                        "_vs_", comparison)) %>% gather(variable, value,
                                                                                                                        -c(rowname, comparison)) %>% dplyr::mutate(variable = recode(variable,
                                                                                                                                                                              logFC = "diff", P.Value = "p.val", qval = "p.adj")) %>%
    unite(temp, comparison, variable) %>% tidyr::spread(temp, value)
  SummarizedExperiment::rowData(se) <- merge(SummarizedExperiment::rowData(se, use.names = FALSE), table,
                       by.x = "name", by.y = "rowname", all.x = TRUE,
                       sort = FALSE)
  return(se)
}



#' Add rejections to a SummarizedExperiment
#'
#' This function adds a logical vector of rejections (TRUE/FALSE) to a SummarizedExperiment object.
#'
#' @param diff A SummarizedExperiment object.
#' @param alpha A numeric value for the alpha level (adj. p value) threshold.
#' @param lfc A numeric value for the log2 fold change threshold.
#'
#' @return A SummarizedExperiment object. The results of this function can be accessed with \code{ SummarizedExperiment::rowData(se)}
#' @export
#'
#'
#' @importFrom assertthat assert_that
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData
#'
prot.add_rejections <- function (diff, alpha = 0.05, lfc = 1)
{
  if (is.integer(alpha))
    alpha <- as.numeric(alpha)
  if (is.integer(lfc))
    lfc <- as.numeric(lfc)
  assertthat::assert_that(inherits(diff, "SummarizedExperiment"),
                          is.numeric(alpha), length(alpha) == 1, is.numeric(lfc),
                          length(lfc) == 1)
  row_data <- SummarizedExperiment::rowData(diff, use.names = FALSE) %>% as.data.frame()
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(diff)), "'\nRun prot.make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) <
      1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(diff)), "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }
  cols_p <- grep("_p.adj", colnames(row_data))
  cols_diff <- grep("_diff", colnames(row_data))
  if (length(cols_p) == 1) {
    SummarizedExperiment::rowData(diff)$significant <- row_data[, cols_p] <= alpha &
      abs(row_data[, cols_diff]) >= lfc
    SummarizedExperiment::rowData(diff)$contrast_significant <- SummarizedExperiment::rowData(diff, use.names = FALSE)$significant
    colnames(SummarizedExperiment::rowData(diff))[ncol(SummarizedExperiment::rowData(diff, use.names = FALSE))] <- gsub("p.adj",
                                                                            "significant", colnames(row_data)[cols_p])
  }
  if (length(cols_p) > 1) {
    p_reject <- row_data[, cols_p] <= alpha
    p_reject[is.na(p_reject)] <- FALSE
    diff_reject <- abs(row_data[, cols_diff]) >= lfc
    diff_reject[is.na(diff_reject)] <- FALSE
    sign_df <- p_reject & diff_reject
    sign_df <- cbind(sign_df, significant = apply(sign_df,
                                                  1, function(x) any(x)))
    colnames(sign_df) <- gsub("_p.adj", "_significant",
                              colnames(sign_df))
    sign_df <- cbind(name = row_data$name, as.data.frame(sign_df))
    SummarizedExperiment::rowData(diff) <- merge(SummarizedExperiment::rowData(diff, use.names = FALSE),
                           sign_df, by = "name")
  }
  diff@metadata["alpha"] <- as.numeric(alpha)
  diff@metadata["lfc"] <- as.numeric(lfc)
  return(diff)
}




#' @title Get Results
#' @description Function to generate a table with results of a differential expression analysis
#' @param dep A SummarizedExperiment object
#' @return A data frame containing the results
#'
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom SummarizedExperiment rowData assay colData
#' @importFrom dplyr group_by summarize mutate
#' @importFrom tidyr spread
#' @importFrom tibble column_to_rownames
prot.get_results <- function (dep)
{
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"))
  row_data <- SummarizedExperiment::rowData(dep, use.names = FALSE)
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(dep)), "'\nRun prot.make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) <
      1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(dep)), "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }
  row_data$mean <- rowMeans(SummarizedExperiment::assay(dep), na.rm = TRUE)
  centered <- SummarizedExperiment::assay(dep) - row_data$mean
  centered <- data.frame(centered) %>% tibble::rownames_to_column() %>%
    gather(ID, val, -rowname) %>% left_join(., data.frame(SummarizedExperiment::colData(dep)),
                                            by = "ID")
  centered <- group_by(centered, rowname, condition) %>% dplyr::summarize(val = mean(val,
                                                                              na.rm = TRUE)) %>% dplyr::mutate(val = signif(val, digits = 3)) %>%
    tidyr::spread(condition, val)
  colnames(centered)[2:ncol(centered)] <- paste(colnames(centered)[2:ncol(centered)],
                                                "_centered", sep = "")
  ratio <- as.data.frame(row_data) %>% tibble::column_to_rownames("name") %>%
    select(ends_with("diff")) %>% signif(., digits = 3) %>%
    tibble::rownames_to_column()
  colnames(ratio)[2:ncol(ratio)] <- gsub("_diff", "_log2fc",
                                         colnames(ratio)[2:ncol(ratio)])
  df <- left_join(ratio, centered, by = "rowname")
  pval <- as.data.frame(row_data) %>% tibble::column_to_rownames("name") %>%
    select(ends_with("p.val"), ends_with("p.adj"), ends_with("imputed"),
           ends_with("significant")) %>% tibble::rownames_to_column()
  pval[, grep("p.adj", colnames(pval))] <- pval[, grep("p.adj",
                                                       colnames(pval))] %>% signif(digits = 3)
  ids <- as.data.frame(row_data) %>% select(name, ID)
  table <- left_join(ids, pval, by = c(name = "rowname"))
  table <- left_join(table, df, by = c(name = "rowname")) %>%
    arrange(desc(significant))
  return(table)
}




#' Constructor Method for HeatmapAnnotation class using ComplexHeatmap::HeatmapAnnotation
#'
#' @param dep SummarizedExperiment object
#' @param indicate Character vector indicating the column names of the annotation
#' @return HeatmapAnnotation object
#' @keywords internal
#' @noRd
get_annotation <- function (dep, indicate)
{
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(indicate))
  col_data <- SummarizedExperiment::colData(dep) %>% as.data.frame()
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




#' @title Report Proteomics Data
#' @description Generates a report with proteomics data in PDF and HTML format
#' @param results A list containing the results of a differential protein expression analysis, generated with \code{prot.workflow}.
#' @param report.dir A character string specifying the directory where the report will be saved
#' @param ... Additional parameters to be passed to markdown to generate the report.
#' currently supported: volcano.adjusted, pathway_enrichment, heatmap.show_all, heatmap.kmeans,volcano.add_names, k = k
#' @export
prot.report <- function(results, report.dir = NULL, ...)
{
  assertthat::assert_that(is.list(results))
  if (any(!c("data", "se", "norm",
             "imputed", "diff", "dep", "results",
             "param") %in% names(results))) {
    stop("run prot.report() with appropriate input generated by prot.workflow",
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
  dep <- results$dep
  param <- results$param
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
    wd <- report.dir
  } else {
    wd <- paste(getwd(), "/Report.prot_", format(Sys.time(),
                                                 "%Y%m%d_%H%M%S"), sep = "")
  }
  dir.create(wd, showWarnings = F)
  for(i in 1:length(.libPaths())){
    VisomX.ndx <- grep("VisomX", list.files(.libPaths()[i]))
    if(length(VisomX.ndx)>0){
      Report.wd <- paste0(.libPaths()[i], "/VisomX")
    }
  }
  file <- paste0(Report.wd, "/Report_Prot.Rmd")
  message("Render reports...")
  rmarkdown::render(file, output_format = "all", output_dir = wd,
                    quiet = TRUE)

  message(paste0("Files saved in: '", wd, "'"))
}





#' @title Pathway Enrichment Analysis
#'
#' @description
#' This function performs pathway enrichment analysis. Uses the DOSE package with the internal function \code{enricher_internal}.
#'
#' @param gene A numeric vector or a character vector of gene symbols that were found enriched in a differential expression analysis.
#' @param organism A character string. The organism to use. Default is "ppu" for Pseudomonas putida KT2440.
#' @param keyType A character string. The type of gene ID. Default is "kegg".
#' @param pvalueCutoff A numeric value. The cutoff of p-value. Default is 0.05.
#' @param pAdjustMethod A character string. The method of p-value adjustment. Default is "BH".
#' @param universe A numeric vector or a character vector of gene symbols.
#' @param minGSSize A numeric value. The minimum gene set size. Default is 10.
#' @param maxGSSize A numeric value. The maximum gene set size. Default is 500.
#' @param qvalueCutoff A numeric value. The cutoff of q-value. Default is 0.2.
#' @param use_internal_kegg A logical value. Whether to use internal KEGG database. Default is FALSE.
#' @param custom_gene_sets A logical value. Whether to use custom gene sets. Default is FALSE.
#' @param custom_pathways A data frame. The custom pathways to use. Needs to have column "Pathway" and corresponding "Accession" column with comma-separated entries of gene/protein identifiers for that pathway.
#'
#' @return A data frame containing the results of the pathway enrichment analysis.
#'
#' @export
#'
#' @importFrom Biobase reverseSplit
pathway_enrich <- function (gene, organism = "ppu", keyType = "kegg",
                                 pvalueCutoff = 0.05, pAdjustMethod = "BH", universe,
                                 minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_kegg = FALSE, custom_gene_sets = FALSE, custom_pathways = NULL)
{
  if (custom_gene_sets) {
    custom_pathways$Accession <-
      custom_pathways$Accession %>% str_replace_all(., " // ", ", ")
    custom_vec <-
      custom_pathways[, match(c("Pathway", "Accession"), colnames(custom_pathways))] %>% tibble::deframe()
    NAME2EXTID <- strsplit(as.character(custom_vec), ", ")
    names(NAME2EXTID) <- names(custom_vec)
    EXTID2NAME <- Biobase::reverseSplit(NAME2EXTID)
    DATA <- new.env()
    DATA$NAME2EXTID  <- NAME2EXTID
    DATA$EXTID2NAME  <- EXTID2NAME
    DATA$PATHID2EXTID <- NAME2EXTID
    res <- enricher_custom(gene, pvalueCutoff = pvalueCutoff,
                                pAdjustMethod = pAdjustMethod, universe = universe, minGSSize = minGSSize,
                                maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff, USER_DATA = DATA)
    if (is.null(res))
      return(res)
  } else {
    species <- clusterProfiler:::organismMapper(organism)
    if (use_internal_kegg) {
      DATA <- clusterProfiler:::get_data_from_KEGG_db(species)
    } else {
      DATA <- clusterProfiler:::prepare_KEGG(species, "KEGG", keyType)
    }
    res <- DOSE:::enricher_internal(gene, pvalueCutoff = pvalueCutoff,
                                    pAdjustMethod = pAdjustMethod, universe = universe, minGSSize = minGSSize,
                                    maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff, USER_DATA = DATA)
    if (is.null(res))
      return(res)
    res@ontology <- "KEGG"
    res@organism <- species
    res@keytype <- keyType
  }
  res <- dplyr::mutate(res, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
  return(res)
}





#' Enricher Custom
#'
#' This internal function performs pathway enrichment analysis using a custom set of pathways and genes/proteins.
#'
#' @param gene A character vector of gene or protein IDs found to be enriched.
#' @param pvalueCutoff The cutoff p-value for enrichment.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param universe background genes/proteins
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of each geneSet for analyzing.
#' @param qvalueCutoff The cutoff q-value (adjusted p value) for enrichment.
#' @param USER_DATA ontology information.
#'
#' @return An enriched results objects of class \code{enrichResult}
#'
#' @noRd
#' @keywords internal
enricher_custom <- function (gene, pvalueCutoff, pAdjustMethod = "BH", universe = NULL,
                                  minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, USER_DATA)
{
  gene <- as.character(unique(gene))

  EXTID2NAME <- get("EXTID2NAME", envir = USER_DATA)
  qExtID2Name <- EXTID2NAME[gene]
  len <- sapply(qExtID2Name, length)
  notZero.idx <- len != 0
  qExtID2TermID <- qExtID2Name[notZero.idx]
  qTermID <- unlist(qExtID2TermID)
  if (is.null(qTermID)) {
    message("--> No gene can be mapped. Either no gene was found as part of any annotated custom pathway or gene IDs have the wrong format.")
    p2e <- get("NAME2EXTID", envir = USER_DATA)
    sg <- unlist(p2e[1:10])
    sg <- sample(sg, min(length(sg), 6))
    message("--> Expected input gene IDs (examples): ", paste0(sg,
                                                               collapse = ","))
    return(NULL)

  }
  qExtID2TermID.df <- data.frame(extID = rep(names(qExtID2TermID),
                                             times = lapply(qExtID2TermID, length)), termID = qTermID)
  qExtID2TermID.df <- unique(qExtID2TermID.df)
  extID <- get("NAME2EXTID", envir = USER_DATA) %>% unlist() %>% unique()

  qTermID2ExtID <- with(qExtID2TermID.df, split(as.character(extID),
                                                as.character(termID)))
  if (missing(universe)){
    universe <- NULL
  }
  if (!is.null(universe)) {
    if (is.character(universe)) {
      extID <- intersect(extID, universe)
    }
    else {
      message("`universe` is not in character and will be ignored...")
    }
  }
  qTermID2ExtID <- lapply(qTermID2ExtID, intersect, extID)
  qTermID <- unique(names(qTermID2ExtID))

  termID2ExtID <- get("PATHID2EXTID", envir = USER_DATA)
  termID2ExtID <- termID2ExtID[qTermID]
  termID2ExtID <- lapply(termID2ExtID, intersect, extID)
  geneSets <- lapply(termID2ExtID, intersect, extID)

  if (is.na(minGSSize) || is.null(minGSSize))
    minGSSize <- 1
  if (is.na(maxGSSize) || is.null(maxGSSize))
    maxGSSize <- Inf
  geneSet_size <- sapply(geneSets, length)
  idx <- minGSSize <= geneSet_size & geneSet_size <= maxGSSize


  if (sum(idx) == 0) {
    msg <- paste("No gene sets have size between",
                 minGSSize, "and", maxGSSize, "...")
    message(msg)
    message("--> return NULL...")
    return(NULL)
  }
  termID2ExtID <- termID2ExtID[idx]
  qTermID2ExtID <- qTermID2ExtID[idx]
  qTermID <- unique(names(qTermID2ExtID))
  k <- sapply(qTermID2ExtID, length)
  k <- k[qTermID]
  M <- sapply(termID2ExtID, length)
  M <- M[qTermID]
  N <- rep(length(extID), length(M))
  n <- rep(length(qExtID2TermID), length(M))
  args.df <- data.frame(numWdrawn = k - 1, numW = M, numB = N -
                          M, numDrawn = n)
  pvalues <- apply(args.df, 1, function(n) stats::phyper(n[1], n[2],
                                                  n[3], n[4], lower.tail = FALSE))
  GeneRatio <- apply(data.frame(a = k, b = n), 1, function(x) paste(x[1],
                                                                    "/", x[2], sep = "", collapse = ""))
  BgRatio <- apply(data.frame(a = M, b = N), 1, function(x) paste(x[1],
                                                                  "/", x[2], sep = "", collapse = ""))
  Over <- data.frame(ID = as.character(qTermID), GeneRatio = GeneRatio,
                     BgRatio = BgRatio, pvalue = pvalues, stringsAsFactors = FALSE)
  p.adj <- stats::p.adjust(Over$pvalue, method = pAdjustMethod)
  qobj <- tryCatch(qvalue::qvalue(p = Over$pvalue, lambda = 0.05, pi0.method = "bootstrap"),
                   error = function(e) NULL)
  if (class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  } else {
    qvalues <- NA
  }
  geneID <- sapply(qTermID2ExtID, function(i) paste(i, collapse = "/"))
  geneID <- geneID[qTermID]
  Over <- data.frame(Over, p.adjust = p.adj, qvalue = qvalues,
                     geneID = geneID, Count = k, stringsAsFactors = FALSE)
  Description <- qTermID
  if (length(qTermID) != length(Description)) {
    idx <- qTermID %in% names(Description)
    Over <- Over[idx, ]
  }
  Over$Description <- Description
  nc <- ncol(Over)
  Over <- Over[, c(1, nc, 2:(nc - 1))]
  Over <- Over[order(pvalues), ]
  Over$ID <- as.character(Over$ID)
  row.names(Over) <- as.character(Over$ID)
  x <- new("enrichResult", result = Over, pvalueCutoff = pvalueCutoff,
           pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,
           gene = as.character(gene), universe = extID, geneSets = geneSets,
           organism = "UNKNOWN", keytype = "UNKNOWN",
           ontology = "UNKNOWN", readable = FALSE)
  return(x)
}





#' Create a SummarizedExperiment object from a data frame
#'
#' This function creates a SummarizedExperiment object from a data frame containing protein expression data.
#'
#' @param proteins_unique A data frame containing the protein expression data
#' @param columns A vector of column numbers containing the numeric expression data
#' @param expdesign A data frame containing the experimental design
#'
#' @return A SummarizedExperiment object
#'
#' @export
#'
make_se <- function (proteins_unique, columns, expdesign)
{
  assertthat::assert_that(is.data.frame(proteins_unique), is.integer(columns),
                          is.data.frame(expdesign))
  if (any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)), "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns",
         "are not present in the experimental design", call. = FALSE)
  }
  if (any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric", "\nRun make_se_parse() with the appropriate columns as argument",
         call. = FALSE)
  }
  if (tibble::is_tibble(proteins_unique))
    proteins_unique <- as.data.frame(proteins_unique)
  if (tibble::is_tibble(expdesign))
    expdesign <- as.data.frame(expdesign)
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)
  expdesign <- dplyr::mutate(expdesign, condition = make.names(condition)) %>%
    tidyr::unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID
  matched <- match(make.names(delete_prefix(expdesign$label)),
                   make.names(delete_prefix(colnames(raw))))
  if (any(is.na(matched))) {
    stop("None of the labels in the experimental design match ",
         "with column names in 'proteins_unique'", "\nRun make_se() with the correct labels in the experimental design ",
         "and/or correct columns specification")
  }
  colnames(raw)[matched] <- expdesign$ID
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$name
  se <- SummarizedExperiment::SummarizedExperiment(assays = as.matrix(raw), colData = expdesign,
                                                   rowData = row_data)
  return(se)
}

#' Delete Prefix
#'
#' @param words A character vector of words
#'
#' @return A character vector of words with a common prefix
#'
#' @examples
#' delete_prefix(c("unhappy", "unfair"))
#'
#' @export
#'
delete_prefix <- function (words)
{
  prefix <- get_prefix(words)
  gsub(paste0("^", prefix), "", words)
}



#' Get the common prefix of words
#'
#' @param words A character vector of words
#'
#' @return A character string that is the common prefix of all words
#'
#' @export
#'
#' @examples
#' get_prefix(c("animal", "antelope", "anteater"))
#'
get_prefix <- function (words)
{
  assertthat::assert_that(is.character(words))
  if (length(words) <= 1) {
    stop("'words' should contain more than one element")
  }
  if (any(is.na(words))) {
    stop("'words' contains NAs")
  }
  minlen <- min(nchar(words))
  truncated <- substr(words, 1, minlen)
  if (minlen < 1) {
    stop("At least one of the elements is too short")
  }
  mat <- data.frame(strsplit(truncated, ""), stringsAsFactors = FALSE)
  identical <- apply(mat, 1, function(x) length(unique(x)) ==
                       1)
  prefix <- as.logical(cumprod(identical))
  paste(mat[prefix, 1], collapse = "")
}




#' Filter SummarizedExperiment object based on missing values
#'
#' Will keep proteins that have a maximum of \code{thr} missing values in at least one condition.
#'
#' @param se A SummarizedExperiment object
#' @param thr A numeric value representing the missing value threshold
#'
#' @return A filtered SummarizedExperiment object
#'
#' @export
#'
#' @details Used with option \code{type = 'condition'} in \code{\link{prot.filter_missing()}}.
#'
#' @importFrom assertthat assert_that
#' @importFrom SummarizedExperiment rowData colData assay
#' @importFrom dplyr summarize
#' @importFrom tidyr spread
#'
#' @export filter_missval
filter_missval <- function (se, thr = 0)
{
  if (is.integer(thr))
    thr <- as.numeric(thr)
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.numeric(thr), length(thr) == 1)
  if (any(!c("name", "ID") %in% colnames(SummarizedExperiment::rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)), "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(SummarizedExperiment::colData(se)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)), "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  max_repl <- max(SummarizedExperiment::colData(se)$replicate)
  if (thr < 0 | thr > max_repl) {
    stop("invalid filter threshold applied", "\nRun filter_missval() with a threshold ranging from 0 to ",
         max_repl)
  }
  bin_data <- SummarizedExperiment::assay(se)
  idx <- is.na(SummarizedExperiment::assay(se))
  bin_data[!idx] <- 1
  bin_data[idx] <- 0
  keep <- bin_data %>% data.frame() %>% tibble::rownames_to_column() %>%
    gather(ID, value, -rowname) %>% left_join(., data.frame(SummarizedExperiment::colData(se)),
                                              by = "ID") %>% group_by(rowname, condition) %>% dplyr::summarize(miss_val = n() -
                                                                                                          sum(value)) %>% filter(miss_val <= thr) %>% tidyr::spread(condition,
                                                                                                                                                             miss_val)
  se_fltrd <- se[keep$rowname, ]
  return(se_fltrd)
}




#' Make Unique Proteins
#'
#' This function takes a data frame of proteins, a column name for names, a column name for IDs, and a delimiter and returns a data frame with unique proteins.
#'
#' @param proteins A data frame of protein abundances.
#' @param names A character string for the name of the column with protein names.
#' @param ids A character string for the name of the column with protein IDs (e.g., Accession numbers).
#' @param delim A character string for the delimiter.
#'
#' @return A data frame with unique proteins.
#'
#' @examples
#' prot.make_unique(proteins, names, ids, delim = ";")
#'
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom tibble is_tibble
#' @importFrom dplyr mutate
prot.make_unique <- function (proteins, names, ids, delim = ";")
{
  assertthat::assert_that(is.data.frame(proteins), is.character(names),
                          length(names) == 1, is.character(ids), length(ids) ==
                            1, is.character(delim), length(delim) == 1)
  #remove columns with empty header
  proteins <- proteins[, !(colnames(proteins)=="")]

  col_names <- colnames(proteins)
  if (!names %in% col_names) {
    stop("'", names, "' is not a column in '", deparse(substitute(proteins)),
         "'", call. = FALSE)
  }
  if (!ids %in% col_names) {
    stop("'", ids, "' is not a column in '", deparse(substitute(proteins)),
         "'", call. = FALSE)
  }
  if (tibble::is_tibble(proteins))
    proteins <- as.data.frame(proteins)
  proteins[, names][proteins[, names]=="NA"] <- NA
  proteins[, ids][proteins[, ids]=="NA"] <- NA
  double_NAs <- apply(proteins[, c(names, ids)], 1, function(x) all(is.na(x)))
  if (any(double_NAs)) {
    answer_satisfied <- "n"
    answer_satisfied <- readline("NAs in both the 'names' and 'ids' columns.\n Do you want to remove the respective protein entries (y/n)?\n\n")
    if ("n" %in% answer_satisfied) {
    stop()
    } else {
      proteins <- proteins[!double_NAs, ]
    }
  }
  proteins_unique <- proteins %>% dplyr::mutate(name = gsub(paste0(delim,
                                                            ".*"), "", get(names)), ID = gsub(paste0(delim, ".*"),
                                                                                              "", get(ids)), name = make.unique(ifelse(name == "" |
                                                                                                                                         is.na(name), ID, name)))
  return(proteins_unique)
}


#' Parse a data frame into a SummarizedExperiment object
#'
#' This function takes a data frame containing proteins unique to each sample and the corresponding expression values, and parses it into a SummarizedExperiment object.
#'
#' @param proteins_unique A data frame containing proteins unique to each sample. Must contain 'name' and 'ID' columns.
#' @param columns The columns containing the expression values. Must be numeric.
#' @param mode A character string specifying the format of the sample labels in the columns (default is "char").
#' @param chars The number of characters to consider when parsing the sample labels (default is 1).
#' @param sep The delimiter used to separate the sample labels (default is "_").
#'
#' @return A SummarizedExperiment object.
#' @noRd
#' @keywords internal
#' @importFrom SummarizedExperiment SummarizedExperiment
make_se_parse <- function (proteins_unique, columns, mode = c("char", "delim"),
          chars = 1, sep = "_")
{
  assertthat::assert_that(is.data.frame(proteins_unique), is.integer(columns),
                          is.character(mode), is.numeric(chars), length(chars) ==
                            1, is.character(sep), length(sep) == 1)
  mode <- match.arg(mode)
  if (any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)), "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if (any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric", "\nRun make_se_parse() with the appropriate columns as argument",
         call. = FALSE)
  }
  if (tibble::is_tibble(proteins_unique))
    proteins_unique <- as.data.frame(proteins_unique)
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)
  colnames(raw) <- delete_prefix(colnames(raw)) %>% make.names()
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$name
  if (mode == "char") {
    col_data <- data.frame(label = colnames(raw), stringsAsFactors = FALSE) %>%
      dplyr::mutate(condition = substr(label, 1, nchar(label) -
                                  chars), replicate = substr(label, nchar(label) +
                                                               1 - chars, nchar(label))) %>% unite(ID, condition,
                                                                                                   replicate, remove = FALSE)
  }
  if (mode == "delim") {
    col_data <- data.frame(label = colnames(raw), stringsAsFactors = FALSE) %>%
      separate(label, c("condition", "replicate"), sep = sep,
               remove = FALSE, extra = "merge") %>% unite(ID,
                                                          condition, replicate, remove = FALSE)
  }
  rownames(col_data) <- col_data$ID
  colnames(raw)[match(col_data$label, colnames(raw))] <- col_data$ID
  raw <- raw[, !is.na(colnames(raw))]
  se <- SummarizedExperiment:SummarizedExperiment(assays = as.matrix(raw), colData = col_data,
                             rowData = row_data)
  return(se)
}

#' Imputation by random draws from a manually defined distribution
#'
#' This function imputes missing values in a proteomics dataset by random draws from a manually defined distribution.
#'
#' @param se SummarizedExperiment, Proteomics data (output from make_se() or make_se_parse()). It is adviced to first remove proteins with too many missing values using filter_missval() and normalize the data using normalize_vsn().
#' @param scale Numeric(1), Sets the width of the distribution relative to the standard deviation of the original distribution.
#' @param shift Numeric(1), Sets the left-shift of the distribution (in standard deviations) from the median of the original distribution.
#'
#' @return A SummarizedExperiment object with imputed values
#'
#' @export
#'
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr filter group_by summarise
#' @importFrom tidyr gather
#' @importFrom stats rnorm
#'
#' @importFrom SummarizedExperiment assay
#'
#' @export manual_impute
manual_impute <- function (se, scale = 0.3, shift = 1.8)
{
  if (is.integer(scale))
    scale <- is.numeric(scale)
  if (is.integer(shift))
    shift <- is.numeric(shift)
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.numeric(scale), length(scale) == 1, is.numeric(shift),
                          length(shift) == 1)
  se_assay <- SummarizedExperiment::assay(se)
  if (!any(is.na(se_assay))) {
    stop("No missing values in '", deparse(substitute(se)),
         "'", call. = FALSE)
  }
  stat <- se_assay %>% data.frame() %>% tibble::rownames_to_column() %>%
    gather(samples, value, -rowname) %>% filter(!is.na(value)) %>%
    group_by(samples) %>% summarise(mean = mean(value), median = median(value),
                                    sd = sd(value), n = n(), infin = nrow(se_assay) - n)
  for (a in seq_len(nrow(stat))) {
    SummarizedExperiment::assay(se)[is.na(SummarizedExperiment::assay(se)[, stat$samples[a]]), stat$samples[a]] <- rnorm(stat$infin[a],
                                                                             mean = stat$median[a] - shift * stat$sd[a], sd = stat$sd[a] *
                                                                               scale)
  }
  return(se)
}
