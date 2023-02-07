#' Reads flux data from a file or data frame
#'
#' @param file.data Path to the file containing the flux data. Alternatively, a data frame can be provided.
#' @param csvsep Separator used in the file. Default is ";".
#' @param dec Decimal separator used in the file. Default is ".".
#' @param sheet.data Sheet number to read from the file. Default is 1.
#' @param file.match Path to the file containing the equation-BiGG name matches.
#' @param sheet.match Sheet number to read from the file matching equations to BiGG reaction names. Default is 1.
#' @param rename.reaction Logical, whether to rename reactions according to BiGG standard. Default is TRUE.
#' @param rescale.reaction Name of a reaction to normalize all fluxes (will be set as 100%). Default is NULL.
#' @param col.id Name of the column containing reaction IDs. Default is "Reaction ID".
#' @param col.eq Name of the column containing carbon transition equations. Default is "Equation".
#' @param col.flux Name of the column containing flux values. Default is "Value (%)".
#' @param col.sd Name of the column containing standard deviation values. Default is "Standard Deviation".
#' @param FBA Logical, whether the data is from an FBA simulation. Default is FALSE.
#' @return A data frame with reaction IDs, equations, flux values, and standard deviation.
#' @export
#'
read_flux <- function(file.data = dat_flux,
                      csvsep = ";",
                      dec = ".",
                      sheet.data = 1,
                      file.match= "eq_bigg.csv",
                      sheet.match = 1,
                      rename.reaction = TRUE,
                      rescale.reaction = NULL,
                      col.id="Reaction ID",
                      col.eq = "Equation",
                      col.flux= "Value (%)",
                      col.sd= "Standard Deviation",
                      FBA = FALSE){

  # Read file that matches equations to BiGG reaction names
  if (rename.reaction == TRUE) {
    if (!is.character(file.match)) {
      eq_bigg <- file.match
    } else {
      # Read table file
      if (file.exists(file.match)) {
        eq_bigg <- read_file(file.match, csvsep=csvsep, dec=dec, sheet=sheet.match)
      } else {
        stop(
          paste0(
            "File \"",
            file.match,
            "\" does not exist.\n Please provide a table file to match equations to BiGG names"
          ),
          call. = F
        )
      }
    }
  }
  # if(is.null(rescale.reaction)) stop("Please provide the name of a reaction to normalize all fluxes (will be set as 100%).")

  # Read data file
  if (!is.character(file.data)) {
    df_data <- file.data
  } else {
    # Read table file
    if (file.exists(file.data)) {
      df_data <- read_file(file.data, csvsep = csvsep, dec = dec, sheet = sheet.data)
    } else {
      stop(paste0("File \"", file.data, "\" does not exist."), call. = F)
    }
  }
  if(!(col.id %in% colnames(df_data))){
    stop(paste0("A column with name ", col.id,
                " does not exist in the provided data file.\nPlease define the header of the column containing reaction IDs with 'col.id ='\n",
                "Valid column names are: ", paste(colnames(df_data), collapse=", ")))
  }
  if(FBA == FALSE){
    if(!(col.eq %in% colnames(df_data))){
      stop(paste0("A column with name ", col.eq,
                  " does not exist in the provided data file.\nPlease define the header of the column containing carbon transition equations with 'col.eq ='\n",
                  "Valid column names are: ", paste(colnames(df_data), collapse=", ")))
    }
  }
  if(!(col.flux %in% colnames(df_data))){
    stop(paste0("A column with name ", col.flux,
                " does not exist in the provided data file.\nPlease define the header of the column containing flux values with 'col.flux ='\n",
                "Valid column names are: ", paste(colnames(df_data), collapse=", ")))
  }
  if(FBA == FALSE){
    if(!(col.sd %in% colnames(df_data))){
      stop(paste0("A column with name ", col.sd,
                  " does not exist in the provided data file.\nPlease define the header of the column containing standard deviations with 'col.sd ='\n",
                  "Valid column names are: ", paste(colnames(df_data), collapse=", ")))
    }
  }
  # Rename reactions according to BiGG standard based on matches of the equations in the reference file.
  if (FBA == FALSE) {
    df_data <- as.data.frame(cbind(reaction= df_data[,col.id], equation = df_data[,col.eq], flux = df_data[,col.flux], SD = df_data[,col.sd]))
    for (i in 1:nrow(df_data)) {
      df_data$reaction[i] <-
        eq_bigg$reaction[match(df_data$equation[i], eq_bigg$equation)]
    }
  } else {
    df_data <- as.data.frame(cbind(reaction= df_data[,col.id], flux = df_data[,col.flux]))
    df_data$SD <- 0
    df_data$equation <- 0
    if(!is.null(rescale.reaction)){
      # Scale all fluxes except the growth rate
      df_rescaled <- df_data
      df_rescaled$flux <- as.numeric(df_data$flux) / as.numeric(df_data[match(rescale.reaction, df_data$reaction), "flux"]) *100
      # Restore original growth rate value
      df_rescaled[grep("BIOMASS_KT2440_WT3", df_rescaled$reaction), "flux"] <-  df_data[grep("BIOMASS_KT2440_WT3", df_data$reaction), "flux"]
      df_data <- df_rescaled
    }

    if (rename.reaction == TRUE) {
      # Filter reactions in results from FBA to only contain reactions listed in the match file
      df_data <- df_data %>% filter(reaction %in% (paste(eq_bigg$reaction)) )
    }
  }

  # convert columns with standard deviation and flux values to numeric
  df_data <- transform(df_data, SD=as.numeric(SD))
  df_data <- transform(df_data, flux=as.numeric(flux))

  # while being treated separately in INCA, the two reaction pairs OAADC/PC and PPC/PPCK are lumped together to visualize the flux
  # between Pyr <-> OAA or PEP <-> OAA, respectively. The standard deviation is calculated via Gaussian error propagation

  # Calculate the difference between the reaction OAADC and PC, then rename the reaction
  df_PC <- df_data[match("PC", df_data$reaction),]
  df_OAADC <- df_data[match("OAADC", df_data$reaction),]
  df_Pyr_OAA <- rbind(df_PC, df_OAADC) %>%
    mutate(flux = -diff(flux))
  # Apply Gaussian Error propagation to StdDev values
  df_Pyr_OAA <- df_Pyr_OAA %>%
    mutate(SD = sqrt(sum((SD)^2)))
  # Renaming
  df_Pyr_OAA <- df_Pyr_OAA[-c(2,0),]
  df_Pyr_OAA["reaction"] <- "Pyr_OAA"
  df_Pyr_OAA["equation"] <- "net flux from Pyr to OAA"
  # difference between the reaction PPCK and PPC
  df_PPC <- df_data[match("PPC", df_data$reaction),]
  df_PPCK <- df_data[match("PPCK", df_data$reaction),]
  df_PEP_OAA <- rbind(df_PPC, df_PPCK) %>%
    mutate(flux = -diff(flux))
  # apply Gaussian Error propagation to StdDev values
  df_PEP_OAA <- df_PEP_OAA %>%
    mutate(SD = sqrt(sum((SD)^2)))
  # Renaming
  df_PEP_OAA <- df_PEP_OAA[-c(2,0),]
  df_PEP_OAA["reaction"] <- "PEP_OAA"

  df_PEP_OAA["equation"] <- "net flux from PEP to OAA"
  # append Pyr <-> OAA and PEP <-> OAA reactions to dataframe
  df_data <- rbind(df_data, df_Pyr_OAA)
  df_data <- rbind(df_data, df_PEP_OAA)

  # remove the initial rows for the reactions OAADC, PC, PPC, and PPCK
  df_data <- df_data[- grep("OAADC|PC|PPC|PPCK", df_data$reaction),]

  # Adjust sign of flux for 2-ketogluconate and gluconate exchange reactions based on their equation
  df_data[match("GLCNtex", df_data$reaction), "flux"] <-
    if_else(df_data[match("GLCNtex", df_data$reaction), "equation"] == "Gluco.per (abcdef) -> Gluco.ext (abcdef)",
            df_data[match("GLCNtex", df_data$reaction), "flux"],
            -df_data[match("GLCNtex", df_data$reaction), "flux"])

  df_data[match("2DHGLCNtex", df_data$reaction), "flux"] <-
    if_else(df_data[match("2DHGLCNtex", df_data$reaction), "equation"] == "Kgluco.per (abcdef) -> Kgluco.ext (abcdef)",
            df_data[match("2DHGLCNtex", df_data$reaction), "flux"],
            -df_data[match("2DHGLCNtex", df_data$reaction), "flux"])

  return(df_data)
}

#' Create a metabolic map with fluxes and export as SVG file
#'
#' \code{flux_to_map} takes a table file or dataframe as well as a template map in SVG format to visualize metabolic flux values in their metabolic context.
#'
#' @param df A dataframe containing fluxes and standard deviations.
#' @param template A template SVG file containing the metabolic map.
#' @param result.nm The name of the output SVG file.
#' @param pal The color palette for the flux arrows. G, green; R, red; O, orange; BYR, blue-yellow-red; BW, black-white (grayscale); YR, yellow-red; PuRd, purple-red.
#' @param title The title of the output plot.
#' @param export A logical value indicating whether to export the plot as a PNG and PDF file.
#' @param inkscape_path The local path to the 'inkscape.exe' file. Required if \code{export = TRUE}.
#' @param export_dpi The dpi of the exported PNG and PDF file.
#' @param export_width The width of the exported PNG and PDF file. The width of the template SVG can be inspected after opening it in InkScape.
#' @param export_height The height of the exported PNG and PDF file. The height of the template SVG can be inspected after opening it in InkScape.
#' @param FBA A logical value indicating whether the dataframe contains FBA results. If TRUE, standard deviations are ignored.
#' @return A metabolic map with fluxes.
#' @export
#'
#' @importFrom scales rescale
#' @importFrom stringr str_replace
#'
flux_to_map <- function (df = "dat_flux",
                         template = NULL,
                         result.nm = NULL,
                         pal = c("G", "R", "O", "B", "BYR", "BW", "YR", "PuRd"),
                         title = "",
                         export = TRUE,
                         inkscape_path = 'C:/Program Files/Inkscape/bin/inkscape.exe',
                         export_dpi = 300,
                         export_width = 2281,
                         export_height = 2166,
                         FBA = FALSE)
{
  pal <- match.arg(pal)
  dir.create(paste0(getwd(),"/Plots"), showWarnings = F)

  # Define palettes for flux arrows and legend
  if (pal == "G") {
    colflux_0   <- '#d9fcba' # light green for a flux value of 0%
      colflux_50 <-  '#83c087' # medium-light green for a flux value of 50%
        colflux_100 <- '#2e8555' # full green for a flux value of 100%
          colflux_max <- '#1d5254' # dark green for fluxes > 100%
  } else if (pal == "R") {
    colflux_0   <- '#fee0d2' # light red for a flux value of 0%
      colflux_50 <-  '#fc9272' # medium-light red for a flux value of 50%
        colflux_100 <- '#de2d26' # full red for a flux value of 100%
          colflux_max <- '#67000d' # dark red for fluxes > 100%
  } else if (pal == "O") {
    colflux_0   <- '#fee6ce' # light red for a flux value of 0%
      colflux_50 <-  '#fd8d3c' # medium-light red for a flux value of 50%
        colflux_100 <- '#d94801' # full red for a flux value of 100%
          colflux_max <- '#7f2704' # dark red for fluxes > 100%
  } else if (pal == "B") {
    colflux_0   <- '#cfdff5' # light blue for a flux value of 0%
      colflux_50 <-  '#6eadcd' # medium-light blue for a flux value of 50%
        colflux_100 <- '#116aab' # full blue for a flux value of 100%
          colflux_max <- '#004678' # dark blue for fluxes > 100%
  } else if (pal == "BYR") {
    colflux_0   <- '#2b83ba' # blue for a flux value of 0%
      colflux_50 <-  '#FFE135' # yellow for a flux value of 50%
        colflux_100 <- '#d7191c' # red for a flux value of 100%
          colflux_max <- '#a50f15' # dark red for fluxes > 100%
  } else if (pal == "gray" | pal ==  "grey") {
    colflux_0   <- '#f0f0f0' # light gray for a flux value of 0%
      colflux_50 <-  '#bdbdbd' # medium-gray for a flux value of 50%
        colflux_100 <- '#636363' # full gray for a flux value of 100%
          colflux_max <- '#252525' # dark gray for fluxes > 100%
  } else if (pal == "BW") {
    colflux_0   <- '#000000' # black for a flux value of 0%
      colflux_50 <-  '#000000' # black for a flux value of 50%
        colflux_100 <- '#000000' # black for a flux value of 100%
          colflux_max <- '#000000' # black for fluxes > 100%
  } else if (pal == "YR") {
    colflux_0   <- '#fff7bc' # black for a flux value of 0%
      colflux_50 <-  '#fec44f' # black for a flux value of 50%
        colflux_100 <- '#d95f0e' # black for a flux value of 100%
          colflux_max <- '#993404' # black for fluxes > 100%
  } else if (pal == "PuRd"){
    colflux_0   <- '#f1eef6' # black for a flux value of 0%
      colflux_50 <-  '#df65b0' # black for a flux value of 50%
        colflux_100 <- '#ce1256' # black for a flux value of 100%
          colflux_max <- '#91003f' # black for fluxes > 100%
  } else {
    stop(
      "No valid color palette provided to draw flux arrows.
      Please define a suitable palette (G, R, B, BYR, gray, YR, or BW) as \"pal =\" argument.",
      call. = FALSE
    )
  }

  # Create color palettes for flux values
  pal_0to100 <- colorRampPalette(c(colflux_0, colflux_50, colflux_100))(201)
  pal_above100 <- colorRampPalette(c(colflux_100,colflux_max))(201)

  # Define palette for Flux legend
  legend_elements <- c(sprintf("rect_flux%s", seq(1:67)))
  pal_legend_gradient <- colorRampPalette(pal_0to100)(67)
  df_fluxlegend <- data.frame(legend_elements,pal_legend_gradient)

  # Store growth rate value (rounded to two decimal digits)
  biomass <- df[grep("BIOMASS", df[,"reaction"]), grep("flux", colnames(df))] %>% round(digits = 2) %>% format(nsmall = 2)
  #round flux and standard deviation values to no decimal digits
  df$flux <- round(df$flux, digits = 0)
  df$SD <- round(df$SD, digits = 0)

  # Scaling and coloring of flux arrows according to flux values.
  # Create dataset only with reactions that have a flux unequal 0
  metabolic_flux <- df %>%
    filter(flux != 0)
  # Separate reactions with flux greater than 100 to ensure a uniform color coding up to 100
  metabolic_flux_0to100 <- metabolic_flux %>% filter(abs(flux) <= 100)
  metabolic_flux_above100 <- metabolic_flux %>% filter( abs(flux) > 100)

  # add three rows corresponding to flux values of 0.01, 50, and 100 to have a uniform color distribution
  df_colref <- data.frame(matrix(NA, nrow = 3, ncol = 4))
  colnames(df_colref) <- colnames(df)
  df_colref$reaction <-  c("flux0", "flux50", "flux100")
  df_colref$equation <- c("lower reference for color palette", "mid reference for color palette", "upper reference for color palette")
  df_colref$SD <- c(0, 0, 0)
  df_colref$flux <- c(0.01, 50, 100)
  metabolic_flux_0to100 <- rbind(metabolic_flux_0to100, df_colref)

  # take the absolute value of fluxes
  df$flux_abs <- abs(df$flux)

  # concatenate flux values and standard deviations into new column (excluding the value for the growth rate)
  df$concat <- paste(df$flux_abs, "\u00B1", df$SD)
  if(FBA == FALSE){
    df$concat <- gsub(" ", "", paste(df$concat)) #this line removes the spaces before and after "+-"
  } else {
    df$concat <- gsub(" \u00B1 .+", "", paste(df$concat)) #this line the standard deviation (not applicable to FBA results)

  }
  # remove values for reactions with zero flux
  df <- df %>% dplyr::mutate(concat = ifelse(flux == 0, " ", concat))

  # As some reactions have very high flux and others no flux at all, we apply a square root transformation to the fluxes so that stroke width is more balanced.
  # add two columns for stroke width and color
  metabolic_flux_0to100 <- metabolic_flux_0to100 %>%
    mutate(
      stroke_width = 1 + 0.45*sqrt(abs(as.numeric(flux))),
      stroke_color = (abs(as.numeric(flux))) %>%
        scales::rescale(to = c(1, 200)) %>% round,
      stroke_color_rgb = pal_0to100[stroke_color])
  # reactions with flux >100 is assigned a dedicated color palette
  metabolic_flux_above100 <- metabolic_flux_above100 %>%
    mutate(
      stroke_width = 1 + 0.45*sqrt(abs(flux)),
      stroke_color = sqrt(abs(flux)) %>%
        scales::rescale(to = c(1, 200)) %>% round,
      stroke_color_rgb = pal_above100[stroke_color])
  # combine the two dataframes with formatted reactions
  metabolic_flux_formatted <- merge(metabolic_flux_0to100, metabolic_flux_above100, all = TRUE)

  message("Reading map template...")

  # read SVG map template
  if (!is.null(template) && is.character(template)) {
    template.nm <- template %>% str_replace(".svg", "")
    template_svg <- paste0(template.nm, ".svg")
    if (!file.exists(template)) {
      stop(paste0("The template file \"", template_svg, "\" does not exist."), .call = FALSE)
    } else {
      MAP <- fluctuator::read_svg(template_svg)
    }
  }
  else{
    stop(
      "No template SVG provided. Templates can be specified with \"template =\".",
      call. = FALSE
    )
  }

  if (title != ""){
    # change the Title of the image
    message("Changing map title...")
    MAP <- fluctuator::set_values(MAP,
                                  node = paste0("Title"),
                                  value = title)
  }

  # get names of reactions in map
  rxn.nm <- MAP@summary$label[!is.na(MAP@summary$label)][!grepl("_text|rect|scalebar", MAP@summary$label[!is.na(MAP@summary$label)])]
  rxn.nm <- gsub("value_", "", rxn.nm)
  # reduce df to contain only entries present in map (for quicker computation)
  df <- df[df[["reaction"]] %in% rxn.nm, ]

  message("Applying flux values to text fields...")

  MAP <- fluctuator::set_values(MAP,
                                node = paste0("value_", df$reaction),
                                value = df$concat)

  message("Applying growth rate value to text field...")
  if(!is.na(charmatch("BIOMASS", df[,"reaction"]))){
    MAP <- fluctuator::set_values(MAP,
                                  node = paste0("value_BIOMASS_KT2440_WT3"),
                                  value = biomass)
  } else {
    cat("The reaction 'BIOMASS' was not found in the data table. The growth rate will be removed from the map.\n")
    MAP <- fluctuator::set_values(MAP,
                                  node = paste0("value_BIOMASS_KT2440_WT3"),
                                  value = "")
    MAP <- fluctuator::set_values(MAP,
                                  node = paste0("?_text"),
                                  value = "")
  }

  message("Applying stroke widths to arrows...")
  MAP <- fluctuator::set_attributes(MAP,
                                    node = metabolic_flux_formatted$reaction,
                                    attr = "style",
                                    pattern = "stroke-width:[0-9]+\\.[0-9]+",
                                    replacement = paste0("stroke-width:", metabolic_flux_formatted$stroke_width))

  message("Applying stroke colors to arrows...")
  MAP <- fluctuator::set_attributes(MAP,
                                    node = metabolic_flux_formatted$reaction,
                                    attr = "style",
                                    pattern = "stroke:#b2b2b2",
                                    replacement = paste0("stroke:", metabolic_flux_formatted$stroke_color_rgb))

  message("Apply color scheme to flux legend...")
  MAP <- fluctuator::set_attributes(MAP,
                                    node = df_fluxlegend$legend_elements,
                                    attr = "style",
                                    pattern = "fill:#b3b3b3",
                                    replacement = paste0("fill:", df_fluxlegend$pal_legend_gradient))

  # adjust arrow head size in the map
  message("Adjusting arrow head sizes...")
  MAP <- fluctuator::set_attributes(MAP,
                                    node = grep("marker", MAP@summary$id, value = TRUE),
                                    node_attr = "id",
                                    attr = "transform",
                                    pattern = "scale\\(0.2\\)",
                                    replacement = "scale(0.3)")

  MAP <- fluctuator::set_attributes(MAP,
                                    node = grep("marker", MAP@summary$id, value = TRUE),
                                    node_attr = "id",
                                    attr = "transform",
                                    pattern = "scale\\(-0.2\\)",
                                    replacement = "scale(-0.1)")

  # display arrows with flux 0.0 (rounded to one decimal digit) as grey, dashed lines
  metabolic_flux_zero <- df %>%
    filter(flux == 0)

  # apply stroke width=1.5 to arrows with flux 0.0
  message("Changing appearance of arrows with zero flux ...")
  MAP <- fluctuator::set_attributes(MAP,
                                    node = metabolic_flux_zero$reaction,
                                    attr = "style",
                                    pattern = "stroke-width:[0-9]+\\.[0-9]+",
                                    replacement = paste0("stroke-width:", "1.5"))

  # apply stroke style  (dashed) to arrows in the map
  MAP <- fluctuator::set_attributes(MAP,
                                    node = metabolic_flux_zero$reaction,
                                    attr = "style",
                                    pattern = "stroke-dasharray:none",
                                    replacement = paste0("stroke-dasharray:2"))


  # create empty arrowheads based on the directionality of the reaction
  #define dataframe with negative flux values
  metabolic_flux_neg <- df %>%
    filter(flux < 0)
  #replace end arrow heads for negative fluxes with empty heads
  message("Replace end arrow heads for negative fluxes with empty heads...")
  MAP <- fluctuator::set_attributes(MAP,
                                    node = metabolic_flux_neg$reaction,
                                    attr = "style",
                                    pattern = "marker-end:url\\(#TriangleOutS.*\\);",
                                    replacement = "marker-end:url\\(#EmptyTriangleOutS\\);")
  #define dataframe with positive flux values
  metabolic_flux_pos <- df %>%
    filter(flux > 0)
  #replace start arrow heads for positive fluxes with empty heads
  message("Replace start arrow heads for positive fluxes with empty heads...")
  MAP <- fluctuator::set_attributes(MAP,
                                    node = metabolic_flux_pos$reaction,
                                    attr = "style",
                                    pattern = "marker-start:url\\(#TriangleInS+\\-[0-9]+",
                                    replacement = "marker-start:url\\(#EmptyTriangleInS")

  df[grep("BIOMASS", df$reaction), grep("flux", colnames(df))] <- biomass

  # Write modified SVG map

  if (!is.null(result.nm)){
    result.nm <- result.nm %>% str_replace(".svg", "")
    result_svg <- paste0("Plots/", result.nm, ".svg")
    message(paste0("Writing results SVG file with defined name: \"",
                   result_svg,"\" ...") )
  } else {
    result.nm <- str_replace(paste0("Plots/", template.nm), "template", "result")
    message(paste0("Writing results SVG file:\n'Plots/",
                   result.nm, "+F",".svg' ..."))
    result_svg <- paste0(result.nm,"+F",".svg")
  }

  fluctuator::write_svg(MAP, file = result_svg)

  if(export == TRUE){
    svg_to_png(svg = result_svg, width=export_width, inkscape_path = inkscape_path, height=export_height, dpi=export_dpi)
    svg_to_pdf(svg = result_svg, width=export_width, inkscape_path = inkscape_path, height=export_height, dpi=export_dpi)
  }

  invisible(df)
}

#' Title
#'
#' @param se \code{SummarizedExperiment} object, proteomics data parsed with \code{\link{prot.read_data}}.
#' @param protein_set Data frame containing a list of proteins to be plotted. The provided names in column \code{} must match the entries in the column provided as \code{col.id} as well as the IDs of objects in the map template.
#' @param col.id (Character) Column name of the protein identifier in the SummarizedExperiment object.
#' @param select.id Column name of the protein identifier in the \code{protein_set} data frame.
#' @param pfx Prefix of the protein identifiers to be plotted. all protein objects in the template map should have a common prefix (e.g., 'prot_')
#' @param imp_fun (Character string)  Function used for data imputation. "SampMin", "man", "bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb", "min", "zero", "mixed", or "nbavg". See (\code{\link{prot.impute}}) for details.
#' @param q (Numeric) q value for imputing missing values with method \code{imp_fun = 'MinProb'}.
#' @param type (Character string) Type of analysis to perform. "contrast" (provide a specific contrast in the \code{contrast} argument to plot log2-fold change values) or "centered" (provide a single condition as \code{condition} to plot centered abundances in that condition).
#' @param condition  Condition to be displayed if type = "centered".
#' @param contrast Tested contrast if type = "contrast" in the format "A_vs_B".
#' @param template Path to the SVG map template.
#' @param RSD_threshRelative standard deviation threshold.
#' @param p_thresh (adjusted) p-value threshold.
#' @param result Name of the output files.
#' @param pal Brewer color palette used for protein squares.
#' @param legend_min Lower limit of the color legend.
#' @param legend_max Upper limit of the color legend.
#' @param title (Character string) Title shown above the map
#' @param export (Logical) A logical value indicating whether to export the plot as a PNG and PDF file.
#' @param inkscape_path The local path to the 'inkscape.exe' file. Required if \code{export = TRUE}.
#' @param export_dpi The dpi of the exported PNG and PDF file.
#' @param export_width The width of the exported PNG and PDF file. The width of the template SVG can be inspected after opening it in InkScape.
#' @param export_height The height of the exported PNG and PDF file. The height of the template SVG can be inspected after opening it in InkScape.
#'
#' @return A metabolic map with protein objects colored according to their abundance or log2 fold changes in a given contrast.
#' @export
#'
prot_to_map <- function (se,
                         protein_set = NULL,
                         col.id = 'ID',
                         select.id = 'Accession',
                         pfx = NULL,
                         imp_fun = c("man", "bpca", "knn", "QRILC", "MLE", "MinDet", # Method for imputing of missing values
                                     "MinProb", "min", "zero", "mixed", "nbavg", "SampMin"),
                         q = 0.01, # q value for imputing missing values with method "fun = 'MinProb'".
                         type = c("contrast", "centered"),
                         condition = NULL, # Condition to be displayed if type = "centered".
                         contrast = NULL, # Tested contrast if type = "contrast" in the format "A_vs_B"
                         template = NULL,
                         RSD_thresh = 0.5,
                         p_thresh = 1,
                         indicate_nonsig = c("gray", "dashed"),
                         result = NULL,
                         pal = "RdYlBu",
                         legend_min = NULL,
                         legend_max = NULL,
                         title = "",
                         export = TRUE,
                         inkscape_path = 'C:/Program Files/Inkscape/bin/inkscape.exe',
                         export_dpi = 300,
                         export_width = 2281,
                         export_height = 2166)
{
  type <- match.arg(type)
  imp_fun <- match.arg(imp_fun)
  indicate_nonsig <- match.arg(indicate_nonsig)
  dir.create(paste0(getwd(),"/Plots"), showWarnings = F)
  # Define palettes for protein squares and legend
  if (pal %in% rownames(RColorBrewer::brewer.pal.info)) {
    pal_prot <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =pal)))(201) # Spectral Palette
  } else {
    stop(
      paste0(
        "No valid color palette provided to indicate protein log2(fold change) values.\nPlease define any brewer palette as \"pal =\" argument:\n",
        "'",
        paste(rownames(RColorBrewer::brewer.pal.info), collapse = "', '"),
        "'"
      ),
      call. = FALSE
    )
  }
  valid_conditions <- se$condition %>% unique() %>% paste(., collapse=", ")
  if (type == "contrast" & is.null(contrast)){
    valid_contrasts <- se$condition %>% unique() %>% paste(., collapse=", ")
    stop(paste0("Type 'contrast' was chosen without providing a valid contrast.\n Please provide a contrast in the format \"A_vs_B\" with two of the conditions: '", valid_conditions, "'"),
         call. = FALSE)
  }
  # Define palette for Proteomics legend
  legend_elements <- c(sprintf("rect_prot%s", seq(1:67)))
  pal_legend_gradient <- colorRampPalette(pal_prot)(67)
  df_protlegend <- data.frame(legend_elements,pal_legend_gradient)

  # Variance stabilization
  prot_norm <- suppressMessages(prot.normalize_vsn(se, plot = FALSE, export = FALSE))
  # Impute missing values
  prot_imp <- prot.impute(prot_norm, fun = imp_fun, q = q)

  message("Reading map template...")

  # read SVG map template
  if (!is.null(template)) {
    template.nm <- template %>% str_replace(".svg", "")
    template_svg <- paste0(template.nm, ".svg")
  } else {
    stop( "No template SVG file provided. Templates can be addressed by specifying \"template =\".",
          call. = FALSE )
  }
  if (!(file.exists(template_svg))){
    stop( "The template file '", template_svg, "' does not exist.\n Please provide a valid template SVG file in the 'template = ' argument.",
          call. = FALSE )
  }

  MAP <- fluctuator::read_svg(template_svg)

  # Create list of proteins in the map based on matches with the defined prefix (pfx)
  if (is.null(protein_set) && !is.null(pfx)){
    assign("protein_set", stats::setNames(data.frame(MAP@summary$label[grep( paste(pfx, collapse="|"), MAP@summary$label)]), paste(select.id)))
  }

  # create a vector for the proteins that were not present in the data set
  selected_ids <- protein_set[,select.id]
  detected  <- SummarizedExperiment::rowData(prot_imp) %>% data.frame(check.names = F) %>% pull(col.id)
  not_detected <- setdiff(selected_ids, detected)

  if (type == "centered" && is.null(condition)){
    stop(paste0( "Type 'centered' was chosen without providing a valid condition.\n Valid conditions are: '", valid_conditions, "'" ),
         call. = FALSE)
  }
  if (type == "centered"){
    if(!(condition %in% se$condition)){
      stop(paste0( "The selected 'condition' is not present in the dataset.\n Valid conditions are: '", valid_conditions, "'" ),
           call. = FALSE)
    }
    df <- data.frame(assay(prot_imp), check.names = F)
    row_data <- data.frame(SummarizedExperiment::rowData(prot_imp), check.names = F)
    df$ID <- row_data[[col.id]]

    #Calculate the relative standard deviation for the chosen condition
    SD <-  df %>% select(contains(condition)) %>% `^`(2, .) %>% as.matrix() %>% rowSds()
    mean <- df %>% select(contains(condition)) %>% `^`(2, .) %>% rowMeans(na.rm = TRUE)
    df$RSD <-  if_else( is.na(SD/mean), 0, (SD/mean) )
    # Center the data
    log2mean <- df %>% select(-(ID:RSD)) %>% rowMeans(na.rm = TRUE)
    df <- df %>% select(-(ID:RSD)) %>% - log2mean %>% cbind(., ID = df$ID, RSD = df$RSD)
    # Calculate the centered average of the chosen condition
    df$avg <- df %>% select(contains(condition)) %>% rowMeans()
    # filter dataframe for defined list of genes (in argument "protein_set")
    df <- filter(df, str_detect(ID, paste(selected_ids, collapse="|")))
    # Filter for proteins with a relative standard deviation below "RSD_thresh = "
    df_sig <- df %>% filter(RSD < RSD_thresh)

    filtered <- setdiff(df$ID, df_sig$ID)
    if(length(filtered)>0){
      cat(paste0("The following proteins exceeded the RSD threshold of 'RSD_thresh =",
                 RSD_thresh,
                 "' and/or the p-value threshold of 'p_thresh = ",
                 p_thresh,
                 ":\n",
                 paste(filtered, collapse = ", "),"\n"))
    }
    df_filtered <- filter(df, str_detect(ID, paste(filtered, collapse="|"))) %>% select(ID, values = avg, RSD)

    # Create data frame with protein IDs and plotted values
    df_plot <- cbind(ID = df_sig$ID, values = df_sig$avg, RSD = df_sig$RSD) %>% as.data.frame(check.names = FALSE)
  }

  if (type == "contrast") {
    # Test for differential expression by empirical Bayes moderation
    # of a linear model and defined contrasts
    prot_diff <- prot.test_diff(prot_imp, type = "manual", test = contrast)
    row_data <- data.frame(SummarizedExperiment::rowData(prot_diff), check.names = F)
    df <- data.frame(assay(prot_diff), check.names = F)
    df$ID <- row_data[[col.id]]
    #add log2(fold change) values from prot_diff
    df$diff <- row_data[,grep("_diff", colnames(row_data))]
    #add p-values values from prot_diff
    df$p.val <- row_data[,grep("_p.adj", colnames(row_data))]
    # Calculate the relative standard deviation for subject and reference conditions
    subject <- str_replace(contrast, "_vs.+", "")
    ref <- str_replace(contrast, ".+vs_", "")
    SD_subject <-  df %>% select(contains(subject)) %>% `^`(2, .) %>% as.matrix() %>% rowSds()
    mean_subject <- df %>% select(contains(subject)) %>% `^`(2, .) %>% rowMeans(na.rm = TRUE)
    SD_ref <-  df %>% select(contains(ref)) %>% `^`(2, .) %>% as.matrix() %>% rowSds()
    mean_ref <- df %>% select(contains(ref)) %>% `^`(2, .) %>% rowMeans(na.rm = TRUE)
    df$RSD_subject <-  SD_subject/mean_subject
    df$RSD_ref <-  if_else(is.na(SD_ref), 0 ,SD_ref)/mean_ref

    # Filter dataframe for defined list of genes (argument "protein_set")
    df <- filter(df, str_detect(ID, paste(selected_ids, collapse="|")))

    # Filter for proteins with a relative standard deviation below "RSD_thresh = " or a p-value < 0.05
    if(indicate_nonsig == "gray"){
      df_sig <- df %>% filter( (RSD_subject <= RSD_thresh & RSD_ref <= RSD_thresh) & (p.val <= p_thresh) )
    }
    else{
      df_sig <- df %>% filter( (RSD_subject <= RSD_thresh & RSD_ref <= RSD_thresh) )
      df_sig_pval <- df %>% filter( (p.val <= p_thresh) )
      filtered_pval <- setdiff(df$ID, df_sig_pval$ID)
      df_filtered_pval <- filter(df, str_detect(ID, paste(filtered_pval, collapse="|"))) %>% select(ID, values = diff,
                                                                                                    RSD_ref, RSD_subject, p.val)
    }
    filtered <- setdiff(df$ID, df_sig$ID)


    if(length(filtered)>0){
      cat(paste0("The following proteins exceeded the RSD threshold of 'RSD_thresh = ",
                 RSD_thresh,
                 "' and/or the p-value threshold of 'p_thresh = ",
                 p_thresh,
                 ":\n",
                 paste(filtered, collapse = ", "), "\n"))

      df_filtered <- filter(df, str_detect(ID, paste(filtered, collapse="|"))) %>% select(ID, values = diff,
                                                                                          RSD_ref, RSD_subject, p.val)
    } else {
      df_filtered <- data.frame()
    }

    # Create data frame with protein IDs and plotted values
    df_plot <- cbind(ID = df_sig$ID, values = as.numeric(df_sig$diff), RSD_ref = as.numeric(df_sig$RSD_ref),
                     RSD_subject = as.numeric(df_sig$RSD_subject), p.val = as.numeric(df_sig$p.val)) %>% as.data.frame(check.names = FALSE)
  }

  # add three rows corresponding to log2(fold change) values "min", 0, and "max" to have a uniform color distribution and
  # adjust the proteomics legend in the figure. min = -max, with their value corresponding to the highest abs(log2fc) observed
  if (!is.null(legend_min)){
    legend_min_value <- legend_min
  } else {
    legend_min_value <- -ceiling(max(abs(as.numeric(df_plot$values))))
  }
  if (!is.null(legend_max)){
    legend_max_value <- legend_max
  } else {
    legend_max_value <- ceiling(max(abs(as.numeric(df_plot$values))))
  }
  prot_0_value <-   0
  colref_min <- c("colref_min",     legend_min_value) # lower reference for color palette
  colref_zero <-   c("colref_zero",    prot_0_value) # mid reference for color palette
  colref_max <- c("colref_max",     legend_max_value)  # upper reference for color palette
  # create dataframe with legend boundaries
  prot_legend_names <- c(colref_min[1],colref_zero[1],colref_max[1])
  prot_legend_values <- c(colref_min[2],colref_zero[2],colref_max[2])
  prot_legend <- data.frame( prot_legend_names, prot_legend_values, check.names = FALSE )
  # add legend boundaries as rows to the dataframe with proteomics data
  if(type == "centered"){
    df_colbound <- data.frame( ID = prot_legend_names, values = prot_legend_values,
                               RSD = 0, check.names = FALSE)
  } else if (type == "contrast"){
    df_colbound <- data.frame( ID = prot_legend_names, values = prot_legend_values,
                               RSD_ref = 0, RSD_subject = 0, p.val = 0, check.names = FALSE)
  }
  df_plot <- rbind.data.frame(df_plot, df_colbound)
  df_plot$values <- as.numeric(df_plot$values)
  # add columns for colors
  df_plot_within_colbounds <-
    df_plot[which(df_plot$values >= legend_min_value & df_plot$values <= legend_max_value), ] %>%
    mutate(stroke_color = values %>%
             scales::rescale(to = c(1, 201)) %>% round,
           fill_color_rgb = pal_prot[stroke_color])

  df_plot_outside_colbounds <-
    df_plot[which(df_plot$values < legend_min_value | df_plot$values > legend_max_value), ] %>%
    mutate(stroke_color = if_else(values < legend_min_value, 1, 201),
           fill_color_rgb = pal_prot[stroke_color])

  df_plot <- bind_rows(df_plot_within_colbounds, df_plot_outside_colbounds)
  # remove rows for the legend to color protein boxes
  df_plot <- df_plot[-(grep("colref", df_plot$ID)),]

  if (title != ""){
    message("Changing map title...")
    MAP <- fluctuator::set_values(MAP,
                                  node = paste0("Title"),
                                  value = title)
  }


  message("Apply fill colors to protein boxes in the map...")
  MAP <- fluctuator::set_attributes(MAP,
                                    node = df_plot$ID, attr = "style",
                                    pattern = "fill:#858585",
                                    replacement = paste0("fill:", df_plot$fill_color_rgb))

  if(indicate_nonsig == "dashed"){
    message("Make borders of protein boxes dashed for nonsignificant changes...")
    MAP <- fluctuator::set_attributes(MAP,
                                      node = filtered_pval, attr = "style",
                                      pattern = "stroke-dasharray:none",
                                      replacement = "stroke-dasharray:0.568075,1.704225")
  }


  message("Apply min and max values to the proteomics legend...")
  MAP <- fluctuator::set_values(MAP,
                                node = prot_legend$prot_legend_names,
                                value = prot_legend$prot_legend_values)

  if (type == "centered"){
    # Change legend to "Log2 centered intensity"
    message("Change legend to 'Log2 (centered intensity)'...")
    MAP <- fluctuator::set_values(MAP,
                                  node = "prot_type",
                                  value = "Log2(centered intensity)")
  }

  message("Apply color scheme to protein legend...")
  MAP <- fluctuator::set_attributes(MAP,
                                    node = df_protlegend$legend_elements, attr = "style",
                                    pattern = "fill:#b3b3b3",
                                    replacement = paste0("fill:", df_protlegend$pal_legend_gradient))

  message("Make boxes for proteins not found in the dataset black...")
  MAP <- fluctuator::set_attributes(MAP,
                                    node = not_detected, attr = "style",
                                    pattern = "fill:#858585",
                                    replacement = "fill:#000000")

  # Write modified SVG map
  result.nm <- str_replace(template.nm, "template", "result")

  if (!is.null(result)){
    result.nm <- result %>% str_replace(".svg", "")
    result_svg <- paste0("Plots/", result.nm, ".svg")
    message(paste0("Writing results SVG file with defined name: \"",
                   result_svg,"\" ...") )
  } else {
    message(paste0("Writing results SVG file: \n'Plots/",
                   result.nm, "+P",".svg' ..."))
    result_svg <- paste0("Plots/", result.nm,"+P",".svg")
  }

  fluctuator::write_svg(MAP, file = result_svg)
  if ( export == TRUE ){
    svg_to_png(svg = result_svg, width=export_width, inkscape_path = inkscape_path, height=export_height, dpi=export_dpi)
    svg_to_pdf(svg = result_svg, width=export_width, inkscape_path = inkscape_path, height=export_height, dpi=export_dpi)
  }
  if(length(filtered)>0){
    df_filtered$stroke_color <- ""
    df_filtered$fill_color_rgb <- ""
  }

  if(type == "centered"){
    df_plot <- rbind.data.frame(df_plot, df_filtered) %>% mutate(values = as.numeric(values),
                                                                 RSD = as.numeric(RSD))
  } else if (type == "contrast"){
    df_plot <- rbind.data.frame(df_plot, df_filtered) %>% mutate(values = as.numeric(values),
                                                                 RSD_ref = as.numeric(RSD_ref),
                                                                 RSD_subject = as.numeric(RSD_subject),
                                                                 p.val = as.numeric(p.val))
  }
  df_plot <- df_plot[(df_plot$stroke_color != ""),]
  invisible(df_plot)
}


#' @title Convert SVG to PNG (requires Python and InkScape installed)
#'
#' @description
#' This function converts an SVG file to a PNG file by calling the respective InkScape function through Python.
#'
#' @param svg A character string specifying the path of the SVG file to be converted.
#' @param inkscape_path The local path to the 'inkscape.exe' file.
#' @param out A character string specifying the desired name of the output file. If `NULL` (the default), the output file will have the same name as the input file, but with a .png extension.
#' @param dpi A numeric value specifying the desired dots per inch of the output file. The default is 300.
#' @param width A numeric value specifying the desired width of the output file. The default is 2281.
#' @param height A numeric value specifying the desired height of the output file. The default is 2166.
#'
#' @return A message specifying the location of the output file.
#'
#' @export
#'
#' @import stringr
#'
#' @importFrom reticulate py_run_string
#'
#'
#' @keywords internal
svg_to_png <-
  function (svg = svg,
            inkscape_path = 'C:/Program Files/Inkscape/bin/inkscape.exe',
            out = NULL,
            dpi = 300,
            width = 2281,
            height = 2166) {
    dir.create(paste0(getwd(),"/Plots"), showWarnings = F)
    if (file.exists(svg)){
      svg <- svg
    } else if (file.exists(paste0(getwd(),"/Plots/",svg))){
      svg <- paste0(getwd(),"/Plots/",svg)
    } else {
      stop("File \"", svg, "\" not found in specified location or in ",getwd(), "/Plots.
           Please provide a complete file path or assure that file exists.", call. = FALSE)
    }

    if (is.null(out)) {
      png.nm <- str_replace_all(svg,".{1,}/","") %>% str_replace(.,".svg", "")
      png <- paste0(getwd(),"/Plots/", png.nm, ".png")
    }
    else {
      png = paste0(getwd(),"/Plots/", out)
    }

    width <- width/300*dpi
    height <- height/300*dpi
    py_run_string(paste0('import subprocess
inkscape = ', '"', inkscape_path, '"'))

    py_run_string(stringr::str_glue('subprocess.run([inkscape, "--export-type=png",
  f"--export-filename={png}",
  f"--export-width={width}",
  f"--export-height={height}",
  "{svg}"
  ])'))

    message(paste0("Exporting PNG file to: ", png))
  }

#' Convert SVG to PDF
#'
#' @param svg character. File path to the svg file.
#' @param inkscape_path The local path to the 'inkscape.exe' file.
#' @param out character. File path to the output pdf file.
#' @param dpi numeric. DPI of the output pdf file.
#' @param width numeric. Width of the output pdf file.
#' @param height numeric. Height of the output pdf file.
#'
#' @export
#'
#' @return A message indicating the file path of the output pdf file.
#'
#' @importFrom stringr str_glue
#'
#' @export
#'
svg_to_pdf <- function (svg = svg,
                        inkscape_path = 'C:/Program Files/Inkscape/bin/inkscape.exe',
                        out = NULL,
                        dpi = 300,
                        width = 2281,
                        height = 2166)
{
    dir.create(paste0(getwd(),"/Plots"), showWarnings = F)
    if (file.exists(svg)){
      svg <- svg
    } else if (file.exists(paste0(getwd(),"/Plots/",svg))){
      svg <- paste0(getwd(),"/Plots/",svg)
    } else {
      stop("File \"", svg, "\" not found in specified location or in ",getwd(), "/Plots.
           Please provide a complete file path or assure that file exists.", call. = FALSE)
    }

    if (is.null(out)) {
      pdf.nm <- str_replace_all(svg,".{1,}/","") %>% str_replace(.,".svg", "")
      pdf <- paste0(getwd(),"/Plots/", pdf.nm, ".pdf")
    }
    else {
      pdf = paste0(getwd(),"/Plots/", out)
    }

    width <- width/300*dpi
    height <- height/300*dpi
    py_run_string(paste0('import subprocess
inkscape = ', '"', inkscape_path, '"'))



    py_run_string(stringr::str_glue('subprocess.run([inkscape, "--export-type=pdf",
  f"--export-filename={pdf}",
  f"--export-width={width}",
  f"--export-height={height}",
  "{svg}"
  ])'))

    message(paste0("Exporting PDF file to: ", pdf ))
  }


#' @title Convert PNGs to an animated GIF
#' @description Converts a vector of PNG files to an animated GIF.
#'
#' @param png Path to a PNG file.
#' @param ... Paths to additional PNG files.
#' @param fps Frames per second for the animated GIF.
#' @param out Path to the output file.
#'
#' @return A list of \code{\link{magick}} objects.
#'
#' @examples
#' \dontrun{
#' png_files <- list.files(path = "Plots/",
#'                         pattern = "*.png",
#'                         full.names = TRUE)
#' png_to_gif(png = png_files,
#'            fps = 1.0,
#'            out = "Plots/animated_gif.gif")
#' }
#'
#' @export
#'
#' @import magick
png_to_gif <-
  function(png, ..., fps = 0.5, out = "Plots/animated_gif.gif")
  {
    call <- match.call()
    call$out <- NULL
    call$fps <- NULL
    arglist <- lapply(call[-1], function(x)
      x)
    var.names <- vapply(arglist, deparse, character(1))
    arglist <- lapply(arglist, eval.parent, n = 2)
    names(arglist) <- var.names
    png_list <- list()
    for (i in 1:length(arglist)) {
      png_list[[i]] <- magick::image_read(arglist[[i]])
    }
    img_joined <- magick::image_join(unlist(png_list))
    img_animated <- magick::image_animate(img_joined, fps = fps)
    magick::image_write(image = img_animated,
                        path = out)
  }
