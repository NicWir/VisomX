## ----library, echo = T, eval = F, results = "hide"----------------------------
# # Load the VisomX package
# library(VisomX)

## ----libraryCode, echo = F, eval = T, include = F, results = "hide"-----------
# Load the VisomX package
suppressPackageStartupMessages(library(VisomX))

## ----SamplePath, echo = T, eval = T, include = T, results = "asis"------------
system.file("prot_KT2440_glc_ac_LB.txt", package = "VisomX")


## ----read, include = T, eval = T----------------------------------------------
# Read data, incl. filtering
data <- prot.read_data(
          data = system.file("prot_KT2440_glc_ac_LB.txt", package = "VisomX"),
          pfx = "abundance.",
          id = "Ensembl Gene ID",
          name = "Gene Symbol",
          filt_type = "fraction",
          filt_min = 0.66)

