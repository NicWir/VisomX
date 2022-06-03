library(tidyr)

setwd("C:/Users/nicwir/Documents/DTU_Biosustain/Scripts_and_Modelling/fluctuator/220111/R_package/VisomX")
devtools::load_all()
setwd("C:/R_test/transcriptomics")

# Add old locus_id
reference.df <- read.csv("KT2440_gene_references.csv")
dataset <- read.csv("RNASeq_raw_complete_names_NICO.csv", sep = ";")
dataset$'locus_tag' <- reference.df$old_locus_tag[match(dataset$gene_id, reference.df$locus_tag)]

se <- prot.read_data(dataset, pfx = "abundance.", id = "locus_tag", name = "SymbolID")


custom_df <- read.table("C:/Users/nicwir/OneDrive - Danmarks Tekniske Universitet/PhD DTU/Scripts and Modelling/fluctuator/220111/proteomics analysis/All-pathways-of-P.-putida-KT2440.txt",
                        sep = "\t", header = T, stringsAsFactors = F, fill = T, na.strings = "", quote = "",
                        comment.char = "",  check.names = F)

results <- prot.workflow(se,
                         imp_fun = "SampMin",
                         type = "control",
                         control = "Glc",
                         # contrast = c(
                         #   "pS6311.PHB_glucose_vs_pS631.PHB_glucose",
                         #   "pS6311.PHB_acetate_vs_pS631.PHB_acetate",
                         #   "pS6310.PHB_acetate_vs_pS631.PHB_acetate"
                         # ),
                         alpha = 0.05,
                         lfc = 2,
                         heatmap.kmeans = T,
                         k = 6,
                         heatmap.show_row_names = T,
                         heatmap.row_font_size = 4.5,
                         heatmap.show_all = F,
                         volcano.adjusted = F,
                         plot = F,
                         export = T,
                         volcano.add_names = T,
                         report = F,
                         pathway_enrichment = T,
                         pathway_kegg = T,
                         kegg_organism = "ppu",
                         custom_pathways = custom_df
)
