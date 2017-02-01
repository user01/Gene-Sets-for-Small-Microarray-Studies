suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})


EMF_GSVA_Relative_Normalilzed_Expression_PATH <- file.path("data", "EMF_GSVA_Relative_Normalilzed_Expression.tsv.tar.gz")
EMF_GSVA_Relative_Normalilzed_Expression <- read_tsv(EMF_GSVA_Relative_Normalilzed_Expression_PATH,
  col_types = cols(.default = col_double(), EMF_GSVA_Relative_Normalilzed_Expression.tsv = col_character(),
    Symbol = col_character(), Entrez = col_integer(), Name = col_character())) %>%
  filter(!is.na(Symbol)) %>% select(-1:-4)

Gautier_Immgen_Sample_Metadata_PATH <- file.path("data", "Gautier_Immgen_Sample_Metadata.tsv")
Gautier_Immgen_Sample_Metadata <- read_tsv(Gautier_Immgen_Sample_Metadata_PATH, col_types = cols(GSM_ID = col_character(),
  Cell_Type = col_character(), General_Cell_Type = col_character()))

gene_data_vs_cell_type <- EMF_GSVA_Relative_Normalilzed_Expression %>% t %>% as.data.frame %>%
  mutate(GSM_ID = colnames(EMF_GSVA_Relative_Normalilzed_Expression)) %>% inner_join(Gautier_Immgen_Sample_Metadata)

write_tsv(gene_data_vs_cell_type, file.path("results", "gene_data_vs_cell_type.tsv"))

