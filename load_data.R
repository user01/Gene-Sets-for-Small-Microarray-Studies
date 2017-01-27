suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

EMF_GSVA_Relative_Normalilzed_Expression_PATH <- file.path("data", "EMF_GSVA_Relative_Normalilzed_Expression.tsv.tar.gz")
EMF_GSVA_Relative_Normalilzed_Expression <- read_tsv(EMF_GSVA_Relative_Normalilzed_Expression_PATH,
  col_types = cols(.default = col_double(), EMF_GSVA_Relative_Normalilzed_Expression.tsv = col_character(),
    Symbol = col_character(), Entrez = col_integer(), Name = col_character())) %>%
  filter(!is.na(Symbol)) %>% select(-1:-4) %>% as.matrix

write_tsv(EMF_GSVA_Relative_Normalilzed_Expression, file.path("results", "EMF_GSVA_Relative_Normalilzed_Expression_MATRIX"), col_names = FALSE)
