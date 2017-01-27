suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
})

EMF_GSVA_Gene_Set_Enrichment_PATH <- file.path("data", "EMF_GSVA_Gene_Set_Enrichment.tsv.tar.gz")
EMF_GSVA_Gene_Set_Enrichment <- read_tsv(EMF_GSVA_Gene_Set_Enrichment_PATH, col_types = cols(.default = col_integer(),
  EMF_GSVA_Gene_Set_Enrichment.tsv = col_character(), Gene_Set = col_character()))

EMF_GSVA_Relative_Normalilzed_Expression_PATH <- file.path("data", "EMF_GSVA_Relative_Normalilzed_Expression.tsv.tar.gz")
EMF_GSVA_Relative_Normalilzed_Expression <- read_tsv(EMF_GSVA_Relative_Normalilzed_Expression_PATH,
  col_types = cols(.default = col_double(), EMF_GSVA_Relative_Normalilzed_Expression.tsv = col_character(),
    Symbol = col_character(), Entrez = col_integer(), Name = col_character()))

Gautier_Immgen_Sample_Metadata_PATH <- file.path("data", "Gautier_Immgen_Sample_Metadata.tsv")
Gautier_Immgen_Sample_Metadata <- read_tsv(Gautier_Immgen_Sample_Metadata_PATH, col_types = cols(GSM_ID = col_character(),
  Cell_Type = col_character(), General_Cell_Type = col_character()))

Gautier_Immgen_Standardized_Norm_Data_PATH <- file.path("data", "Gautier_Immgen_Standardized_Norm_Data.tsv.tar.gz")
Gautier_Immgen_Standardized_Norm_Data <- read_tsv(Gautier_Immgen_Standardized_Norm_Data_PATH,
  col_types = cols(.default = col_double(), Gautier_Immgen_Standardized_Norm_Data.tsv = col_character(),
    Symbol = col_character(), Entrez = col_integer(), Name = col_character()))

EMF_GSVA_Gene_Set_Enrichment
EMF_GSVA_Relative_Normalilzed_Expression
# EMF_GSVA_Relative_Normalilzed_Expression$GSM854326 %>% sd(na.rm=T)
Gautier_Immgen_Sample_Metadata
Gautier_Immgen_Standardized_Norm_Data
