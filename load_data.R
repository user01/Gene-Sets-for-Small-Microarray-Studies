library(readr)

EMF_GSVA_Gene_Set_Enrichment_PATH <- file.path("data", "EMF_GSVA_Gene_Set_Enrichment.tsv.tar.gz")
EMF_GSVA_Gene_Set_Enrichment <- read_tsv(EMF_GSVA_Gene_Set_Enrichment_PATH, col_types = cols(.default = col_integer(),
  EMF_GSVA_Gene_Set_Enrichment.tsv = col_character(), Gene_Set = col_character()))

# EMF_GSVA_Gene_Set_Enrichment
