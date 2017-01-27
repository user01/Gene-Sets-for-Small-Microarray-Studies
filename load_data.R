library(readr)

EMF_GSVA_Gene_Set_Enrichment_PATH <- file.path("data", "EMF_GSVA_Gene_Set_Enrichment.tsv.tar.gz")
EMF_GSVA_Gene_Set_Enrichment <- read_tsv(EMF_GSVA_Gene_Set_Enrichment_PATH, col_types = cols(.default = col_integer(),
  EMF_GSVA_Gene_Set_Enrichment.tsv = col_character(), Gene_Set = col_character()))
# EMF_GSVA_Gene_Set_Enrichment

EMF_GSVA_Gene_Expression_PATH <- file.path("data", "EMF_GSVA_Relative_Normalilzed_Expression.tsv.tar.gz")
EMF_GSVA_Gene_Expression <- read_tsv(EMF_GSVA_Gene_Expression_PATH, col_types = cols(.default = col_character(),
  EMF_GSVA_Relative_Normalilzed_Expression.tsv = col_character(), Gene_Set = col_character()))
EMF_GSVA_Gene_Expression <- EMF_GSVA_Gene_Expression[-20271,]
#EMF_GSVA_Gene_Expression

numcols <- c(5:141)
EMF_GSVA_Gene_Expression[, numcols] = apply(EMF_GSVA_Gene_Expression[,numcols], 2, 
  function(x) as.numeric(as.character(x)))
#Changed gene expression columns to numeric for EMF_GSVA_Gene_Expression

