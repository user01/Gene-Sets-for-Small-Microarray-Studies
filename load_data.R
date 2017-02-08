suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})


"EMF_GSVA_Relative_Normalilzed_Expression.tsv.tar.gz" %>%
  file.path("data", .) %>%
  read_tsv(
    col_types = cols(
      .default = col_double(),
      EMF_GSVA_Relative_Normalilzed_Expression.tsv = col_character(),
      Symbol = col_character(),
      Entrez = col_integer(),
      Name = col_character()
    )
  ) %>%
  filter(!is.na(Symbol)) ->
  EMF_GSVA_Relative_Normalilzed_Expression

EMF_GSVA_Relative_Normalilzed_Expression %>%
  get("EMF_GSVA_Relative_Normalilzed_Expression.tsv", .) %>%
  c(., "GSM_ID", "Cell_Type", "General_Cell_Type")->
  gene_names

"Gautier_Immgen_Sample_Metadata.tsv" %>%
  file.path("data", .) %>%
  read_tsv(col_types = cols(
      GSM_ID = col_character(),
      Cell_Type = col_character(),
      General_Cell_Type = col_character()
    )
  ) ->
  Gautier_Immgen_Sample_Metadata

EMF_GSVA_Relative_Normalilzed_Expression %>%
  t %>%
  as.data.frame %>%
  mutate(GSM_ID = colnames(EMF_GSVA_Relative_Normalilzed_Expression)) %>%
  inner_join(Gautier_Immgen_Sample_Metadata) %>%
  `colnames<-`(gene_names) ->
  gene_data_vs_cell_type

write_tsv(gene_data_vs_cell_type, file.path("results", "gene_data_vs_cell_type.tsv"))
