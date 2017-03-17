suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})


"Gautier_Immgen_Norm_Data.tsv" %>%
  file.path("data", .) %>%
  read_tsv(
    col_types = cols(
      Ensembl = col_character(),
      .default = col_double()
    )
  ) ->
  Gautier_Immgen_Norm_Data

Gautier_Immgen_Norm_Data %>%
  get("Ensembl", .) %>%
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

Gautier_Immgen_Norm_Data %>%
  t %>%
  as.data.frame %>%
  mutate(GSM_ID = colnames(Gautier_Immgen_Norm_Data)) %>%
  inner_join(Gautier_Immgen_Sample_Metadata) %>%
  `colnames<-`(gene_names) ->
  unscaled_gene_data_vs_cell_type

unscaled_gene_data_vs_cell_type %>%
  select(-GSM_ID, -Cell_Type, -General_Cell_Type) %>%
  apply(2,function(x) as.numeric(as.character(x))) %>%
  scale %>%
  cbind(unscaled_gene_data_vs_cell_type %>%
          select(GSM_ID, Cell_Type, General_Cell_Type)) ->
  gene_data_vs_cell_type

write_tsv(gene_data_vs_cell_type, file.path("results", "gene_data_vs_cell_type.tsv"))
