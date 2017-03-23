suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(argparse)
})

parser <- ArgumentParser()

parser$add_argument("-n", "--inputnormal", type="character",
    default=file.path("data", "Gautier_Immgen_Norm_Data.tsv"),
    help="Path to input normal data")

parser$add_argument("-m", "--inputmeta", type="character",
    default=file.path("data", "Gautier_Immgen_Sample_Metadata.tsv"),
    help="Path to input meta data")

parser$add_argument("-o", "--output", type="character",
    default=file.path("results", "gene_data_vs_cell_type.tsv"),
    help="Path to output gentic data")

args <- parser$parse_args()

inputnormal_path <- args$inputnormal
# inputnormal_path <- file.path("data", "Gautier_Immgen_Norm_Data.tsv")
# inputnormal_path <- file.path("data", "Amit_Norm_Counts.tsv")
inputmeta_path <- args$inputmeta
# inputmeta_path <- file.path("data", "Gautier_Immgen_Sample_Metadata.tsv")
# inputmeta_path <- file.path("data", "Amit_Sample_Metadata.csv")
output_path <- args$output
# output_path <- file.path("results", "gene_data_vs_cell_type.tsv")
# output_path <- file.path("results", "gene_data_vs_cell_type.validation.tsv")


inputnormal_path %>%
  read_tsv(
    col_types = cols(
      Ensembl = col_character(),
      .default = col_double()
    )
  ) ->
  Norm_Data

Norm_Data %>%
  get("Ensembl", .) %>%
  c(., "GSM_ID", "Cell_Type", "General_Cell_Type")->
  gene_names

inputmeta_path %>%
  read_tsv(col_types = cols(
      GSM_ID = col_character(),
      Cell_Type = col_character(),
      General_Cell_Type = col_character()
    )
  ) ->
  Sample_Metadata

Norm_Data %>%
  t %>%
  as.data.frame %>%
  mutate(GSM_ID = colnames(Norm_Data)) %>%
  inner_join(Sample_Metadata) %>%
  `colnames<-`(gene_names) ->
  gene_data_vs_cell_type

write_tsv(gene_data_vs_cell_type, output_path)
