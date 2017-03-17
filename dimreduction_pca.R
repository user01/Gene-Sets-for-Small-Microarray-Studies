suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("-d", "--dimensions", type="integer", default=2,
    help="Number of dimensions to recover")

parser$add_argument("-i", "--input", type="character",
    default=file.path("results", "gene_data_vs_cell_type.tsv"),
    help="Path to input data")

parser$add_argument("-o", "--output", type="character",
    default=file.path("results", "results_pca.tsv"),
    help="Path to output reduced dimension data")

args <- parser$parse_args()

dimensions <- args$dimensions
# dimensions <- 2
input_path <- args$input
# input_path <- file.path("results", "gene_data_vs_cell_type.tsv")
output_path <- args$output
# output_path <- file.path("results", "results_pca.tsv")


gene_data_vs_cell_type <- input_path %>%
  read_tsv(col_types=cols(
    .default = col_double(),
    GSM_ID = col_character(),
    Cell_Type = col_character(),
    General_Cell_Type = col_character()
  ))

expression_data <- gene_data_vs_cell_type %>%
                   select(-GSM_ID,-Cell_Type,-General_Cell_Type) %>%
                   as.matrix

pca_results <- expression_data %>%
               scale(center = TRUE, scale = FALSE) %>%
               t %>%
               scale(center = FALSE, scale = TRUE) %>%
               t %>%
               prcomp(center = FALSE)


pca_results$x[, 1:args$dimensions] %>%
  as.data.frame %>%
  cbind(gene_data_vs_cell_type %>%
  select(Cell_Type,General_Cell_Type)) %>%
  write_tsv(output_path)
