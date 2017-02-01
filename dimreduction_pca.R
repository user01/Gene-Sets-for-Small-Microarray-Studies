suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("-d", "--dimensions", type="integer", default=2,
    help="Number of dimensions to recover")

parser$add_argument("-n", "--name", type="character", default="standard",
    help="Name of data results")
args <- parser$parse_args()

gene_data_vs_cell_type <- file.path("results", "gene_data_vs_cell_type.tsv") %>%
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
# pca_results %>% glimpse
# plot(pca_results$x[,1:2])

path_target <- file.path("results",
  paste0("dimreduced_pca_",
         args$name,
         "_",
         args$dimensions,
         ".tsv"
         )
       )

pca_results$x[, 1:args$dimensions] %>%
  as.data.frame %>%
  cbind(gene_data_vs_cell_type %>%
  select(Cell_Type,General_Cell_Type)) %>%
  write_tsv(path_target)
