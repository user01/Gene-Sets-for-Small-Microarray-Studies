suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(argparse)
  library(Rtsne)
})


parser <- ArgumentParser()
parser$add_argument("-s", "--seed", type="integer", default=451,
    help="Random seed value")

parser$add_argument("-p", "--perplexity", type="integer", default=30,
    help="T-SNE perplexity number")

parser$add_argument("-m", "--max_iter", type="integer", default=1000,
    help="Maximum number of iterations")

parser$add_argument("-c", "--pca", type="character", default="false",
    help="Use PCA in first pass")

parser$add_argument("-i", "--input", type="character",
    default=file.path("results", "gene_data_vs_cell_type.tsv"),
    help="Path to input data")

parser$add_argument("-o", "--output", type="character",
    default=file.path("results", "results_tsne.tsv"),
    help="Path to output reduced dimension data")

args <- parser$parse_args()


seed <- args$seed
# seed <- 451
perplexity <- args$perplexity
# perplexity <- 2
max_iter <- args$max_iter
# max_iter <- 2
pca <- args$pca == "true"
# pca <- TRUE
input_path <- args$input
# input_path <- file.path("results", "gene_data_vs_cell_type.tsv")
output_path <- args$output
# output_path <- file.path("results", "results_tsne.tsv")


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

set.seed(seed) # Sets seed for reproducibility
tsne_out <- Rtsne(expression_data,
                  perplexity = perplexity,
                  max_iter = max_iter,
                  pca = pca)

tsne_out %>%
  get("Y", .) %>%
  as.data.frame %>%
  cbind(gene_data_vs_cell_type %>%
  select(Cell_Type,General_Cell_Type)) %>%
  write_tsv(output_path)
