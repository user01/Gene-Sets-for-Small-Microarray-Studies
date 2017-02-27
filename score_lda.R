suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("-s", "--seed", type="integer", default=451,
    help="Random seed")

parser$add_argument("-p", "--pairs", type="integer", default=200,
    help="Number of gene pairs to pick as relevant from a component")

parser$add_argument("-b", "--bootstrap", type="double", default=0.66,
    help="Fraction of genes to bootstrap")

parser$add_argument("-n", "--name", type="character",
    help="Name of cell type to target")

parser$add_argument("-t", "--type", type="character",
    help="General or specific cell type")


parser$add_argument("-i", "--input", type="character",
    default=file.path("results", "gene_data_vs_cell_type.tsv"),
    help="Path to input data")

parser$add_argument("-o", "--output", type="character",
    default=file.path("results", "results_pairs.tsv"),
    help="Path to output pair/scores data")

args <- parser$parse_args()

random_seed <- args$seed
# random_seed <- 451
fraction_bootstrap <- args$bootstrap
# fraction_bootstrap <- 0.66
top_genes <- ceiling(1/2 * (sqrt(8*args$pairs + 1) + 1)) # this uses nC2 to pick the top set to create ~number of pairs
# top_genes <- ceiling(1/2 * (sqrt(8*2000 + 1) + 1))
input_path <- args$input
# input_path <- file.path("results", "gene_data_vs_cell_type.tsv")
cell_name <- args$name
# cell_name <- "NK cell"
cell_type <- args$type
# cell_name <- "general"
input_path <- args$input
# input_path <- file.path("results", "gene_data_vs_cell_type.tsv")
output_path <- args$output
# output_path <- file.path("results", "results.tsv")



# genes <- input_path %>%
#   read_tsv(col_types=cols(
#     .default = col_double(),
#     GSM_ID = col_character(),
#     Cell_Type = col_character(),
#     General_Cell_Type = col_character()
#   ))
#
# genes %>% select(-GSM_ID, -Cell_Type, -General_Cell_Type) -> gene_data
# genes %>% select(Cell_Type, General_Cell_Type) -> gene_labels

set.seed(random_seed)

# TODO: Fill in LDA code

place_holder <- data.frame(low=c('a','b','c'),high=c('b','c','d'),frequency=1:3,score=1:3)

place_holder %>%
  write_tsv(output_path)
