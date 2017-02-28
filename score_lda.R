suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(argparse)
  library(MASS)
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

# random_seed <- args$seed
random_seed <- 451
# fraction_bootstrap <- args$bootstrap
fraction_bootstrap <- 0.66
# top_genes <- ceiling(1/2 * (sqrt(8*args$pairs + 1) + 1)) # this uses nC2 to pick the top set to create ~number of pairs
top_genes <- ceiling(1/2 * (sqrt(8*2000 + 1) + 1))
# input_path <- args$input
input_path <- file.path("results", "gene_data_vs_cell_type.tsv")
# cell_name <- args$name
cell_name <- "NK cell"
cell_name <- "Microglia"
cell_name <- "Macrophage"
cell_name <- "Monocyte"
# cell_type <- args$type
cell_type <- "general"
# input_path <- args$input
input_path <- file.path("results", "gene_data_vs_cell_type.tsv")
# output_path <- args$output
output_path <- file.path("results", "results.tsv")

if (cell_type == "general") {
  label_col <- "General_Cell_Type"
} else {
  label_col <- "Cell_Type"
}

genes <- input_path %>%
  read_tsv(col_types=cols(
    .default = col_double(),
    GSM_ID = col_character(),
    Cell_Type = col_character(),
    General_Cell_Type = col_character()
  ))

genes %>% dplyr:: select(-GSM_ID, -Cell_Type, -General_Cell_Type) -> gene_data
genes %>% dplyr:: select(Cell_Type, General_Cell_Type) -> gene_labels

set.seed(random_seed)

# bootstrap genes, generate gene indexes 
gene_data_bootstraped <- gene_data[sample(1: floor(ncol(gene_data) * fraction_bootstrap))]

# balance the number of control vs. target sample
target_data <- gene_data_bootstraped[gene_labels[[label_col]] == cell_name ,] 

control_data <- gene_data_bootstraped[gene_labels[[label_col]] != cell_name ,] %>%
  sample_n(nrow(target_data)) 

training_data <- rbind(target_data, control_data)
training_labels <- c(rep(cell_name, nrow(target_data)), rep("others", nrow(target_data)))

# train LDA model 
lda.model <- lda(training_data, training_labels,CV=FALSE) 

# LDA prediction score 
gene_data_bootstraped %>%
  predict(lda.model, newdata = .) -> prediction

prediction_result <- prediction$class

expected_result <- gene_labels[[label_col]]
expected_result[expected_result != cell_name] = "others"

accuracy <- sum(prediction_result == expected_result)/length(expected_result)

# accuracy
# TODO: sort LDA loadings/scalings, combine into pairs score with lda.predict and freq = 1

# NOTE: loadings/scaling are stored in lda.model$scaling

place_holder <- data.frame(
  low=c('a','b','c'),
  high=c('b','c','d'),
  cell_name = cell_name,
  cell_type = cell_type,
  frequency=1:3,
  score=1:3)

place_holder %>%
  write_tsv(output_path)
