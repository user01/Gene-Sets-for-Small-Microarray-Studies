suppressPackageStartupMessages({
  library(MASS)
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
# top_genes <- ceiling(1/2 * (sqrt(8*200 + 1) + 1))
input_path <- args$input
# input_path <- file.path("results", "gene_data_vs_cell_type.tsv")
cell_name <- args$name
# cell_name <- "NK cell"
cell_type <- args$type
# cell_type <- "General_Cell_Type"
input_path <- args$input
# input_path <- file.path("results", "gene_data_vs_cell_type.tsv")
output_path <- args$output
# output_path <- file.path("results", "results.tsv")

# Read base data file and transform for type
input_path %>%
  read_tsv(col_types=cols(
    .default = col_double(),
    GSM_ID = col_character(),
    Cell_Type = col_character(),
    General_Cell_Type = col_character()
  )) %>%
  mutate_(.dots = setNames(cell_type, "type_truth")) %>%
  select(-GSM_ID, -Cell_Type, -General_Cell_Type) ->
  genes

genes %>% select(-type_truth) -> gene_data
genes %>%
  select(type_truth) %>%
  unlist %>%
  unname %>%
  map_chr(~ if (. != cell_name) { "other" } else { cell_name }) ->
  gene_labels

set.seed(random_seed)

# Bootstrap genes, generate gene indexes
gene_data %>%
  ncol %>%
  `*`(fraction_bootstrap) %>%
  floor ->
  genes_count
gene_data %>%
  ncol %>%
  sample(replace = TRUE) %>%
  head(genes_count) %>%
  select(gene_data, .) %>%
  cbind(data.frame(type_truth = gene_labels)) ->
  genes_bootstrapped

# Perform downsampling to match given set
genes_bootstrapped %>%
  filter(type_truth == cell_name) ->
  genes_bootstrapped_truth
genes_bootstrapped %>%
  filter(gene_labels != cell_name) %>%
  head(nrow(genes_bootstrapped_truth)) %>%
  rbind(genes_bootstrapped_truth) ->
  genes_bootstrapped_downsampled

# Seperate out for LDA
genes_bootstrapped_downsampled %>%
  select(-type_truth) ->
  genes_bootstrapped_downsampled_values
genes_bootstrapped_downsampled %>%
  select(type_truth) %>%
  unlist %>%
  unname %>%
  as.character ->
  genes_bootstrapped_downsampled_labels

# Perform fit for components
lda_model <- lda(genes_bootstrapped_downsampled_values, genes_bootstrapped_downsampled_labels, CV = FALSE)

# Select top loading genes from components
lda_model %>%
  get("scaling", .) %>%
  (function(m) {
    data.frame(
      loading = m %>% unlist %>% unname,
      gene = m %@% "dimnames" %>% nth(1)
    )
  }) %>%
  mutate(
    loading = abs(loading),
    gene = as.character(gene)
  ) %>%
  arrange(-loading) %>%
  head(top_genes) %>%
  get("gene", .) %>%
  unlist ->
  top_gene_names



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
