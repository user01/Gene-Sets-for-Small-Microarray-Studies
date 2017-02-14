suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("-s", "--seed", type="integer", default=451,
    help="Random seed")

parser$add_argument("-t", "--toppairs", type="integer", default=200,
    help="Number of gene pairs to pick as relevant from a component")

parser$add_argument("-b", "--bootstrap", type="double", default=0.66,
    help="Fraction of genes to bootstrap")

parser$add_argument("-f", "--fraction", type="double", default=0.5,
    help="Fraction of variation to include in relevant components")


parser$add_argument("-i", "--input", type="character",
    default=file.path("results", "gene_data_vs_cell_type.tsv"),
    help="Path to input data")

parser$add_argument("-c", "--outputcomponents", type="character",
    default=file.path("results", "results_components.tsv"),
    help="Path to output components data")

parser$add_argument("-p", "--outputpairs", type="character",
    default=file.path("results", "results_pairs.tsv"),
    help="Path to output pairs data")

args <- parser$parse_args()

random_seed <- args$seed
# random_seed <- 451
fraction_bootstrap <- args$bootstrap
# fraction_bootstrap <- 0.66
fraction_variation <- args$fraction
# fraction_variation <- 0.5
top_genes <- ceiling(1/2 * (sqrt(8*args$toppairs + 1) + 1)) # this uses nC2 to pick the top set to create ~number of pairs
# top_genes <- ceiling(1/2 * (sqrt(8*2000 + 1) + 1))
input_path <- args$input
# input_path <- file.path("results", "gene_data_vs_cell_type.tsv")
output_components_path <- args$outputcomponents
# output_components_path <- file.path("results", "results_components.tsv")
output_pairs_path <- args$outputpairs
# output_pairs_path <- file.path("results", "results_pairs.tsv")



genes <- input_path %>%
  read_tsv(col_types=cols(
    .default = col_double(),
    GSM_ID = col_character(),
    Cell_Type = col_character(),
    General_Cell_Type = col_character()
  ))

genes %>% select(-GSM_ID, -Cell_Type, -General_Cell_Type) -> gene_data
genes %>% select(Cell_Type, General_Cell_Type) -> gene_labels


set.seed(random_seed)
gene_data %>%
  ncol %>%
  `*`(fraction_bootstrap) %>%
  floor %>%
  sample.int(replace = TRUE) ->
  selected_genes_idx

gene_data %>%
  select(selected_genes_idx) %>%
  prcomp(scale = TRUE) ->
  pca_results

pca_results %>%
  get('sdev', .) %>%
  `^`(2) %>%
  (function(x) {x / sum(x) }) %>%
  cumsum %>% # Cumulative Proportion of Variance Explained
  discard(~ . > fraction_variation) %>%
  length ->
  pca_relevant_components_length

pca_results$x[, 1:pca_relevant_components_length] %>%
  as.data.frame %>%
  cbind(gene_labels) %>%
  write_tsv(output_components_path)

pca_results %>%
  get('rotation', .) %>%
  abs %>%
  (function(aload) { sweep(aload, 2, colSums(aload), "/") }) %>% # proportional contribution to the each principal component
  (function(m) { # Pick top top_genes genes from each
    m %>%
      colnames %>%
      head(pca_relevant_components_length) ->
      cols
    map(cols, function(c) {
      m[,c] %>%
        sort(decreasing = TRUE) %>%
        head(top_genes) %>%
        names %>%
        data.frame
    }) %>%
    reduce(cbind) %>%
    `colnames<-`(cols)
  }) ->
  genes_relevant


pair_ <- function(v) {
  if (length(v) < 2) {
    return(data.frame())
  }
  v %>% head(1) -> h
  v %>% tail(-1) -> t
  h %>%
    rep(length(t)) %>%
    data.frame(low=., high=t, frequency=as.integer(1)) %>%
    rbind(pair_(t))
}

pair <- function(v) {
  v %>%
    as.character %>%
    sort %>%
    pair_
}

count_pairs <- function(df) {
  df %>%
    colnames %>%
    map(function(colname){
      df %>%
      get(colname, .) %>%
      pair
    }) %>%
    reduce(rbind) %>%
    group_by(low,high) %>%
    summarise(frequency = sum(frequency)) %>%
    ungroup %>%
    arrange(-frequency)
}


genes_relevant %>%
  count_pairs %>%
  write_tsv(output_pairs_path)
