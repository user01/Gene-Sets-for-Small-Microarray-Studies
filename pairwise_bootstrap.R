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

args <- parser$parse_args()

random_seed <- args$seed
# random_seed <- 451
fraction_bootstrap <- args$bootstrap
# fraction_bootstrap <- 0.66
fraction_variation <- args$fraction
# fraction_variation <- 0.5
top_genes <- ceiling(1/2 * (sqrt(8*args$toppairs + 1) + 1)) # this uses nC2 to pick the top set to create ~number of pairs
# top_genes <- ceiling(1/2 * (sqrt(8*200 + 1) + 1))


genes <- file.path("results", "gene_data_vs_cell_type.tsv") %>%
  read_tsv(col_types=cols(
    .default = col_double(),
    GSM_ID = col_character(),
    Cell_Type = col_character(),
    General_Cell_Type = col_character()
  ))

genes %>% select(-GSM_ID, -Cell_Type, -General_Cell_Type) -> gene_data
genes %>% select( GSM_ID,  Cell_Type,  General_Cell_Type) -> gene_labels

gene_data %>% glimpse

set.seed(random_seed)
gene_data %>%
  ncol %>%
  `*`(fraction_bootstrap) %>%
  floor %>%
  sample.int(replace = TRUE) ->
  selected_genes_idx

gene_data %>%
  select(selected_genes_idx) ->
  gene_selected

gene_selected %>% glimpse

# gene_data_bootstrapped %>% glimpse
#
# expression_data <- gene_data_bootstrapped %>%
#                    select(-GSM_ID,-Cell_Type,-General_Cell_Type)
#
# expression_data %>% glimpse

pca_results <- prcomp(gene_selected, scale = TRUE)
pca_results %>% glimpse

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
  write_tsv(path_target)



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
  res

res %>%
  glimpse
  # as.data.frame %>%
  # select(1:pca_relevant_components_length) %>%
  # glimpse


pair_ <- function(v) {
  if (length(v) < 2) {
    return(list())
  }
  v %>% head(1) -> h
  v %>% tail(-1) %>% sort -> t
  h %>%
    rep(length(t)) %>%
    list(., t) %>%
    transpose %>%
    map(unlist) %>%
    c(pair_(t))
}

pair <- function(v) {
  v %>%
    as.character %>%
    pair_
}

# res %>%
#   get("PC1", .) %>%
#   pair %>%
#   as.data.frame


sweep(aload, 2, colSums(aload), "/")

pca_results <- expression_data %>%
               scale(center = TRUE, scale = FALSE) %>%
               t %>%
               scale(center = FALSE, scale = TRUE) %>%
               t %>%
               prcomp(center = FALSE)
# pca_results %>% glimpse
# plot(pca_results$x[,1:2])

# dimreduction_pca_dimensions_2.tsv
path_target <- file.path("results",
  paste0("dimreduction_pca_dimensions_",
         args$dimensions,
         ".tsv"
         )
       )

pca_results$x[, 1:args$dimensions] %>%
  as.data.frame %>%
  cbind(gene_data_vs_cell_type %>%
  select(Cell_Type,General_Cell_Type)) %>%
  write_tsv(path_target)
