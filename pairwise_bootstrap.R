suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("-s", "--seed", type="integer", default=451,
    help="Random seed")

parser$add_argument("-t", "--topgenes", type="integer", default=200,
    help="Number of genes to pick as relevant from a component")

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
top_genes <- args$topgenes
# top_genes <- 200


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
  select(selected_genes) ->
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
  pca_relevant_components


pca_results$sdev^2 %>% glimpse
prvar <- pca_results$sdev^2
# prvar
pve <- prvar/sum(prvar)
pve %>% glimpse
cumsum(pve)
plot(cumsum(pve), xlab="Principal Component ", ylab=" Cumulative Proportion of Variance Explained ", ylim=c(0,1), type='b')


pca_results$x[, 1:4] %>% as.data.frame %>% slice(1:5)


pca_results$rotation %>%
  glimpse

pca_results %>%
  get('rotation', .) %>%
  abs %>%
  (function(aload) { sweep(aload, 2, colSums(aload), "/") }) %>% # proportional contribution to the each principal component
  (function(m) { m[,1:pca_relevant_components] }) %>%
  (function(m) { m[,3] }) %>%
  sort(decreasing = TRUE) %>%
  head(5) %>%
  names
  # as.data.frame %>%
  # select(1:pca_relevant_components) %>%
  # glimpse


1:50 %>% head(4)

pca_results %>%
  get('rotation', .) %>%
  abs %>%
  (function(aload) { sweep(aload, 2, colSums(aload), "/") }) %>% # proportional contribution to the each principal component
  (function(m) { # Pick top top_genes genes from each
    m %>%
      colnames %>%
      head(pca_relevant_components) ->
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
  }) %>%
  glimpse
  # as.data.frame %>%
  # select(1:pca_relevant_components) %>%
  # glimpse


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
