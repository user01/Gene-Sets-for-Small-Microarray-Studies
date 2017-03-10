suppressPackageStartupMessages({
  library(MASS)
  library(readr)
  library(stringr)
  library(dplyr)
  library(purrr)
  library(lubridate)
  library(argparse)
})

timestamp_started <- now()

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

parser$add_argument("-f", "--feedback", type="character",
    default=file.path("results", "results_feedback.tsv"),
    help="Path to output feedback on run. Errors, time, etc")

args <- parser$parse_args()

random_seed <- args$seed
# random_seed <- 451
fraction_bootstrap <- args$bootstrap
# fraction_bootstrap <- 0.66
top_genes <- ceiling(1/2 * (sqrt(8*args$pairs + 1) + 1)) # this uses nC2 to pick the top set to create ~number of pairs
# top_genes <- ceiling(1/2 * (sqrt(8*200 + 1) + 1))
input_path <- args$input
# input_path <- file.path("results", "gene_data_vs_cell_type.tsv")
cell_name <- args$name %>% str_replace_all('"','')
# cell_name <- '"Monocyte"' %>% str_replace_all('"','')
cell_type <- args$type
# cell_type <- "General_Cell_Type"
input_path <- args$input
# input_path <- file.path("results", "gene_data_vs_cell_type.tsv")
output_path <- args$output
# output_path <- file.path("results", "results.tsv")
feedback_path <- args$feedback
# feedback_path <- file.path("results", "feedback.tsv")



# function to note statistics on done
# this will stop the script
note_complete <- function(success, message) {
  data.frame(
      start = timestamp_started,
      end = now(),
      cell_name = cell_name,
      cell_type = cell_type,
      random_seed = random_seed,
      fraction_bootstrap = fraction_bootstrap,
      top_genes = top_genes,
      input_path = input_path,
      output_path = output_path,
      success = success,
      message = message
    ) %>%
    write_tsv(feedback_path)
  quit()
}

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
  sample_n(nrow(genes_bootstrapped_truth)) %>%
  rbind(genes_bootstrapped_truth) ->
  genes_bootstrapped_downsampled

# LDA fit for components
lda_fit <- function(df) {
  tryCatch({
    results <- lda(df %>%
                     select(-type_truth) %>%
                     as.matrix,
                   df %>%
                     select(type_truth) %>%
                     unlist %>%
                     unname %>%
                     as.character,
                   CV = FALSE)
    results
  }, error = function(errorCondition) {
    # Note this will force the script to stop
    note_complete(FALSE, errorCondition$message)
  })
}

# Select top loading genes from components
genes_bootstrapped_downsampled %>%
  lda_fit %>%
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

# Perform Leave One Out
data_range <- 1:nrow(genes_bootstrapped_downsampled)
data_range %>%
  map(function(idx){
    # Create data sets
    data_train <- genes_bootstrapped_downsampled %>% filter(idx != data_range)
    data_test <- genes_bootstrapped_downsampled %>% filter(idx == data_range)

    lda_model <- lda_fit(data_train)

    data_test %>%
      select(type_truth) %>%
      unlist %>%
      unname %>%
      as.character ->
      type_truth

    data_test %>%
      select(-type_truth) %>%
      predict(lda_model, newdata = .) %>%
      get("class", .) ->
      type_predicted

    data.frame(
      idx = idx,
      type_truth = type_truth,
      type_predicted = type_predicted
    )
  }) %>%
  reduce(rbind) %>%
  as.data.frame ->
  data_results

# Score is simply accuracy of the LOO predictions between the balanced classes
#  of the type and other
data_results %>%
  mutate(correct = type_truth == type_predicted) %>%
  select(correct) %>%
  unlist %>%
  sum %>%
  `/`(nrow(data_results)) ->
  score

# Pairs a vector into low/high pairs, with a frequency of one
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

top_gene_names %>%
  pair %>%
  mutate(
    score = score,
    cell_name = cell_name,
    cell_type = cell_type
  ) %>%
  write_tsv(output_path)

note_complete(TRUE, NA)
