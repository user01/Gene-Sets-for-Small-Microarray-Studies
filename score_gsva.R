suppressPackageStartupMessages({
  library(GSEABase)
  library(GSVA)
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(argparse)
  library(lazyeval)
})


parser <- ArgumentParser()

parser$add_argument("-i", "--input", type="character",
    default=file.path("results", "base", "base.set.feedback.General_Cell_Type.Neutrophil.tsv"),
    help="Path to input feedback set data")

parser$add_argument("-d", "--basedata", type="character",
    default=file.path("results", "base", "gene_data_vs_cell_type.tsv"),
    help="Path to base data set")

parser$add_argument("-v", "--validationdata", type="character",
    default=file.path("results", "base", "gene_data_vs_cell_type.validation.tsv"),
    help="Path to validation data set")

parser$add_argument("-s", "--source", type="character",
    default=file.path("results", "base"),
    help="Path to source of sets")

parser$add_argument("-o", "--output", type="character",
    default=file.path("results","base","base.set.results.General_Cell_Type.Neutrophil.tsv"),
    help="Path to output results")

args <- parser$parse_args()

input_path <- args$input
# input_path <- "results/full/full.set.feedback.Cell_Type.SI_Serosal_Mf.tsv"
basedata_path <- args$basedata
# basedata_path <- "results/full/gene_data_vs_cell_type.tsv"
validationdata_path <- args$validationdata
# validationdata_path <- "results/full/gene_data_vs_cell_type.tsv"
source_path <- args$source
# source_path <- "results/full"
output_path <- args$output
# output_path <- "full.set.results.Cell_Type.SI_Serosal_Mf.tsv"

print("Started")

input_path %>%
  read_tsv(col_types = cols(
    index = col_integer(),
    genes_considered = col_integer(),
    genes_in_set = col_integer(),
    score_avg = col_double(),
    score_max = col_double(),
    score_min = col_double(),
    score_std = col_double(),
    cell_name = col_character(),
    cell_type = col_character(),
    title = col_character()
  )) ->
  set_data

cell_type <- set_data %>% get("cell_type", .) %>% first
cell_name <- set_data %>% get("cell_name", .) %>% first
title <- set_data %>% get("title", .) %>% first

print("Set Data")

data_read <- function(data_path){
  data_path %>%
    read_tsv(col_types = cols(
      .default = col_double(),
      GSM_ID = col_character(),
      Cell_Type = col_character(),
      General_Cell_Type = col_character()
    ))
}

data_base_raw <- data_read(basedata_path)
data_validation_raw <- data_read(validationdata_path)


data_process <- function(data, target_name, target_type){
  data %>%
    # select_(cell_type) %>%
    select_(target_type) %>%
    unlist %>%
    unname %>%
    map_chr(~ if (. == target_name) { "target" } else { "other" }) %>%
    list(., 1:length(.)) %>%
    transpose %>%
    map_chr(~ paste0(.[1],.[2])) ->
    col_names

  data %>%
    select(-GSM_ID, -Cell_Type, -General_Cell_Type) %>%
    as.matrix %>%
    t %>%
    `colnames<-`(col_names)
}

cell_name_fixed <- str_replace_all(cell_name, "_", " ")
data_sets_base <- data_process(data_base_raw, cell_name_fixed, cell_type)
data_sets_validation <- data_process(data_validation_raw, cell_name_fixed, cell_type)


c("Macrophage","Microglia","Neutrophil","Monocyte") %>%
  map(function(current_type){
    list(
      name=current_type,
      mat_train=data_process(data_base_raw, current_type, "General_Cell_Type"),
      mat_validation=data_process(data_validation_raw, current_type, "General_Cell_Type")
    )
  }) ->
  four_vs_others


set_data %>%
  get('index', .) %>%
  map(function(idx){
    idx %>%
      sprintf("%05d", .) %>%
      paste0(title, ".set.", cell_type, ".", cell_name, ".", ., ".tsv") ->
      set_name

    set_name %>%
      file.path(source_path, .) %>%
      read_tsv(col_types = cols(
        gene = col_character()
      )) %>%
      get("gene", .) %>%
      unlist %>%
      GeneSet(.,
              setName = set_name,
              setIdentifier = set_name,
              geneIdType = ENSEMBLIdentifier())
  }) %>%
  GeneSetCollection ->
  gene_sets


enrichment_values <- function(gene_data) {
  gsva(gene_data,
       gene_sets,
       verbose = FALSE,
       parallel.sz = 1) %>%
    get('es.obs', .) %>%
    as.data.frame
}


base_enrichments <- enrichment_values(data_sets_base)
validation_enrichments <- enrichment_values(data_sets_validation)


enrichment_mean <- function(df, prefix) {
  df %>%
    select(starts_with(prefix)) %>%
    rowMeans %>%
    abs
}
score_sets <- function(df) {
  # compare target vs other
  scores_target <- enrichment_mean(df, "target")
  scores_other <- enrichment_mean(df, "other")
  # score of the target vs other samples
  scores_target - scores_other
}

base_scores <- score_sets(base_enrichments)


data_sets_validation %>%
  as.data.frame %>%
  select(starts_with("target")) %>%
  ncol ->
  count_of_targets

# stop(paste("Count of targets: ", count_of_targets))

validation_scores <- if(count_of_targets > 0) {
  score_sets(validation_enrichments)
} else {
  warning("Insufficent data for validation")
  1:length(gene_sets) * 0
}


four_vs_others %>%
  map(function(payload){
    name <- payload %>% get('name', .)
    mat_train <- payload %>% get('mat_train', .)
    mat_validation <- payload %>% get('mat_validation', .)

    enrichment_train <- enrichment_values(mat_train)
    enrichment_validation <- enrichment_values(mat_validation)
    score_train <- score_sets(enrichment_train)
    score_validation <- score_sets(enrichment_validation)
    data.frame(
      score_train = score_train,
      score_validation = score_validation
    ) %>%
    rename_(.dots = setNames("score_train", paste0("score_train_", name))) %>%
    rename_(.dots = setNames("score_validation", paste0("score_validation_", name)))
  }) %>%
  reduce(cbind) ->
  four_vs_others_df


set_data %>%
  select(-score_avg, -score_max, -score_min, -score_std) %>%
  mutate(
    base_score = base_scores,
    validation_score = validation_scores
  ) %>%
  cbind(
    .,
    four_vs_others_df
  ) %>%
  write_tsv(output_path)
