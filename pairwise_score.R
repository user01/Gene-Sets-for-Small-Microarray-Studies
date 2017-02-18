suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(argparse)
})

parser <- ArgumentParser()

parser$add_argument("-i", "--inputpairs", type="character",
    default=file.path("pairwise", "pairs_0.tsv"),
    help="Path to input pair data")

parser$add_argument("-p", "--inputpredictions", type="character",
    default=file.path("pairwise", "predict_0_randomforest.tsv"),
    help="Path to input data, string delimited paths.")

parser$add_argument("-g", "--outputgeneral", type="character",
    default=file.path("results", "score_general.tsv"),
    help="Path to output general score data")

parser$add_argument("-s", "--outputspecific", type="character",
    default=file.path("results", "score_specific.tsv"),
    help="Path to output specific score data")

args <- parser$parse_args()

input_pairs_path <- args$inputpairs
# input_pairs_path <- file.path("pairwise", "pairs_0.tsv")
input_predictions_path <- args$inputpredictions %>% str_split(",") %>% unlist
# input_predictions_path <- file.path("pairwise", "predict_0_randomforest.tsv") %>% str_split(",") %>% unlist
output_generalscore_path <- args$outputgeneral
# output_generalscore_path <- file.path("pairwise", "score_general.tsv")
output_specificscore_path <- args$outputspecific
# output_specificscore_path <- file.path("pairwise", "score_specific.tsv")


"Gautier_Immgen_Sample_Metadata.tsv" %>%
  file.path("data", .) %>%
  read_tsv(col_types = cols(
      GSM_ID = col_character(),
      Cell_Type = col_character(),
      General_Cell_Type = col_character()
    )
  ) ->
  Gautier_Immgen_Sample_Metadata


input_predictions_path %>%
  map(function(path) {
    path %>%
      read_tsv(col_types = cols(
        idx = col_integer(),
        Cell_Type = col_character(),
        Cell_Type_Predicted = col_character(),
        General_Cell_Type = col_character(),
        General_Cell_Type_Predicted = col_character()
      ))
  }) ->
  predictions
# predictions %>% glimpse

input_pairs_path %>%
  read_tsv(col_types = cols(
    low = col_character(),
    high = col_character(),
    frequency = col_integer()
  )) ->
  pairs
# pairs %>% glimpse



# Gautier_Immgen_Sample_Metadata %>% glimpse

scores_general <- function(general_cell_type) {
  predictions %>%
    map_dbl(function(df) {
      df %>%
        filter(General_Cell_Type == general_cell_type) %>%
        mutate(correct = General_Cell_Type == General_Cell_Type_Predicted) ->
        df_perf
      sum(df_perf$correct) / nrow(df_perf)
    })
}
scores_specific <- function(specific_cell_type) {
  predictions %>%
    map_dbl(function(df) {
      df %>%
        filter(Cell_Type == specific_cell_type) %>%
        mutate(correct = Cell_Type == Cell_Type_Predicted) ->
        df_perf
      sum(df_perf$correct) / nrow(df_perf)
    })
}


Gautier_Immgen_Sample_Metadata %>%
  get("General_Cell_Type", .) %>%
  unique %>%
  map(function(general_cell_type) {
    # compute score for this general cell type based on predictions
    # TODO: Currently makes use of the mean of prediction scores
    general_cell_type %>%
      scores_general %>%
      mean -> type_score
    pairs %>%
      mutate(
        General_Cell_Type = general_cell_type,
        score = type_score)
  }) %>%
  reduce(rbind) ->
  output_generalscore

output_generalscore %>%
  write_tsv(output_generalscore_path)

Gautier_Immgen_Sample_Metadata %>%
  get("Cell_Type", .) %>%
  unique %>%
  map(function(specific_cell_type) {
    # compute score for this general cell type based on predictions
    # TODO: Currently makes use of the mean of prediction scores
    specific_cell_type %>%
      scores_general %>%
      mean -> type_score
    pairs %>%
      mutate(
        Cell_Type = specific_cell_type,
        score = type_score)
  }) %>%
  reduce(rbind) ->
  output_specificscore

output_specificscore %>%
  write_tsv(output_specificscore_path)
