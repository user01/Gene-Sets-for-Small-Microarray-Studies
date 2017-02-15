suppressPackageStartupMessages({
  library(randomForest)
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(argparse)
})


parser <- ArgumentParser()
parser$add_argument("-s", "--seed", type="integer", default=451,
    help="Random seed offset")

parser$add_argument("-t", "--trees", type="integer", default=64,
    help="Number of trees to grow")

parser$add_argument("-i", "--input", type="character", default=file.path("results", "results_components.tsv"),
    help="Path to input data")

parser$add_argument("-o", "--output", type="character", default=file.path("results", "classification_randomforest.tsv"),
    help="Path to output data")


args <- parser$parse_args()

random_seed <- args$seed
# random_seed <- 451
name <- args$name
# name <- "pca_dimensions_2"
ntree <- args$trees
# ntree <- 64
input_path <- args$input
# input_path <- file.path("results", "pca_results.tsv")
output_path <- args$output
# output_path <- file.path("results", "cluster_results.tsv")



data_target <- input_path %>%
  read_tsv(col_types = cols(
    .default = col_double(),
    Cell_Type = col_character(),
    General_Cell_Type = col_character()
  )) %>%
  mutate(
    Cell_Type = as.factor(Cell_Type),
    General_Cell_Type = as.factor(General_Cell_Type)
  )


data_target %>%
  select(-Cell_Type, -General_Cell_Type) %>%
  colnames %>%
  paste(collapse = " + ") -> formula_columns

formula_generate <- function(cols,y_name) {
  cols %>%
    c(y_name, .) %>%
    paste(collapse = " ~ ") %>%
    as.formula
}

formula_Cell_Type <- formula_generate(formula_columns, "Cell_Type")
formula_General_Cell_Type <- formula_generate(formula_columns, "General_Cell_Type")


data_range <- 1:nrow(data_target)

data_range %>%
  map(function(idx){
    set.seed(random_seed + idx)

    # Create data sets
    data_train <- data_target %>% filter(idx != data_range)
    data_test <- data_target %>% filter(idx == data_range)

    # Perform Random Forests
    rf_Cell_Type <- randomForest(formula_Cell_Type, data = data_train,
      importance = TRUE,
      ntree = ntree)
    rf_General_Cell_Type <- randomForest(formula_General_Cell_Type, data = data_train,
      importance = TRUE,
      ntree = ntree)

    rf_Cell_Type %>%
      predict(data_test) %>%
      unname -> predicted_class_cell_type
    rf_General_Cell_Type %>%
      predict(data_test) %>%
      unname -> predicted_class_general_cell_type

    data.frame(
      idx = idx,
      Cell_Type = data_test$Cell_Type,
      Cell_Type_Predicted = predicted_class_cell_type,
      General_Cell_Type = data_test$General_Cell_Type,
      General_Cell_Type_Predicted = predicted_class_general_cell_type
    )
  }) %>%
  reduce(rbind) %>%
  as.data.frame ->
  data_results

data_results %>%
  write_tsv(output_path)
