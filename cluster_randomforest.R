suppressPackageStartupMessages({
  library(randomForest)
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(argparse)
})


parser <- ArgumentParser()
parser$add_argument("-t", "--trees", type="integer", default=64,
    help="Number of trees to grow")

parser$add_argument("-n", "--name", type="character", required=TRUE,
    help="Name of dimension reduced data set. Used to locate input TSV file and write output TSV. Input TSV must conform to table with Cell_Type, General_Cell_Type, and any number of float fields.")

args <- parser$parse_args()

name <- args$name
# name <- "pca_dimensions_2"
ntree <- args$trees
# ntree <- 64

data_target <- name %>%
  paste0("dimreduction_", ., ".tsv") %>%
  file.path("results", .) %>%
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


results_str <- function(truth, predicted) {
  cnt <- sum(as.character(truth) == as.character(predicted), na.rm = TRUE)
  total = length(truth)
  paste0(cnt,"/",total," (",round(100 * cnt/total), "%)")
}
gct_res <- results_str(data_results$General_Cell_Type, data_results$General_Cell_Type_Predicted)
ct_res <- results_str(data_results$Cell_Type, data_results$Cell_Type_Predicted)
paste0("For ", name," with ", ntree,
       " trees, General Cell Type success was ", gct_res,
       " and Cell Type was ", ct_res, ".") %>%
       print()

data_results_path <- paste0("cluster_randomforest_", name, "_trees_", ntree, ".tsv") %>%
  file.path("results", .)

data_results %>%
  write_tsv(data_results_path)
