suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("-c", "--clusters", type="integer", default=8,
    help="Number of clusters to create")

parser$add_argument("-i", "--input", type="character",
    default=file.path("results", "results_pca.tsv"),
    help="Name of dimension reduced data set. Used to locate input TSV file and write output TSV. Input TSV must conform to table with Cell_Type, General_Cell_Type, and any number of float fields.")

parser$add_argument("-o", "--output", type="character",
    default=file.path("results", "results_kmeans.tsv"),
    help="Path to output prediction results")

args <- parser$parse_args()

input_path <- args$input
# input_path <- file.path("results", "results_pca.tsv")
output_path <- args$output
# output_path <- file.path("results", "results_kmeans.tsv")
clusters <- args$clusters
# clusters <- 8


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

data_floats <- data_target %>% select(-Cell_Type, -General_Cell_Type)
data_labels <- data_target %>% select(Cell_Type, General_Cell_Type)
data_range <- 1:nrow(data_target)
data_floats %>%
  colnames %>%
  length -> dimensions_available

get_top <- function(v) {
  v %>%
    table %>%
    as.data.frame -> df

  if (nrow(df) < 1) {
    return(NA)
  }

  df %>%
    get(".",.) %>%
    nth(1)
}

classify <- function(training_floats,
                     training_labels,
                     testing_floats,
                     testing_labels,
                     clusters) {

  training_floats %>%
    aggregate(by=list(clusters), FUN=mean) ->
    agg_clusters

  agg_clusters %>%
    select(-`Group.1`) %>%
    rbind(testing_floats) %>%
    data.frame %>%
    dist(method = 'euclidean') %>%
    as.matrix %>%
    as.data.frame %>%
    slice(1:(nrow(.)-1)) %>%
    select(ncol(.)) %>%
    unname %>%
    unlist ->
    distances # distances that match with the cluster centers

  data.frame(distances=distances) %>%
    cbind(agg_clusters) %>%
    arrange(distances) %>%
    slice(1) %>%
    select(`Group.1`) %>%
    unlist %>%
    unname ->
    picked_cluster


  data.frame(cluster=clusters) %>%
    cbind(training_labels) %>%
    filter(cluster == picked_cluster) %>%
    select(-cluster) ->
    labels_in_cluster

  labels_in_cluster %>%
    select(Cell_Type) %>%
    unlist %>%
    get_top ->
    cell_type_predicted
  labels_in_cluster %>%
    select(General_Cell_Type) %>%
    unlist %>%
    get_top ->
    general_cell_type_predicted

  data.frame(
    Cell_Type_Predicted=cell_type_predicted,
    General_Cell_Type_Predicted=general_cell_type_predicted)
}

data_range %>%
  map(function(idx){

    # Create data sets
    data_floats_train <- data_floats %>% filter(idx != data_range)
    data_labels_train <- data_labels %>% filter(idx != data_range)
    data_floats_test <- data_floats %>% filter(idx == data_range)
    data_labels_test <- data_labels %>% filter(idx == data_range)

    # Perform cluster
    fit <- kmeans(data_floats_train, clusters)
    fit_clusters <- fit$cluster

    # Classify remaining point
    classified <- classify(data_floats_train,
                           data_labels_train,
                           data_floats_test,
                           data_labels_test,
                           fit_clusters)

    data.frame(
        idx = idx,
        Cell_Type = data_labels_test$Cell_Type,
        General_Cell_Type = data_labels_test$General_Cell_Type
      ) %>%
      cbind(classified)
  }) %>%
  reduce(rbind) %>%
  as.data.frame %>%
  select(idx, Cell_Type, Cell_Type_Predicted, General_Cell_Type, General_Cell_Type_Predicted) ->
  data_results

results_str <- function(truth, predicted) {
  cnt <- sum(as.character(truth) == as.character(predicted), na.rm = TRUE)
  total = length(truth)
  paste0(cnt,"/",total," (",round(100 * cnt/total), "%)")
}
gct_res <- results_str(data_results$General_Cell_Type, data_results$General_Cell_Type_Predicted)
ct_res <- results_str(data_results$Cell_Type, data_results$Cell_Type_Predicted)
paste0("For ", name," with ", clusters,
       " clusters, General Cell Type success was ", gct_res,
       " and Cell Type was ", ct_res, ".") %>%
       print()


data_results %>%
  write_tsv(output_path)
