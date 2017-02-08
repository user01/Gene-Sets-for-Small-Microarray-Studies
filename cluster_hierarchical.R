suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(argparse)
})


parser <- ArgumentParser()
parser$add_argument("-m", "--methodcluster", type="character", default="complete",
    help="Type of cluster to create ('ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid')")

parser$add_argument("-d", "--methoddistance", type="character", default="euclidean",
    help="Type of distance to compute 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'")

parser$add_argument("-c", "--clusters", type="integer", default=8,
    help="Number of clusters to create")

parser$add_argument("-k", "--neighbors", type="integer", default=1,
    help="Number of nieghbors to consider for classifying cluster")

parser$add_argument("-n", "--name", type="character", required=TRUE,
    help="Name of dimension reduced data set. Used to locate input TSV file and write output TSV. Input TSV must conform to table with Cell_Type, General_Cell_Type, and any number of float fields.")

args <- parser$parse_args()

name <- args$name
# name <- "pca_dimensions_2"
method_cluster <- args$methodcluster
# method_cluster <- "complete"
method_distance <- args$methoddistance
# method_distance <- "euclidean"
clusters <- args$clusters
# clusters <- 8
neighbors <- args$neighbors
# neighbors <- 3

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
    dist(method = method_distance) %>%
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
    distances <- dist(data_floats_train, method = method_distance) # distance matrix
    fit <- hclust(distances, method=method_cluster)
    fit_clusters <- cutree(fit, k=clusters) # cut tree into clusters

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

data_results_path <- paste0("cluster_hierarchical_", name, "_clusters_", clusters, "_methodcluster_", method_cluster, "_methoddistance_", method_distance, "_neighbors_", neighbors, ".tsv") %>%
  file.path("results", .)

data_results %>%
  write_tsv(data_results_path)
