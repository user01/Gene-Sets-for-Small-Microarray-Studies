suppressPackageStartupMessages({
  library(EMCluster)
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(argparse)
})


parser <- ArgumentParser()
parser$add_argument("-c", "--clusters", type="integer", default=8,
    help="Number of clusters to create")

parser$add_argument("-n", "--name", type="character", required=TRUE,
    help="Name of dimension reduced data set. Used to locate input TSV file and write output TSV. Input TSV must conform to table with Cell_Type, General_Cell_Type, and any number of float fields.")
args <- parser$parse_args()


data_target <- args$name %>%
  paste0("dimreduced_", ., ".tsv") %>%
  file.path("results", .) %>%
  read_tsv(col_types = cols(
    .default = col_double(),
    Cell_Type = col_character(),
    General_Cell_Type = col_character()
  ))
# data_target %>% glimpse


data_pca <- data_target %>% select(-Cell_Type, -General_Cell_Type)
data_labels <- data_target %>% select(Cell_Type,General_Cell_Type)
# data_pca %>% glimpse
# data_labels %>% glimpse

get_top <- function(v) {
  v %>%
    table %>%
    as.data.frame %>%
    get(".",.) %>%
    nth(1)
}

data_range <- 1:nrow(data_target)

data_range %>%
  map(function(idx){

    # Create data sets
    data_pca_current <- data_pca %>% filter(idx != data_range)
    data_labels_current <- data_labels %>% filter(idx != data_range)
    data_pca_target <- data_pca %>% filter(idx == data_range)
    data_labels_target <- data_labels %>% filter(idx == data_range)

    # Perform cluster
    emobj <- emgroup(data_pca_current, nclass = args$clusters)
    ret <- emcluster(data_pca_current, emobj, assign.class = TRUE)

    data_pca_target %>%
      assign.class(emobj) %>%
      get("class",.) -> pca_predicted_class

    data_labels_current %>%
      mutate(class = emobj %>% get("class",.)) %>%
      filter(class == pca_predicted_class) -> predicted_class_training

    predicted_class_cell_type <- predicted_class_training %>%
      get("Cell_Type",.) %>%
      get_top
    predicted_class_general_cell_type <- predicted_class_training %>%
      get("General_Cell_Type",.) %>%
      get_top

    # # Render cluster plot
    # png(filename=paste0("results/plot",idx,".png"))
    # plotem(ret, data_pca_current, main=paste0("Cluster ", idx, " Results"))
    # points(data_pca_target, pch=23, col="red", bg="red")
    # dev.off()

    data.frame(
      idx = idx,
      Cell_Type = data_labels_target$Cell_Type,
      Cell_Type_Predicted = predicted_class_cell_type,
      General_Cell_Type = data_labels_target$General_Cell_Type,
      General_Cell_Type_Predicted = predicted_class_general_cell_type
    )
  }) %>%
  reduce(rbind) %>%
  as.data.frame ->
  data_results


data_results_path <- args$name %>%
  paste0("cluster_em_",.,".tsv") %>%
  file.path("results", .)

data_results %>%
  write_tsv(data_results_path)
