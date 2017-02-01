suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
})


file_paths <- list.files("results", recursive=TRUE, pattern="cluster_(\\w+)_(\\w+)_(.+)\\.tsv")

file_paths %>%
  str_match_all("cluster_([^_]+)_([^_]+)_(.+)\\.tsv") %>%
  reduce(rbind) %>%
  `colnames<-`(c("filename","cluster_method","dimension_reduction_method","notes")) %>%
  as.data.frame ->
  file_data

file_data %>%
  get("filename", .) %>%
  map(function(filename) {
    file.path("results", filename) %>%
      read_tsv(col_types = cols(
        idx = col_integer(),
        Cell_Type = col_character(),
        Cell_Type_Predicted = col_character(),
        General_Cell_Type = col_character(),
        General_Cell_Type_Predicted = col_character()
      )) %>%
      mutate(
        Cell_Type_Correct = Cell_Type == Cell_Type_Predicted,
        General_Cell_Type_Correct = General_Cell_Type == General_Cell_Type_Predicted
      ) ->
      results
    data.frame(
      filename = filename,
      Cell_Type_Success = sum(results$Cell_Type_Correct) / nrow(results),
      General_Cell_Type_Success = sum(results$General_Cell_Type_Correct) / nrow(results)
    )
  }) %>%
  reduce(rbind) %>%
  inner_join(file_data) %>%
  select(cluster_method, dimension_reduction_method, notes, Cell_Type_Success, General_Cell_Type_Success) ->
  results


results_path <- file.path("results", "results_all.csv")

results %>%
  write_csv(results_path)
