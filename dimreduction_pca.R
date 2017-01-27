suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

expression_data <- file.path("results", "expression_matrix.tsv") %>% read_tsv(col_names = FALSE) %>%
  as.matrix

pca_results <- expression_data %>% t %>% scale(center = TRUE, scale = FALSE) %>%
  t %>% scale(center = FALSE, scale = TRUE) %>% t %>% prcomp(center = FALSE)

# pca_results %>% glimpse
# plot(pca_results$x[,1:2])

pca_results$x[, 1:2] %>%
  as.data.frame %>%
  write_tsv(file.path("results", "dimreduced_matrix_pca_1.2.tsv"))
