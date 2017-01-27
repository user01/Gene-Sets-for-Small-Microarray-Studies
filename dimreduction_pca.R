suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

gene_data_vs_cell_type <- file.path("results", "gene_data_vs_cell_type.tsv") %>% read_tsv(col_types=cols(
  .default = col_double(),
  GSM_ID = col_character(),
  Cell_Type = col_character(),
  General_Cell_Type = col_character()
))

expression_data <- gene_data_vs_cell_type %>%
                   select(-GSM_ID,-Cell_Type,-General_Cell_Type) %>% as.matrix

pca_results <- expression_data %>%
               scale(center = TRUE, scale = FALSE) %>%
               t %>%
               scale(center = FALSE, scale = TRUE) %>%
               t %>%
               prcomp(center = FALSE)

# pca_results %>% glimpse
# plot(pca_results$x[,1:2])

pca_results$x[, 1:2] %>%
  as.data.frame %>%
  cbind(gene_data_vs_cell_type %>% select(Cell_Type,General_Cell_Type)) %>%
  write_tsv(file.path("results", "dimreduced_matrix_pca_1.2.tsv"))
