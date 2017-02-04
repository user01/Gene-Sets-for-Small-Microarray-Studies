suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(argparse)
  library(Rtsne)
})


parser <- ArgumentParser()
parser$add_argument("-s", "--seed", type="integer", default=451,
    help="Random seed value")

parser$add_argument("-p", "--perplexity", type="integer", default=30,
    help="T-SNE perplexity number")

parser$add_argument("-i", "--max_iter", type="integer", default=1000,
    help="Maximum number of iterations")

parser$add_argument("-c", "--pca", type="character", default="false",
    help="Use PCA in first pass")

args <- parser$parse_args()

gene_data_vs_cell_type <- file.path("results", "gene_data_vs_cell_type.tsv") %>%
  read_tsv(col_types=cols(
    .default = col_double(),
    GSM_ID = col_character(),
    Cell_Type = col_character(),
    General_Cell_Type = col_character()
  ))

expression_data <- gene_data_vs_cell_type %>%
                   select(-GSM_ID,-Cell_Type,-General_Cell_Type) %>%
                   as.matrix


set.seed(args$seed) # Sets seed for reproducibility
tsne_out <- Rtsne(expression_data,
                  perplexity = args$perplexity,
                  max_iter = args$max_iter,
                  pca = args$pca == "true")

# dimreduction_tsne_perplexity_30_pca_false.tsv
path_target <- file.path("results",
  paste0("dimreduction_tsne_perplexity_",
         args$perplexity,
         "_pca_",
         args$pca,
         ".tsv"
         )
       )

tsne_out %>%
  get("Y", .) %>%
  as.data.frame %>%
  cbind(gene_data_vs_cell_type %>%
  select(Cell_Type,General_Cell_Type)) %>%
  write_tsv(path_target)
