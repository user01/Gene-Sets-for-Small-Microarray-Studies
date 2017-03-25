suppressPackageStartupMessages({
  library(GSEABase)
  library(GSVA)
  library(readr)
  library(dplyr)
  library(purrr)
  library(argparse)
  library(lazyeval)
})


parser <- ArgumentParser()

parser$add_argument("-i", "--input", type="character",
    default=file.path("results", "base", "base.set.feedback.General_Cell_Type.Neutrophil.tsv"),
    help="Path to input feedback set data")

parser$add_argument("-d", "--basedata", type="character",
    default=file.path("results", "base", "gene_data_vs_cell_type.tsv"),
    help="Path to base data set")

parser$add_argument("-v", "--validationdata", type="character",
    default=file.path("results", "base", "gene_data_vs_cell_type.validation.tsv"),
    help="Path to validation data set")

parser$add_argument("-s", "--source", type="character",
    default=file.path("results", "base"),
    help="Path to source of sets")

parser$add_argument("-o", "--output", type="character",
    default=file.path("results","base","base.set.results.General_Cell_Type.Neutrophil.tsv"),
    help="Path to output results")

args <- parser$parse_args()

input_path <- args$input
# input_path <- file.path("results", "base", "base.set.feedback.General_Cell_Type.Neutrophil.tsv")
basedata_path <- args$basedata
# basedata_path <- file.path("results", "base", "gene_data_vs_cell_type.tsv")
validationdata_path <- args$validationdata
# validationdata_path <- file.path("results", "base", "gene_data_vs_cell_type.validation.tsv")
source_path <- args$source
# source_path <- file.path("results", "base")
output_path <- args$output
# output_path <- file.path("results","base","base.set.results.General_Cell_Type.Neutrophil.tsv")

print("Started")

input_path %>%
  read_tsv(col_types = cols(
    index = col_integer(),
    genes_considered = col_integer(),
    genes_in_set = col_integer(),
    score_avg = col_double(),
    score_max = col_double(),
    score_min = col_double(),
    score_std = col_double(),
    cell_name = col_character(),
    cell_type = col_character(),
    title = col_character()
  )) ->
  set_data

cell_type <- set_data %>% get("cell_type", .) %>% first
cell_name <- set_data %>% get("cell_name", .) %>% first
title <- set_data %>% get("title", .) %>% first

print("Set Data")


data_read <- function(data_path){
  data_path %>%
    read_tsv(col_types = cols(
      .default = col_double(),
      GSM_ID = col_character(),
      Cell_Type = col_character(),
      General_Cell_Type = col_character()
    )) ->
    data

  data %>%
    select_(cell_type) %>%
    unlist %>%
    unname %>%
    map_chr(~ if (. == cell_name) { "target" } else { "other" }) %>%
    list(., 1:length(.)) %>%
    transpose %>%
    map_chr(~ paste0(.[1],.[2])) ->
    col_names

  data %>%
    select(-GSM_ID, -Cell_Type, -General_Cell_Type) %>%
    as.matrix %>%
    t %>%
    `colnames<-`(col_names)
}

data_sets_base <- data_read(basedata_path)
data_sets_validation <- data_read(validationdata_path)

print("main Data")


set_data %>%
  get('index', .) %>%
  map(function(idx){
    idx %>%
      sprintf("%05d", .) %>%
      paste0(title, ".set.", cell_type, ".", cell_name, ".", ., ".tsv") ->
      set_name

    set_name %>%
      file.path(source_path, .) %>%
      read_tsv(col_types = cols(
        gene = col_character()
      )) %>%
      get("gene", .) %>%
      unlist %>%
      GeneSet(.,
              setName = set_name,
              setIdentifier = set_name,
              geneIdType = ENSEMBLIdentifier())
  }) %>%
  GeneSetCollection ->
  gene_sets

summary(gene_sets)

print("Gene Sets")


enrichment_values <- function(gene_data) {
  gsva(gene_data,
       gene_sets,
       verbose = FALSE,
       parallel.sz = 1) %>%
    get('es.obs', .) %>%
    as.data.frame
}


base_enrichments <- enrichment_values(data_sets_base)
validation_enrichments <- enrichment_values(data_sets_validation)

print("enrichment_values")

# read_csv("something.csv") -> df
# read_tsv("results/base/base.set.results.General_Cell_Type.Microglia.tsv") -> df
#
# df %>% glimpse
# df %>%
#   select(starts_with("other")) %>%
#   summarise_each(funs(mean)) %>%
#   glimpse
# df %>%
#   select(starts_with("other")) %>%
#   rowMeans %>%
#   abs %>%
#   glimpse
#
# df %>%
#   select(starts_with("target")) %>%
#   rowMeans %>%
#   abs %>%
#   glimpse
#
#
# df %>%
#   select(starts_with("other")) %>%
#   mutate(x = rowMeans(.))

enrichment_mean <- function(df, prefix) {
  df %>%
    select(starts_with(prefix)) %>%
    rowMeans %>%
    abs
}
score_sets <- function(df) {
  # compare target vs other
  scores_target <- enrichment_mean(df, "target")
  scores_other <- enrichment_mean(df, "other")
  # score of the target vs other samples
  scores_target - scores_other
}

base_scores <- score_sets(base_enrichments)
validation_scores <- score_sets(validation_enrichments)



ci <- function(vect) {
  #confidence interval
  n = length(vect)
  z = 1.96 #this represents alpha = 0.95
  std_dev <- sd(vect)

  c(
    #upper bound
    mean(vect) + z*std_dev/(n^0.5),
    #lower bound
    mean(vect) - z*std_dev/(n^0.5)
  )
}
ci_upper <- function(vect) { ci(vect)[1] }
ci_lower <- function(vect) { ci(vect)[2] }
#
# df %>%
#   t %>%
#   as.data.frame %>%
#   summarise_each(funs(ci_lower)) %>%
#   unlist %>%
#   unname ->
#   lower_ci
#
# df %>%
#   t %>%
#   as.data.frame %>%
#   summarise_each(funs(ci_lower)) %>%
#   unlist %>%
#   unname ->
#   upper_ci




set_data %>%
  mutate(
    base_scores = base_scores,
    validation_scores = validation_scores
  ) %>%
  write_tsv(output_path)


#
# gsva(data_sets_base$data_other, gene_sets, mx.diff=1, parallel.sz = 1) %>%
#   get('es.obs', .) ->
#   enrichment_gautier_other

#
# # results/base/base.set.feedback.General_Cell_Type.Neutrophil.tsv
#
# p <- 10 ## number of genes
# n <- 30 ## number of samples
# nGrp1 <- 15 ## number of samples in group 1
# nGrp2 <- n - nGrp1 ## number of samples in group 2
#
# matrix(rnorm(3*4), nrow=4, ncol=3,dimnames=list(paste("g", 1:4, sep="") , paste("s", 1:3, sep="")))
#
# y <- matrix(rnorm(n*p), nrow=p, ncol=n,dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
#
# geneSets <- list(set1=paste("g", 1:3, sep=""),
#                  set2=paste("g", 4:6, sep=""),
#                  set3=paste("g", 7:10, sep=""))
#
# res <- gsva(y, geneSets, mx.diff=1)
#
#
# gene_data <- matrix(rnorm(n*p), nrow=p, ncol=n,dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))
#
# gene_sets <- list(set1=paste("g", 1:3, sep=""),
#                   set2=paste("g", 4:6, sep=""),
#                   set3=paste("g", 7:10, sep=""))
#
# res <- gsva(gene_data, geneSets, mx.diff=1)
