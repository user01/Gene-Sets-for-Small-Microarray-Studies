suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(argparse)
})


parser <- ArgumentParser()

parser$add_argument("-s", "--scorespath", type="character",
    default="pairwise",
    help="Path to input score files")
parser$add_argument("-p", "--scorespattern", type="character",
    default="score_general_\\d+.tsv",
    help="Pattern for score files")

parser$add_argument("-o", "--output", type="character",
    default=file.path("pairwise", "set_results.gmt"),
    help="Path to output set data")

args <- parser$parse_args()

scorespath <- args$scorespath
# scorespath <- "pairwise"
scorespattern <- args$scorespattern
# scorespattern <- "score_general_\\d+.tsv"
output <- args$output
# output <- file.path("pairwise", "set_results.gmt")


scorespath %>%
  list.files(recursive=TRUE, pattern=scorespattern) %>%
  map(function(path) {
    path %>%
      file.path(scorespath, .) %>%
      read_tsv(col_types = cols(
        .default = col_character(),
        low = col_character(),
        high = col_character(),
        frequency = col_integer(),
        score = col_double()
      ))
  }) %>%
  reduce(rbind) ->
  scores


scores %>%
  colnames %>%
  keep(~ str_detect(., "Cell_Type")) ->
  cell_type


scores %>%
  group_by_("low", "high", cell_type) %>%
  summarise(
    frequency_sum = sum(frequency),
    frequency_mean = mean(frequency),
    score_sum = sum(score),
    score_mean = mean(score),
    score = sum(frequency * score)
  ) %>%
  ungroup ->
  scores_grouped


scores_grouped %>%
  select_(cell_type) %>%
  unlist %>%
  unique ->
  cell_types

# TODO: Update threshold for score, currently 3rd Quartile
score_threshold <- function(scores) {
  summary(scores)[5] %>% unname
}
add_pair_to_sets <- function(pair, sets){
  map(sets, function(set) {
    is.element(pair, set) %>%
      sum %>%
      `>`(0) ->
      pair_matches
    if (pair_matches) {
      return(unique(c(set,pair)))
    }
    set
  })
}
contains_intersection <- function(sets, pattern) {
  sets %>% unlist %>% intersect(pattern) %>% length %>% `>`(0)
}
pairs_to_set <- function(pairs) {
  c(pairs$low,pairs$high) %>% unique
}

collapse_sets <- function(sets) {
  collapse_sets_(sets, list())
}
collapse_sets_ <- function(unfinished_sets, known_sets) {
  if (length(unfinished_sets) < 1) {
    return(known_sets)
  }
  top_set <- unfinished_sets %>% head(1) %>% unlist
  other_sets <- unfinished_sets %>% tail(-1)

  if (contains_intersection(other_sets, top_set)) {
    known_clean <- known_sets
    remaining_sets <- map(other_sets, function(set) {
      if (intersect(set, top_set) %>% length %>% `>`(0)) {
        return(unique(c(set,top_set)))
      }
      set
    })
  } else {
    known_clean <- c(known_sets, list(top_set))
    remaining_sets <- other_sets
  }
  collapse_sets_(remaining_sets, known_clean)
}


pairs_to_sets <- function(pairs) {
  pairs_to_sets_(pairs, list())
}
pairs_to_sets_ <- function(pairs, sets) {
  if (nrow(pairs) < 1) {
    return(sets)
  }

  # Take the first pair as a char vector
  top_pair <- pairs %>% head(1) %>% pairs_to_set

  # add pair to any set that contains any of the pair
  new_sets <- if (contains_intersection(sets, top_pair)) {
    add_pair_to_sets(top_pair, sets)
  } else {
    c(sets, list(top_pair))
  }

  # collapse sets that have any intersection
  collapsed_sets <- collapse_sets(new_sets)

  new_pairs <- pairs %>% tail(-1)
  pairs_to_sets_(new_pairs, collapsed_sets)
}



cell_types %>%
  map(function(type){
    scores_grouped %>%
      filter(scores_grouped[cell_type] == type) %>%
      select(low,high,score) ->
      pairs
    threshold <- pairs %>% select(score) %>% unlist %>% score_threshold
    pairs_relevant <- pairs %>% filter(score > threshold)
    sets <- pairs_relevant %>% pairs_to_sets

    list(
      type = type,
      sets = sets,
      count = pairs_relevant %>% nrow
    )
  }) ->
  cell_sets


cell_sets %>%
  discard(function(results) {
    results %>% get("sets", .) %>% length %>% `<`(1)
  }) %>%
  map(function(results) {
    results %>%
      get("sets", .) %>%
      map_chr(function(genes){
        str_c(genes, collapse = '\t')
      }) %>%
      str_c(results$type, "\t\t", .)
  }) %>%
  unlist %>%
  write_lines(output)
