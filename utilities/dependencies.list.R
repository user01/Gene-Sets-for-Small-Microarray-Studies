suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

c("MASS", "readr", "stringr", "dplyr", "purrr", "lubridate", "argparse", "readr",
  "dplyr", "stringr", "purrr", "argparse", "Rtsne", "MASS", "lubridate", "randomForest",
  "EMCluster", "sparsediscrim") %>% unique %>% map(function(lib) {
  lib %>% packageVersion %>% paste0(" * `", lib, " :: v", .,"`")
}) %>% paste(collapse = "-")
