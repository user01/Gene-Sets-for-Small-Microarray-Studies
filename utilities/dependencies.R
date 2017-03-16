
list.of.packages <- c("readr", "dplyr", "stringr", "purrr", "argparse", "Rtsne",
  "MASS", "lubridate", "randomForest", "EMCluster")

# http://stackoverflow.com/a/4090208/2601448
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages, repos="http://cran.rstudio.com/")
