suppressPackageStartupMessages({
  library(EMCluster)
  library(readr)
  library(dplyr)
  library(stringr)
})


args <- commandArgs(trailingOnly = TRUE)

target_file <- if (str_detect(args[1], "\\.tsv$")) {
  args[1]
} else {
  file.path("results", "dimreduced_matrix_pca_1.2.tsv")
}


target_data <- read_tsv(target_file, col_types = cols(PC1 = col_double(), PC2 = col_double()))
target_data %>% glimpse


demo(allinit, 'EMCluster', ask = F, echo = F)

assign.class(target_data)

set.seed(1234)

da2$da
x2 <- da2$da
emobj <- emgroup(x2, nclass = 5)
ret_2 <- emcluster(x2, emobj, assign.class = TRUE)

emobj %>% glimpse

plotem(ret_2, x2)

data.frame(x=5,y=14) %>%
  assign.class(emobj) %>%
  glimpse

ret_2 %>% glimpse
ret_2$class %>% glimpse
ret_2$class %>% unique %>% glimpse

cbind(x2,ret_2$class) %>% glimpse

ret <- init.EM(x2, nclass = 2)
ret.new <- assign.class(x2, ret, return.all = FALSE)
str(ret.new)

emcluster(da2)
