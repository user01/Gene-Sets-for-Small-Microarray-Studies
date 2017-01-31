#Getting the Gene_normalized_expression formatted so column names are genes, row names are samples
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

EMF_GSVA_Relative_Normalilzed_Expression_PATH <- file.path("data", "EMF_GSVA_Relative_Normalilzed_Expression.tsv.tar.gz")
gene_normalized <- read_tsv(EMF_GSVA_Relative_Normalilzed_Expression_PATH,
                                                     col_types = cols(.default = col_double(), EMF_GSVA_Relative_Normalilzed_Expression.tsv = col_character(),
                                                                      Symbol = col_character(), Entrez = col_integer(), Name = col_character())) %>%
  filter(!is.na(Symbol)) %>% select(c(-1,-3,-4))

gene_normalized <- t(gene_normalized)

#making first row as header names:
colnames(gene_normalized) = gene_normalized[1, ] # the first row will be the header
gene_normalized = gene_normalized[-1, ]          # removing the first row.
class(gene_normalized) <- "numeric"              # make values numeric

#PCA
pr.out <- prcomp(gene_normalized, scale = TRUE)

#plotting cumulative explanation of variance of PCs
# pr.out$rotation = - pr.out$rotation
# pr.out$x = -pr.out$x
pr.out$sdev
pr.var = pr.out$sdev^2
pr.var
pve = pr.var/sum(pr.var)
pve
plot(cumsum(pve), xlab="Principal Component ", ylab=" Cumulative Proportion of Variance Explained ", ylim=c(0,1), type='b')

#Finding the top 40 Genes (based on the top loadings of the first 20 PCAs)
top40_genes = c()
for(i in 1:20){
  top40_genes <- c(top40_genes, names(tail(sort(pr.out$rotation[,i]))[1]), 
    names(head(sort(pr.out$rotation[,i]))[1]))
}
top40_genes

head(sort(pr.out$rotation[,1]))
