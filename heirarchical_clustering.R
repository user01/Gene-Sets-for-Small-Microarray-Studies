#assume gene_data_vs_cell_type dataframe exists after running load_data.R file

#run rtsne to get coordinates of extract features
rtsne_out <- Rtsne(as.matrix(gene_data_vs_cell_type[, 1:20270]), pca = FALSE, verbose = TRUE, perplexity = 20)

# (optional, might not need in this file). Changing the cell classes to factors in order to 
# visualize our data. plotting our data. 
gene_data_vs_cell_type$General_Cell_Type<-as.factor(gene_data_vs_cell_type$General_Cell_Type)
gene_data_vs_cell_type$Cell_Type<-as.factor(gene_data_vs_cell_type$Cell_Type)
plot(rtsne_out$Y,col=gene_data_vs_cell_type$General_Cell_Type)
plot(rtsne_out$Y,col=gene_data_vs_cell_type$Cell_Type)

#Heirarchical clustering
hc.complete=hclust(dist(rtsne_out$Y), method="complete")
hc.average=hclust(dist(rtsne_out$Y), method="average")
hc.single=hclust(dist(rtsne_out$Y), method="single")

#(optional, might not need in this file). Plotting our dendrogram.
par(mfrow=c(1,3))
plot(hc.complete,main="Complete Linkage", xlab="", sub="",
     cex =.9)
plot(hc.average , main="Average Linkage", xlab="", sub="",
     cex =.9)
plot(hc.single , main="Single Linkage", xlab="", sub="",
     cex =.9)

#Making a new dataframe that will contain our predicted classes
gene.v.cell.tsne = Gautier_Immgen_Sample_Metadata[,2:3]

#combining classes of three different link functions of HC to our known classes
gene.v.cell.tsne = cbind(gene.v.cell.tsne, cutree(hc.single,8))
gene.v.cell.tsne = cbind(gene.v.cell.tsne, cutree(hc.average,8))
gene.v.cell.tsne = cbind(gene.v.cell.tsne, cutree(hc.complete,8))

#combining X, Y coordinates of rTSNE to our classification
gene.v.cell.tsne = cbind(gene.v.cell.tsne, rtsne_out$Y)

#renaming column names of gene.v.cell.tsne
names(gene.v.cell.tsne)[names(gene.v.cell.tsne) == 'cutree(hc.single, 8)'] <- 's.link.class'
names(gene.v.cell.tsne)[names(gene.v.cell.tsne) == 'cutree(hc.average, 8)'] <- 'a.link.class'
names(gene.v.cell.tsne)[names(gene.v.cell.tsne) == 'cutree(hc.complete, 8)'] <- 'c.link.class'
names(gene.v.cell.tsne)[names(gene.v.cell.tsne) == '1'] <- 'x.tsne'
names(gene.v.cell.tsne)[names(gene.v.cell.tsne) == '2'] <- 'y.tsne'

#storing the mean x.tsne and y.tsne for each class for each link 
s.link.mean <- aggregate(gene.v.cell.tsne[,6:7], list(gene.v.cell.tsne$s.link.class), mean)
a.link.mean <- aggregate(gene.v.cell.tsne[,6:7], list(gene.v.cell.tsne$a.link.class), mean)
c.link.mean <- aggregate(gene.v.cell.tsne[,6:7], list(gene.v.cell.tsne$c.link.class), mean)

#renaming Group.1 column for the mean dataframe for each link
names(s.link.mean)[names(s.link.mean) == 'Group.1'] <- 'class'
names(a.link.mean)[names(a.link.mean) == 'Group.1'] <- 'class'
names(c.link.mean)[names(c.link.mean) == 'Group.1'] <- 'class'

#writing Euclidean distance function
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

#(Code below is very naive, there exists better ways to optimize code)
#example classification of a new testing example s.link
test = c(-3,3)
s.link.dist = c()
for (i in 1: nrow(s.link.mean)){
  s.link.dist <- c(s.link.dist, euc.dist(test, s.link.mean[i,2:3]))
}

#storing the predicted class based on s.link
predicted.class.s.link = which(s.link.dist == min(s.link.dist))

#example classification of a new testing example a.link
test = c(-3,3)
a.link.dist = c()
for (i in 1: nrow(a.link.mean)){
  a.link.dist <- c(a.link.dist, euc.dist(test, a.link.mean[i,2:3]))
}

#storing the predicted class based on a.link
predicted.class.a.link = which(a.link.dist == min(a.link.dist))

#example classification of a new testing example c.link
test = c(-3,3)
c.link.dist = c()
for (i in 1: nrow(c.link.mean)){
  c.link.dist <- c(c.link.dist, euc.dist(test, c.link.mean[i,2:3]))
}

#storing the predicted class based on c.link
predicted.class.c.link = which(c.link.dist == min(c.link.dist))

#NOTE: There are different ways to classify a testing point based on a testing point.
#It is also better to associate a cell type for each class number of each link first, since
#the cell type for class 1 for s.link might be a different cell type for class 1 for a.link.
#We may want to discuss how to classify these as a team before proceeding.