suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(MASS)
})
#Initiallizing the data. Assume that load_data.R has already been run.
EMF_GSVA_Relative_Normalilzed_Expression_PATH <- file.path("data", "EMF_GSVA_Relative_Normalilzed_Expression.tsv.tar.gz")
gene_data_vs_cell_type3 <- read_tsv(EMF_GSVA_Relative_Normalilzed_Expression_PATH,
                                    col_types = cols(.default = col_double(), EMF_GSVA_Relative_Normalilzed_Expression.tsv = col_character(),
                                                     Symbol = col_character(), Entrez = col_integer(), Name = col_character())) %>%
  filter(!is.na(Symbol)) %>% dplyr:: select(c(-2,-3,-4))

gene_data_vs_cell_type3 <- t(gene_data_vs_cell_type3)

#making first row as header names:
colnames(gene_data_vs_cell_type3) = gene_data_vs_cell_type3[1, ] # the first row will be the header
gene_data_vs_cell_type3 = gene_data_vs_cell_type3[-1, ]          # removing the first row.
class(gene_data_vs_cell_type3) <- "numeric"              # make values numeric

rownames(gene_data_vs_cell_type3) <- gene_data_vs_cell_type$General_Cell_Type

gene_data_vs_cell_type3 <- t(gene_data_vs_cell_type3)

### #Creating colors and labels for graphing purposes 
gene.v.cell.colors <- c(rep("darkgreen", 3), rep("red", 35), rep("blue",23), rep("brown", 7), rep("orange",34), 
                        rep("darkgray", 15), rep("cyan", 15), rep("violet",3), rep("orange", 2))
gene.v.cell.labels <- c(rep("Mi", 3), rep("Ma", 35), rep("Mo", 23), rep("Ne", 7), rep("De", 34), rep("B", 15)
                        ,rep("T", 15), rep("Nk", 3), rep("De", 2))

### LDA with all classes, with no cross validation. Takes time to run.
lda.fit.all <-lda(t(gene_data_vs_cell_type3),gene.v.cell.labels,CV=FALSE) 
scalings.m <- lda.fit.all$scaling
#NOTE, the absolute values of the scaling is similar to the loadings for PCA.

#Creating new features by LDA for all classes case
new.features.lda.all <- t(scalings.m)%*%gene_data_vs_cell_type3
#Matrix multiplication of Betas(Scalings) and original gene predictors

###Graphing new feature space for LDA all classes, using first two LDA components
x <- as.vector(new.features.lda.all['LD1',])
y<- as.vector(new.features.lda.all['LD2',])
plot(x,y,
     col=gene.v.cell.colors,
     type='n',
     panel.first=grid(col='black'), 
     main="First two LDA components", 
     xlab='LD1', 
     ylab='LD2')
text(x, y,labels=gene.v.cell.labels,col=gene.v.cell.colors,pch=0.5)
#Able to discriminate Neutrophils well

###Graphing new feature space for LDA all classes, using 2nd and 3rd LDA components
x <- as.vector(new.features.lda.all['LD2',])
y<- as.vector(new.features.lda.all['LD3',])
plot(x,y,
     col=gene.v.cell.colors,
     type='n',
     panel.first=grid(col='black'), 
     main="Second Pair LDA components", 
     xlab='LD2', 
     ylab='LD3')
text(x, y,labels=gene.v.cell.labels,col=gene.v.cell.colors,pch=0.5)
#Able to discriminate B cells pretty well

###Graphing new feature space for LDA all classes, using 3rd and 4th LDA components
x <- as.vector(new.features.lda.all['LD3',])
y<- as.vector(new.features.lda.all['LD4',])
plot(x,y,
     col=gene.v.cell.colors,
     type='n',
     panel.first=grid(col='black'), 
     main="Third pair LDA components", 
     xlab='LD3', 
     ylab='LD4')
text(x, y,labels=gene.v.cell.labels,col=gene.v.cell.colors,pch=0.5)
#Able to discriminate Monocytes pretty well

###Graphing new feature space for LDA all classes, using 4th and 5th LDA components
x <- as.vector(new.features.lda.all['LD4',])
y<- as.vector(new.features.lda.all['LD5',])
plot(x,y,
     col=gene.v.cell.colors,
     type='n',
     panel.first=grid(col='black'), 
     main="Fourth Pair LDA components", 
     xlab='LD4', 
     ylab='LD5')
text(x, y,labels=gene.v.cell.labels,col=gene.v.cell.colors,pch=0.5)
#Discriminates Dendric cells fairly well-ish

###Graphing new feature space for LDA all classes, using 5th and 6th LDA components
x <- as.vector(new.features.lda.all['LD5',])
y<- as.vector(new.features.lda.all['LD6',])
plot(x,y,
     col=gene.v.cell.colors,
     type='n',
     panel.first=grid(col='black'), 
     main="Fifth Pair LDA components", 
     xlab='LD5', 
     ylab='LD6')
text(x, y,labels=gene.v.cell.labels,col=gene.v.cell.colors,pch=0.5)
#Discriminate Microglia pretty well 

###Graphing new feature space for LDA all classes, using 6th and 7th LDA components
x <- as.vector(new.features.lda.all['LD6',])
y<- as.vector(new.features.lda.all['LD7',])
plot(x,y,
     col=gene.v.cell.colors,
     type='n',
     panel.first=grid(col='black'), 
     main="Sixth Pair LDA components", 
     xlab='LD6', 
     ylab='LD7')
text(x, y,labels=gene.v.cell.labels,col=gene.v.cell.colors,pch=0.5)
#Discriminate both Microglia and NK cells pretty well.

#None of these determines T cells that well

###### ONE vs. ALL Approaches #######
#Microglia vs. ALL#
group.of.interest <- "Mi"
one.vs.others<- gene.v.cell.labels
one.vs.others[gene.v.cell.labels != group.of.interest] <- "o"
print(table(one.vs.others))
#Correct table

#Making an LDA classifer
one.vs.others.lda.mi <- lda(t(gene_data_vs_cell_type3),one.vs.others,CV=FALSE) 

#Creating dataframe which we will use to look at abs(scalings) to depict significant genes 
LDA.scalings.genes <- as.data.frame(one.vs.others.lda.mi$scaling)
LDA.scalings.genes[,1] <- abs(LDA.scalings.genes[,1])
colnames(LDA.scalings.genes)[1] <- "Mi.LD1"

#Creating new feature matrix by multiple scalings with original features
new.features.lda.mi <- t(one.vs.others.lda.mi$scaling)%*%gene_data_vs_cell_type3

#Plotting this new feature
x <- as.vector(new.features.lda.mi['LD1',])
boxplot(x ~ gene.v.cell.labels, las=1,
        #       col=group.colors,
        horizontal=TRUE,
        main='LD1.Mi', col="#BBBBBB")
### Boxplot separates Microglia very well


######################
##Macrophage vs. ALL##
######################
group.of.interest <- "Ma"
one.vs.others<- gene.v.cell.labels
one.vs.others[gene.v.cell.labels != group.of.interest] <- "o"
print(table(one.vs.others))

#Making an LDA classifer
one.vs.others.lda.ma <- lda(t(gene_data_vs_cell_type3),one.vs.others,CV=FALSE) 

#Adding to the dataframe we looked at before for abs(scalings to depict significant genes)
LDA.scalings.genes <- cbind(LDA.scalings.genes, as.data.frame(one.vs.others.lda.ma$scaling))
LDA.scalings.genes[,2] <- abs(LDA.scalings.genes[,2])
colnames(LDA.scalings.genes)[2] <- "Ma.LD1"

#Creating new feature matrix by multiple scalings with original features
new.features.lda.ma <- t(one.vs.others.lda.ma$scaling)%*%gene_data_vs_cell_type3
#Plotting this new feature
x <- as.vector(new.features.lda.ma['LD1',])
boxplot(x ~ gene.v.cell.labels, las=1,
        #       col=group.colors,
        horizontal=TRUE,
        main='LD1.Ma', col="#BBBBBB")
#Discriminates Macrophages pretty well, despite two false positives which are Dendritic cells 


######################
##Monocytes vs. ALL##
######################
group.of.interest <- "Mo"
one.vs.others<- gene.v.cell.labels
one.vs.others[gene.v.cell.labels != group.of.interest] <- "o"
print(table(one.vs.others))

#Making an LDA classifer
one.vs.others.lda.mo <- lda(t(gene_data_vs_cell_type3),one.vs.others,CV=FALSE) 

#Adding to the dataframe we looked at before for abs(scalings to depict significant genes)
LDA.scalings.genes <- cbind(LDA.scalings.genes, as.data.frame(one.vs.others.lda.mo$scaling))
LDA.scalings.genes[,3] <- abs(LDA.scalings.genes[,3])
colnames(LDA.scalings.genes)[3] <- "Mo.LD1"

#Creating new feature matrix by multiple scalings with original features
new.features.lda.mo <- t(one.vs.others.lda.mo$scaling)%*%gene_data_vs_cell_type3
#Plotting this new feature
x <- as.vector(new.features.lda.mo['LD1',])
boxplot(x ~ gene.v.cell.labels, las=1,
        #       col=group.colors,
        horizontal=TRUE,
        main='LD1.Mo', col="#BBBBBB")
#Discriminates Monocytes pretty well


######################
##Neutrophils vs. ALL##
######################
group.of.interest <- "Ne"
one.vs.others<- gene.v.cell.labels
one.vs.others[gene.v.cell.labels != group.of.interest] <- "o"
print(table(one.vs.others))

#Making an LDA classifer
one.vs.others.lda.ne <- lda(t(gene_data_vs_cell_type3),one.vs.others,CV=FALSE) 

#Adding to the dataframe we looked at before for abs(scalings to depict significant genes)
LDA.scalings.genes <- cbind(LDA.scalings.genes, as.data.frame(one.vs.others.lda.ne$scaling))
LDA.scalings.genes[,4] <- abs(LDA.scalings.genes[,4])
colnames(LDA.scalings.genes)[4] <- "Ne.LD1"

#Creating new feature matrix by multiple scalings with original features
new.features.lda.ne <- t(one.vs.others.lda.ne$scaling)%*%gene_data_vs_cell_type3
#Plotting this new feature
x <- as.vector(new.features.lda.ne['LD1',])
boxplot(x ~ gene.v.cell.labels, las=1,
        #       col=group.colors,
        horizontal=TRUE,
        main='LD1.Ne', col="#BBBBBB")
#Discriminates Neutrophils very well



######################
##Dendritic Cells vs. ALL##
######################
group.of.interest <- "De"
one.vs.others<- gene.v.cell.labels
one.vs.others[gene.v.cell.labels != group.of.interest] <- "o"
print(table(one.vs.others))

#Making an LDA classifer
one.vs.others.lda.de <- lda(t(gene_data_vs_cell_type3),one.vs.others,CV=FALSE) 

#Adding to the dataframe we looked at before for abs(scalings to depict significant genes)
LDA.scalings.genes <- cbind(LDA.scalings.genes, as.data.frame(one.vs.others.lda.de$scaling))
LDA.scalings.genes[,5] <- abs(LDA.scalings.genes[,5])
colnames(LDA.scalings.genes)[5] <- "De.LD1"

#Creating new feature matrix by multiple scalings with original features
new.features.lda.de <- t(one.vs.others.lda.de$scaling)%*%gene_data_vs_cell_type3
#Plotting this new feature
x <- as.vector(new.features.lda.de['LD1',])
boxplot(x ~ gene.v.cell.labels, las=1,
        #       col=group.colors,
        horizontal=TRUE,
        main='LD1.De', col="#BBBBBB")
#Discriminates Dendritic Cells pretty well, but there are three outlier Macrophages. 
#Makes sense because Macrophages had some outlier Dendritic Cells. 



######################
## B Cells vs. ALL##
######################
group.of.interest <- "B"
one.vs.others<- gene.v.cell.labels
one.vs.others[gene.v.cell.labels != group.of.interest] <- "o"
print(table(one.vs.others))

#Making an LDA classifer
one.vs.others.lda.b <- lda(t(gene_data_vs_cell_type3),one.vs.others,CV=FALSE) 

#Adding to the dataframe we looked at before for abs(scalings to depict significant genes)
LDA.scalings.genes <- cbind(LDA.scalings.genes, as.data.frame(one.vs.others.lda.b$scaling))
LDA.scalings.genes[,6] <- abs(LDA.scalings.genes[,6])
colnames(LDA.scalings.genes)[6] <- "B.LD1"

#Creating new feature matrix by multiple scalings with original features
new.features.lda.b <- t(one.vs.others.lda.b$scaling)%*%gene_data_vs_cell_type3
#Plotting this new feature
x <- as.vector(new.features.lda.b['LD1',])
boxplot(x ~ gene.v.cell.labels, las=1,
        #       col=group.colors,
        horizontal=TRUE,
        main='LD1.B', col="#BBBBBB")
#Discriminates B Cells pretty well


######################
## T Cells vs. ALL##
######################
group.of.interest <- "T"
one.vs.others<- gene.v.cell.labels
one.vs.others[gene.v.cell.labels != group.of.interest] <- "o"
print(table(one.vs.others))

#Making an LDA classifer
one.vs.others.lda.t <- lda(t(gene_data_vs_cell_type3),one.vs.others,CV=FALSE) 

#Adding to the dataframe we looked at before for abs(scalings to depict significant genes)
LDA.scalings.genes <- cbind(LDA.scalings.genes, as.data.frame(one.vs.others.lda.t$scaling))
LDA.scalings.genes[,7] <- abs(LDA.scalings.genes[,7])
colnames(LDA.scalings.genes)[7] <- "T.LD1"

#Creating new feature matrix by multiple scalings with original features
new.features.lda.t <- t(one.vs.others.lda.t$scaling)%*%gene_data_vs_cell_type3
#Plotting this new feature
x <- as.vector(new.features.lda.t['LD1',])
boxplot(x ~ gene.v.cell.labels, las=1,
        #       col=group.colors,
        horizontal=TRUE,
        main='LD1.T', col="#BBBBBB")
#Does not discriminate T from NK cells that well. Might be an issue.

######################
## NK Cells vs. ALL##
######################
group.of.interest <- "Nk"
one.vs.others<- gene.v.cell.labels
one.vs.others[gene.v.cell.labels != group.of.interest] <- "o"
print(table(one.vs.others))

#Making an LDA classifer
one.vs.others.lda.nk <- lda(t(gene_data_vs_cell_type3),one.vs.others,CV=FALSE) 

#Adding to the dataframe we looked at before for abs(scalings to depict significant genes)
LDA.scalings.genes <- cbind(LDA.scalings.genes, as.data.frame(one.vs.others.lda.nk$scaling))
LDA.scalings.genes[,8] <- abs(LDA.scalings.genes[,8])
colnames(LDA.scalings.genes)[8] <- "NK.LD1"

#Creating new feature matrix by multiple scalings with original features
new.features.lda.nk <- t(one.vs.others.lda.nk$scaling)%*%gene_data_vs_cell_type3
#Plotting this new feature
x <- as.vector(new.features.lda.nk['LD1',])
boxplot(x ~ gene.v.cell.labels, las=1,
        #       col=group.colors,
        horizontal=TRUE,
        main='LD1.NK', col="#BBBBBB")
#Surprising NK cells is discriminated well, even against T cells

write.csv(LDA.scalings.genes, "LDA_scalings_genes.csv", row.names = TRUE)


#### EXAMPLE OF OVERFITTING OUR DATA ###
#Because we do not have a lot of observations for certain cell types
#it is very easy to overfit our data. Also we have a very imbalanced data set
#Below shows an example of overfitting microglia

#IN this case, we will get 100% prediction accuracy for Microglia if we
#use training data and testing data as the same for LDA classifications

#Microglia vs. ALL#
group.of.interest <- "Mi"
one.vs.others<- gene.v.cell.labels
one.vs.others[gene.v.cell.labels != group.of.interest] <- "o"
print(table(one.vs.others))
#Correct table

#We already made a classifier for microglia
#one.vs.others.lda.mi <- lda(t(gene_data_vs_cell_type3),one.vs.others,CV=FALSE) 

predict.lda.mi <- predict(object=one.vs.others.lda.mi, newdata = t(gene_data_vs_cell_type3))

## Posterior probabilities (print the 10 first rows only)
print(predict.lda.mi$posterior[1:10,])

## Predicted classes
print(predict.lda.mi$class)

## Print the contingency table
table(one.vs.others, predict.lda.mi$class)

## Compute the hit rate
hits <- sum(one.vs.others == predict.lda.mi$class)
errors <- sum(one.vs.others != predict.lda.mi$class)
total <- hits + errors
(hit.rate.internal <- hits / total)
#100% hit rate

#Compute the error rate
(error.rate.internal <- errors / total)
#0% Error Rate

###USING a LOO cross-validation approach for Mi###
## Run a Leave-one-out (LOO) cross-validation
one.vs.others.lda.loo.mi <- lda(t(gene_data_vs_cell_type3),one.vs.others,CV=TRUE) 

table(one.vs.others, one.vs.others.lda.loo.mi$class)

## Compute the hit rate
(hit.rate.loo.mi <- sum(one.vs.others == one.vs.others.lda.loo.mi$class) / total)
###Only 55.47% now 
(error.rate.loo.mi <- 1 - hit.rate.loo.mi)
###Error rate is 44.52%

##Compute the random hit rate ###
## Run 10000 permutation tests to estimate the random expectation for the hit rate
random.hit.rates <- vector()
for (rep in 1:10000) {
  random.hit.rates <- append(random.hit.rates, sum(one.vs.others == sample(one.vs.others)) / total)
}
(random.hit.rates.mean <- mean(random.hit.rates))

## Compute the theoretical value for the random expectation
prior <- as.vector(table(one.vs.others))/length(one.vs.others)
(hit.rate.expect <- sum(prior^2))
#Theoretical baseline hit rate should be 95.71634% 

## Draw the histogram
hist(random.hit.rates, breaks=(0:total)/total, col="lightgrey", 
     freq=TRUE,
     main="Hit rate analysis",
     xlab="Hit rate",
     ylab="Frequency")
arrows(x0=hit.rate.loo.mi, y0 = 1000, x1=hit.rate.loo.mi, y1=100, 
       col="darkgreen", lwd=2, code=2, , length=0.2, angle=20)
arrows(x0=hit.rate.internal, y0 = 1000, x1=hit.rate.internal, y1=100, 
       col="red", lwd=2, code=2, , length=0.2, angle=20)
arrows(x0=hit.rate.expect, y0 = 1000, x1=hit.rate.expect, y1=100, 
       col="darkblue", lwd=2, code=2, , length=0.2, angle=20)

legend("topleft", legend=c("random", "random expectation", "LOO", "internal"), 
       col=c("grey", "darkblue", "darkgreen", "red"), 
       lwd=c(4,2,2,2))
