library(devtools)
library(Biobase)
library(genstats)
#> BiocManager::install("jtleek/genstats",ref="gh-pages")

con=url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset#an expression set
#table1: genomics data
exp_data = exprs(bm)
dim(exp_data)
#table2: the phenotype data
pheno_data = pData(bm)
dim(pheno_data)
#table3: the feature data
feature_data = fData(bm)
dim(feature_data)

#randomization is a good practice to address the confounding variable, to break the confounding
#relationships with the treatments
# it is better to have balanced experimental design
#data exploratory
install_github('andreacirilloac/updateR')
library(updateR)
updateR(admin_password = '032309')#update R to 3.6

install.packages("caTools")
library(gplots)
library(devtools)
library(Biobase)
library(RSkittleBrewer)#devtools::install_github('alyssafrazee/RSkittleBrewer')
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)

trop = RSkittleBrewer("tropical")#choose the colors
palette(trop)#Load the library and set the color palette with the palette function. Now when I type col = 1 it will look for the first color in the trop colors. 
#We also set the character to be a filled dot with par(pch=19)
par(pch=19)
#load data, bm = bodymap.eset
table(pheno_data$gender)
table(pData(bm)$gender,pData(bm)$race)
summary(exp_data)
# Use option useNA to include NA's in table
table(pheno_data$age,useNA="ifany")
# is.na checks for NA values
table(is.na(pheno_data$age))
# Check for other common missing names
sum(pheno_data$age==" ")
gene_na = is.na(exp_data)
table(rowSums(gene_na))

'po[k˚˚l'
boxplot(log2(exp_data+1),col=1,range =0)
par(mfrow=c(2,3))#arrange plots into 2 rows and 2 columns
hist(log2(exp_data[,1]+1),col=1)
hist(log2(exp_data[,2]+1),col=2)
hist(log2(exp_data[,3]+1),col=3)
hist(log2(exp_data[,4]+1),col=4)
hist(log2(exp_data[,5]+1),col=5)
#Or with density plots
edata <- exp_data
plot(density(log2(edata[,1]+1)),col=2,ylim=c(0,3))
lines(density(log2(edata[,2]+1)),col=3)#use line() to overlay the plots
lines(density(log2(edata[,3]+1)),col=4)
for(i in 4:100){lines(density(log2(edata[,i]+1)),col=i,ylim=c(0,3))}
#A very common task is to compare distributions of measurements (say before normalization). You can do this with a qq-plot
qqplot(log2(edata[,1]+1), log2(edata[,2]+1),col=3,pch = 19)#how the two distributions compares t each other, the quantile plot
abline(c(0,1))
#A very widely used plot is what is known as a M-A plot, 
#sometimes called a Bland Altman plot. 
#The basic idea is to plot the sum of the two values on the x-axis and the difference on the y-axis. 
#This can be used to see any difference between the (samples, averages, etc.) and to see if there is any intensity-specific biases.
sum <- log2(edata[,1]+1) + log2(edata[,2]+1)
diff <- log2(edata[,1]+1) - log2(edata[,2]+1)
plot(sum,diff,col=2,pch=19)

edata = as.data.frame(exp_data)#package dplyr
filt_edata = filter(edata,rowMeans(edata) > 1)
filt_edata = edata[rowMeans(edata)>1,]
head(filt_edata)
boxplot(as.matrix(log2(filt_edata+1)),col=2,pch = 19)
boxplot(log2(filt_edata+1),col=3,pch = 19)
#check data mixups
#Get the chromosomes for each gene using the feature data.

aeid = as.character(fdata[,1])
chr = AnnotationDbi::select(org.Hs.eg.db,keys=aeid,keytype="ENSEMBL",columns="CHR")
head(chr)
dim(chr)#52607
dim(edata)#52580
#features have been mapped to more than two chromosomes,to deduplicate
# Take non-duplicated chromsomes
chr = chr[!duplicated(chr[,1]),]

# Confirm that the annotation still is in the right order
all(chr[,1] == rownames(edata))
# Select the chromosome Y samples
edatay = dplyr::filter(edata,chr$CHR=="Y")
# Males have Y chromsome expression as expected
boxplot(colSums(edatay) ~ pdata$gender)
points(colSums(edatay) ~ jitter(as.numeric(pdata$gender)),#overlay the datapoints,jitter means the dots don't overlay each other
       col=as.numeric(pdata$gender),
       pch=19)
#heatmaps
ematrix = as.matrix(edata)[rowMeans(edata) > 10000,]
heatmap(ematrix)
dim(ematrix)
colramp = colorRampPalette(c(3,"white",2))(9)
colramp = colorRampPalette(c("navy","white","firebrick3"))(90)
heatmap(ematrix,col=colramp)
heatmap(ematrix, col = colorRampPalette(c( "white","navy","firebrick3"))(50), fontsize=9, fontsize_row=6) #自定义颜色
#You might have noticed some automatic clustering here, you can turn that off (we’ll learn more about it in a later lecture)
heatmap(ematrix,col=colramp,Rowv=NA,Colv=NA)
#If you load the gplots package you can add a color scale with the heatmap.2 package. Here we have to add some options to make the dendogram disappear, scale the data by rows, and remove a tracing plot.
heatmap.2(ematrix,col=colramp,scale="row",trace="none")
#or remove the dendrogram
heatmap.2(ematrix,col=colramp,
          dendrogram="none", scale="row",trace="none")
hist(log2(edata[,1]+1),col=2,breaks = 100,xlim = c(1,15),ylim = c(0,400))
hist(rowSums(edata==0),xlim = c(0,20))
high_set <- filter(edata,rowMeans(edata) >=5)#object should be data.frame
#k-means clustering sets the starting values randomly
library(dendextend)
#First we log transform and remove lowly expressed genes, then calculate Euclidean distance
high_exp <- filter(as.data.frame(edata),rowMeans(edata) >5000)
#or 
high_exp2<- edata[rowMeans(edata)>5000,]
high_exp <- log2(high_exp+1)
# By default calculates the distance between rows,so use the t() to transpose
dist1 = dist(t(high_exp))
## Look at distance matrix
colramp = colorRampPalette(c(3,"white",2))(9)
heatmap(as.matrix(dist1),col=colramp,Colv=NA,Rowv=NA)
#Now cluster the samples
#Here we use the distance we previously calculated to perform a hierarchical clustering and plot the dendrogram:
hclust1 = hclust(dist1)
plot(hclust1)
#We can also force all of the leaves to terminate at the same spot
plot(hclust1,hang=-1)
#We can also color the dendrogram either into a fixed number of groups
dend = as.dendrogram(hclust1)
dend = color_labels(hclust1,4,col=1:4)#split into 4 clusters
plot(dend)
#Or you can color them directly
labels_colors(dend) = c(rep(1,10),rep(2,9))
plot(dend)
#Now we can perform k-means clustering. By default, the rows are clustered. You can either input the cluster means (often unknown) or the number of clusters (here I input 3 clusters)
kmeans1 = kmeans(high_exp,centers=3)
names(kmeans1)
par(mfrow=c(1,1))
par()
matplot(t(kmeans1$centers),col=c(3,1,2),type="l",lwd=3)
table(kmeans1$cluster)
heatmap(as.matrix(high_exp)[order(kmeans1$cluster),],col=colramp,Colv=NA,Rowv=NA)
dist1

kmeans2 = kmeans(t(high_exp),centers=3)
names(kmeans2)
par(mfrow=c(1,1))
par()
matplot(kmeans2$centers,col=1:3,type="l",lwd=3)
table(kmeans2$cluster)
heatmap(as.matrix(high_exp)[order(kmeans2$cluster),],col=colramp,Colv=NA,Rowv=NA)


