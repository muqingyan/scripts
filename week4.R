library(devtools)
library(Biobase)
library(DESeq2)
library(goseq)
library(MatrixEQTL)
head(supportedGenomes())
head(supportedGeneIDs())
temp_data =read.table(system.file("extdata","Li_sum.txt",
                                  package="goseq"),sep="\t",
                      header=TRUE,
                      stringsAsFactors=FALSE)
expr= temp_data[,-1]#exclude the first row
#or use the row.names =1 to set the rownames for the data.frame in the read.table
rownames(expr) = temp_data[,1]
expr = expr[rowMeans(expr) > 5,]
grp=factor(rep(c("Control","Treated"),times=c(4,3)))
pdata  = data.frame(grp)
#Perform a differential expression analysis
#Now we perform a differential expression analysis for the group variable with DESeq2

de = DESeqDataSetFromMatrix(expr, pdata, ~grp)
de_fit = DESeq(de)
de_results = results(de_fit)
genes <- expr[de_results$padj<0.05,]

#Get the differentially expressed genes after FDR correction

genes2 = as.integer(de_results$padj < 0.05)
not_na = !is.na(genes2)
names(genes2) = rownames(expr)
genes2 = genes2[not_na]
length(genes2[genes2==1])#2702
pwf=nullp(genes2,"hg19","ensGene")#need to adjust for the gene length
#Perform the enrichment analysis parametrically
library(org.Hs.eg.db)
GO.wall=goseq(pwf,"hg19","ensGene")#annotate the genes
#Fetching GO annotations...
#For 2046 genes, we could not find any categories. These genes will be excluded.
#To force their use, please run with use_genes_without_cat=TRUE (see documentation).
#This was the default behavior for version 1.15.1 and earlier.
#Calculating the p-values...
#'select()' returned 1:1 mapping between keys and columns
#Limiting yourself to a single category you are interested in
#Suppose there is a particular category or function you are interested in. You can limit to just that category

GO.MF=goseq(pwf,"hg19","ensGene",test.cats=c("GO:MF"))#in the GO.wall$ontology

head(GO.MF)

library(MatrixEQTL)
library(edgeR)
library(devtools)
library(Biobase)
library(dplyr)
library(snpStats)
library(broom)
#Here we are going to follow along with the tutorial on MatrixEQTL. First we find the files

base.dir = find.package("MatrixEQTL")
SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");
expression_file_name = paste(base.dir, "/data/GE.txt", sep="")
covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="")
output_file_name = tempfile()
#Next we load the data so we can see it

expr = read.table(expression_file_name,sep="\t",
                  header=T,row.names=1)
expr[1,]
snps = read.table(SNP_file_name,sep="\t",
                  header=T,row.names=1)
snps[1,]

cvrt = read.table(covariates_file_name,sep="\t",
                  header=T,row.names=1)
#eQTL is linear regession
#The simplest eQTL analysis just computes linear regression models for each SNP/gene pair.

e1 = as.numeric(expr[1,])
s1 = as.numeric(snps[1,])
lm1 = lm(e1 ~ s1)
tidy(lm1)
set.seed(111)
plot(e1~jitter(s1),col=as.numeric(s1+1))
lines(lm1$fitted.values~s1)
plot(e1 ~ jitter(s1),
     col=(s1+1),xaxt="n",xlab="Genotype",ylab="Expression")
axis(1,at=c(0:2),labels=c("AA","Aa","aa"))
lines(lm1$fitted ~ s1,type="b",pch=15,col="darkgrey")
#Fitting many eQTL models with MatrixEQTL
#Set general parameters
#We need to set up the p-value cutoff and the error model (in this case assuming independent errors)

pvOutputThreshold = 1e-2
errorCovariance = numeric()
useModel = modelLINEAR
#Now we need to set up the snp and gene expression data in the special format required by the MatrixEQTL package

snps = SlicedData$new()
snps$fileDelimiter = "\t"     # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values;
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 2000     # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name )

gene = SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1      # one column of row labels
gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
gene$LoadFile(expression_file_name)

cvrt = SlicedData$new()
#Running MatrixEQTL
#We can now run the code to calculate the eQTL that we are interested in

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = NULL,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)
#We can also figure look at the number and type of eQTL
names(me)
names(me$all)
me$all$neqtls
## [1] 1
me$all$eqtls