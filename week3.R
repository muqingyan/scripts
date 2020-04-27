library(devtools)
library(Biobase)
library(snpStats)
library(broom)
library(MASS)
library(DESeq2)
library(DESeq)
library(dplyr)
#Logistic regression
#Load the data
#Here we use example SNP data from a case-control genome-wide association study Load an example data set and take a smaller subset of samples for computational efficiency
par(pch=19)
data(for.exercise)#1000 samples and 28501 snps
use <- seq(1, ncol(snps.10), 10)#subset 2851 snps
sub.10 <- snps.10[,use]#2851 snps for 1000 samples
#Calculate the PCs

xxmat <- xxt(sub.10, correct.for.missing=FALSE)
evv <- eigen(xxmat, symmetric=TRUE)
pcs <- evv$vectors[,1:5]
#A single logistic regression
#First we do an unadjusted logistic regression assuming an additive model.The coefficient is the change in log-odds for a one unit decrease (because homozygous major allele is coded 1) in the number of copies of the minor allele.

snpdata = sub.10@.Data
status = subject.support$cc#get the status for all the 1000 samples
snp1 = as.numeric(snpdata[,1])#get the snp data for all 1000 peaple for the first SNP
snp1[snp1==0] = NA#zero means NA
glm1 = glm(status ~ snp1,family="binomial")
tidy(glm1)
#We can also easily code in other models. For example suppose we want to code a dominant model (so only an association of risk with the two copies of the common allele, now the coefficient on snp1_dom is the increase in log odds associated with two copies of the major allele).

snp1_dom = (snp1 == 1)
glm1_dom = glm(status ~ snp1_dom,family="binomial")
tidy(glm1_dom)
#We can also easily adjust for other variables.

glm2 = glm(status ~ snp1 + pcs[,1:5],family="binomial")
tidy(glm2)
#Fit many glms at once
#For logistic regression modeling of many SNPs at once we can use the snps.rhs.tests function which computes an asymptotic chi-squared statistic. This isn’t quite the same thing as the F-statistics we have been calculating but can be used in the same way for significance calculations.

glm_all = snp.rhs.tests(status ~ 1,snp.data=sub.10)# ~1 menas using the unadjusted model
slotNames(glm_all)
## [1] "snp.names" "var.names" "chisq"     "df"        "N"
qq.chisq(chi.squared(glm_all),df=1)
#We can also adjust for variables like principal components

glm_all_adj = snp.rhs.tests(status ~ pcs,snp.data=sub.10)
qq.chisq(chi.squared(glm_all_adj),df=1)


#Poisson/negative binomial regression
#that is a comparative RNA-seq analysis of different mouse strains.

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))
fdata = fData(bot)
edata_3 = edata[rowMeans(edata) >10,]
edata_2 <- filter(as.data.frame(edata),rowMeans(edata) >10)
edata_2 <- as.matrix(edata_2)
#A single Poisson regression
#The coefficient in this case is the increase in the log number of counts comparing one strain to the other. If you exponentiate the coefficients then exp(intercept) is expected count for C57BL/6J and exp(intercept)exp(strainDBA/2J) is the expected count for the strain DBA/2J.

glm3 = glm(edata[1, ] ~ pdata$strain,family="poisson")
tidy(glm3)
#You can also fit a negative binomial regression one at a time in R with:
  
  glm.nb1 = glm.nb(edata[1, ] ~ pdata$strain)
tidy(glm.nb1)
#Multiple negative binomial regressions
#We can use the DESeq2 package to perform many (moderated) negative binomial regressions at once. We first need to create a DESeq data set.

de = DESeqDataSetFromMatrix(edata_3, pdata, ~strain)
glm_all_nb = DESeq(de)
result_nb = results(glm_all_nb)
hist(result_nb$stat)

library(bladderbatch)
library(sva)
data(bladderdata)
pheno = pData(bladderEset)
edata = exprs(bladderEset)
#Plot data with two different model fits
## Setting seed so the jitter will be the same
set.seed(123)

cancerjit = jitter(as.numeric(pheno$cancer))
lm1 = lm(edata[1,] ~ 1)
lm2 = lm(edata[1,] ~ pheno$cancer)

par(mfrow=c(1,2))

plot(edata[1,] ~ cancerjit,type="p",xaxt="n",xlab="Cancer Status",ylab="Expression",pch=19,col=as.numeric(pheno$cancer))
axis(1,at=c(1,2,3),c("Biopsy","Cancer","Normal"))
abline(lm1,col="darkgrey",lwd=5)


plot(edata[1,] ~ cancerjit,type="p",xaxt="n",xlab="Cancer Status",ylab="Expression",pch=19,col=as.numeric(pheno$cancer))
axis(1,at=c(1,2,3),c("Biopsy","Cancer","Normal"))
boxplot(lm2$fitted~pheno$cancer,add=T,border=1:3)

#Plot the residuals
par(mfrow=c(1,2))
plot(lm1$residuals ~ cancerjit,type="p",xaxt="n",xlab="Cancer Status",ylab="Expression",pch=19,col=as.numeric(pheno$cancer),ylim=c(-1.1,1.1))
axis(1,at=c(1,2,3),c("Biopsy","Cancer","Normal"))

plot(lm2$residuals ~ cancerjit,type="p",xaxt="n",xlab="Cancer Status",ylab="Residuals",pch=19,col=as.numeric(pheno$cancer),ylim=c(-1.1,1.1))
axis(1,at=c(1,2,3),c("Biopsy","Cancer","Normal"))

edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]

#Calculate t- or F-statistics rapidly
#The genefilter package lets you compute statistics rapidly for very simple cases (two or multi-group) comparisons. These are not moderated in any way.

tstats_obj = rowttests(edata,pdata$strain)
names(tstats_obj)
hist(tstats_obj$statistic,col=2)

#For multi-group comparisons use the F-statistic

fstats_obj = rowFtests(edata,as.factor(pdata$lane.number))
names(fstats_obj)
#Fit many statistics with limma
#This approach fits many moderated statistics simultaneously

mod = model.matrix(~ pdata$strain)
fit_limma = lmFit(edata,mod)
ebayes_limma = eBayes(fit_limma)
head(ebayes_limma$t)
plot(ebayes_limma$t[,2],-tstats_obj$statistic,col=4,
     xlab="Moderated T-stat",ylab="T-stat")
abline(c(0,1),col="darkgrey",lwd=3)
#Fit many adjusted statistics with limma
#Here we adjust for the lane number, now the test-statistic for the strain is adjusted for the lane number (a surrogate for a batch effect).

mod_adj = model.matrix(~ pdata$strain + as.factor(pdata$lane.number))
fit_limma_adj = lmFit(edata,mod_adj)
ebayes_limma_adj = eBayes(fit_limma_adj)
head(ebayes_limma_adj$t)

plot(ebayes_limma_adj$t[,2],-tstats_obj$statistic,col=4,
     xlab="Moderated T-stat",ylab="T-stat")
abline(c(0,1),col="darkgrey",lwd=3)
#Calculating a nested model comparison with limma
#Sometimes we want to compare the null model to the alternative model with some additional covariates. Here we have to know which coefficients we want to test in the alternative model.

#Suppose we wanted to find lane effects then we can fit a limma model and find which coefficients belong to the lane variable.

mod_lane = model.matrix(~ as.factor(pdata$lane.number))
fit_limma_lane = lmFit(edata,mod_lane)
ebayes_limma_lane = eBayes(fit_limma_lane) 
head(ebayes_limma_lane$t)


top_lane = topTable(ebayes_limma_lane, coef=2:7,
                    number=dim(edata)[1],sort.by="none")
head(top_lane)

plot(top_lane$F,fstats_obj$statistic,
     xlab="Moderated F-statistic",ylab="F-statistic",col=3)
#Calculating a nested comparison with edge
#We can also perform the unmoderated nested comparisons in the edge package, which also has functions for calculating a more powerful odp statistic

edge_study = build_study(edata, grp = as.factor(pdata$lane.number))
de_obj = lrt(edge_study)
qval = qvalueObj(de_obj)
plot(qval$stat,fstats_obj$statistic,col=4,
     xlab="F-stat from edge",ylab="F-stat from genefilter")

#We can easily adjust for variables by passing arguments to the adj.var variable.

edge_study2 = build_study(edata, grp = as.factor(pdata$lane.number),
                          adj.var=pdata$strain)
de_obj2 = lrt(edge_study2)
qval2 = qvalueObj(de_obj2)
plot(qval2$stat,fstats_obj$statistic,col=4,
     xlab="F-stat from edge",ylab="F-stat from genefilter")

library(genefilter)
edata <- log2(edata+1)
edata <- edata[rowMeans(edata)>10,]
tstats_obj <- rowttests(edata,pdata$strain)
#We can now permute the sample labels using the sample function in R.

set.seed(135)
strain = pdata$strain
strain0 = sample(strain)
tstats_obj0 = rowttests(edata,strain0)
hist(tstats_obj0$statistic,col=2,xlim=c(-5,2))
#We can now compare the observed and permuted statistics

quantile(tstats_obj0$statistic)
##         0%        25%        50%        75%       100% 
## -0.9629544  0.7488187  1.2858778  1.7910909  4.2376976
quantile(tstats_obj$statistic)

library(qvalue)
hist(tstats_obj$p.value)

#Adjusting for variables with edge
#If you want to adjust for variables you need to use edge

edge_study = build_study(edata, grp = pdata$strain, 
                         adj.var = as.factor(pdata$lane.number))
de_obj = lrt(edge_study)
qval = qvalueObj(de_obj)
hist(qval$pvalues,col=3)

#Multiple testing
#To correct for multiple testing you can use the Bonferroni correction or different FDR corrections.

#Bonferroni and Benjamini-Hochberg FDR correction with p.adjust
#You can use the p.adjust function to get “multiple testing corrected” p-values which you can then use to control error rates.
fstats_obj <- rowFtests(edata,pdata$strain)
fp_bonf = p.adjust(fstats_obj$p.value,method="bonferroni")
hist(fp_bonf,col=3)

