library(devtools)
library(Biobase)
colors <- c("blue2","firebrick3","green4","darkorange2","gray")
palette(colors)
matr <- matrix(1:10,ncol=2,byrow = T)
for(i in 1:5){plot(matr[,1],matr[,2],col =i)}
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
#Here we calculate the singular vectors:
edata = edata[rowMeans(edata) > 100, ]#filter low gene expression
edata = log2(edata + 1)#transform
edata_centered = edata - rowMeans(edata)#need to remove the mean variation
svd1 = svd(edata_centered)
names(svd1)
#D is the diagnal matrix == ncol(edata); v variance across genes;u across samples
#svd: singular value decomposition
#Look at the percent variance explained
#The percent of variance explained is given by dii^2∑jd^2jj
plot(svd1$d,ylab="Singular value",col=2)
plot(svd1$d^2/sum(svd1$d^2),ylab="Percent Variance Explained",col="navy",xlim=c(1,3))
par(pch=19)
#Plot top two principal components
par(mfrow=c(1,2))
plot(svd1$v[,1],col=2,ylab="1st PC")
plot(svd1$v[,2],col=2,ylab="2nd PC")
#Plot PC1 vs. PC2
plot(svd1$v[,1],svd1$v[,2],ylab="2nd PC",
     xlab="1st PC",col=as.numeric(pdata$study))

#Another common plot is to make boxplots comparing the PC for different levels of known covariates (don’t forget to show the actual data!).
boxplot(svd1$v[,1] ~ pdata$study,border=c(1,2))
points(svd1$v[,1] ~ jitter(as.numeric(pdata$study)),col=as.numeric(pdata$study))
#What we have been plotting is not exactly the principal components.
pc1 = prcomp(edata)
plot(pc1$rotation[,1],svd1$v[,1])
#To get the actual PCs you have to subtract the column means rather than the row means when normalizing.
edata_centered2 = t(t(edata) - colMeans(edata))
svd2 = svd(edata_centered2)
plot(pc1$rotation[,1],svd2$v[,1],col=2)
plot(pc1$rotation[,1],pc1$rotation[,2],col=as.numeric(pdata$study))
plot(svd2$v[,1],svd2$v[,2],ylab="2nd PC",
     xlab="1st PC",col=as.numeric(pdata$study))
#What happens if we introduce a single outlying gene

edata_outlier = edata_centered# edata-rowMeans(edata)
edata_outlier[6,] = edata_centered[6,] * 10000# introduce one outlier
svd3 = svd(edata_outlier)
par(mfrow=c(1,2))
plot(svd1$v[,1],col=1,main="Without outlier")
plot(svd3$v[,1],col=2,main="With outlier")
plot(svd1$v[,1],svd3$v[,1],col=3)
#It turns out the new top singular vector is perfectly correlated with the outlying gene

plot(svd3$v[,1],edata_outlier[6,],col=4)
#preprocessing and normalization
#quantile normalization forces the distribution among samples to be the same
#if the variance within one group is high then quantile normalization should be used to address the techical issues

library(preprocessCore)
#show distributions for counts across samples
edata = log2(edata + 1)
edata = edata[rowMeans(edata) > 3, ]
colramp = colorRampPalette(c(2,"white",1))(20)
plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:20){lines(density(edata[,i]),lwd=3,col=colramp[i])}
#Now we perform quantile normalization to make the distributions the same across samples. Note that near the tail the distributions aren’t perfectly the same, but for the most part the distributions land right on top of each other.
norm_edata = normalize.quantiles(as.matrix(edata))
plot(density(norm_edata[,1]),col=colramp[1],lwd=2,ylim=c(0,0.3))
plot(density(norm_edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.20))
for(i in 2:20){lines(density(norm_edata[,i]),lwd=3,col=colramp[i])}
#Normalization removes bulk differences due to technology. But there still may be differences you don’t want after normalization. The only way to figure this out is to check. For example if we plot the quantile normalized data with the first
plot(norm_edata[1,],col=as.numeric(pdata$study))
boxplot(colSums(norm_edata) ~ pdata$study,col=as.numeric(pdata$study),range = 0)
points(colSums(norm_edata) ~ jitter(as.numeric(pdata$study)),#overlay the datapoints,jitter means the dots don't overlay each other
       col=as.numeric(pdata$study),
       pch=19)
svd4 <-svd(norm_edata - rowMeans(norm_edata))
plot(svd4$v[,1],svd4$v[,2],col=as.numeric(pdata$study), xlab = "PC1",ylab="PC2")

pc2 <- prcomp(norm_edata)
plot(pc2$rotation[,1],pc2$rotation[,2],col=as.numeric(pdata$study))

svd4_2 <- svd(t(t(norm_edata)-colMeans(norm_edata)))
plot(svd4_2$v[,1],svd4_2$v[,2],col=as.numeric(pdata$study))
# before norm
svd4_3 <- svd(t(t(edata)-colMeans(edata)))
plot(svd4_3$v[,1],svd4_3$v[,2],col=as.numeric(pdata$study))

library(UsingR)
data("galton")
par(mfrow = c(1,1))
hist(galton$child,breaks = 80,col=2)
hist(galton$parent,breaks = 80,col=1)
meanchild <- mean(galton$child)
lines(rep(meanchild,200),seq(0,200,length=200),col=3,lwd =5)
plot(galton$parent,galton$child,col=2)
near65 <- galton[abs(galton$parent - 65)<1, ]
points(near65$parent,near65$child,pch=19,col="red")
lines(seq(64,66,length=100),rep(mean(near65$child),100),col="red",lwd=4)
points(near65$parent,near65$child,col=3)
lines(seq(64,66,length=2),rep(meanchild,2),col=3,lwd=5)
#Fitting a line
plot(galton$parent,galton$child,pch=19,col=2)
lm1 <- lm(galton$child ~ galton$parent)
lines(galton$parent,lm1$fitted,col=3,lwd=3)
#all>=70
hi70_p <- galton[galton$parent>=70,]
hi70 <- hi70_p[hi70_p$child>=70,]
plot(galton$parent,galton$child,col="blue",pch=19)
points(hi70$parent,hi70$child,col="green",pch=2)
lines(seq(63,73,length=10),rep(70,10),col="red",pch=8,lwd=6)
lines(rep(70,10),seq(68,75,length=10),col="red",lwd=6)
#Plot what is leftover,aka residuals
par(mfrow=c(1,2))
plot(galton$parent,galton$child,pch=19,col="blue")
lines(galton$parent,lm1$fitted,col="red",lwd=3)
plot(galton$parent,lm1$residuals,col="blue",pch=19)
abline(h=0,col="red",lwd=3)
abline(v=mean(galton$parent),col="red",lwd=3,lty=3)

hunger_ori = read.csv("http://apps.who.int/gho/athena/data/GHO/WHOSIS_000008.csv?profile=text&filter=COUNTRY:*;SEX:*")
hunger <- hunger[hunger$Sex!="Both sexes",]
plot(hunger$Year,hunger$Numeric,col=2,pch =19)
head(hunger)
lm2 <-lm(hunger$Numeric ~ hunger$Year)
lines(hunger$Year,lm2$fitted.values,col=3)
hunger = read.csv("http://apps.who.int/gho/athena/data/GHO/WHOSIS_000008.csv?profile=text&filter=COUNTRY:*;SEX:*")
hunger <- hunger[hunger$Sex!="Both sexes",]
hunger <- hunger[!grepl("-",hunger$Year),]
hunger <- hunger[as.matrix(hunger$Year)>=1991,]
head(hunger)
lm1 <- lm(hunger$Numeric ~ hunger$Year)
plot(hunger$Year,hunger$Numeric,pch=19,col="darkgrey")
points(hunger$Year,hunger$Numeric,pch=19,col=((hunger$Sex=="Male")*1+1))
points(hunger$Year,hunger$Numeric,col=as.numeric(hunger$Sex))
lmM <- lm(hunger$Numeric[hunger$Sex=="Male"] ~ hunger$Year[hunger$Sex=="Male"])
lines(hunger$Year[hunger$Sex=="Male"],lmM$fitted,col="black",lwd=3)
lines(hunger$Year[hunger$Sex=="Female"],lmF$fitted,col="red",lwd=3)


lmBoth <- lm(hunger$Numeric ~ hunger$Year + hunger$Sex)
plot(hunger$Year,hunger$Numeric,pch=19)
points(hunger$Year,hunger$Numeric,pch=19,col=((hunger$Sex=="Male")*1+1))
abline(c(lmBoth$coeff[1],lmBoth$coeff[2]),col="red",lwd=3)
abline(c(lmBoth$coeff[1] + lmBoth$coeff[3],lmBoth$coeff[2] ),col="black",lwd=3)

library(broom)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata=pData(bm)
edata=as.data.frame(exprs(bm))
edata <- edata[,!is.na(pdata$age)]
edata = as.matrix(edata)
fdata = fData(bm)
pdata_age <- as.matrix(filter(as.data.frame(pdata$age),!is.na(pdata$age)))
#Fit a simple linear regression
#Here we regress the counts for the first gene on the age of the people from whom the tissue was drawn.
lm1 <- lm(as.matrix(edata)[1,] ~ as.matrix(pdata_age)) 
plot(pdata_age,edata[1,])
lines(pdata_age,lm1$fitted.values)
abline(lm1$coefficients[1],lm1$coefficients[2],col=2)
edata = as.matrix(edata)
lm1 = lm(edata[1,] ~ pdata$age)
tidy(lm1)
plot(pdata$age,edata[1,], col=1)
abline(lm1$coeff[1],lm1$coeff[2], col=2,lwd=3)

#Visualize the difference between genders
boxplot(edata[1,] ~ pdata$gender)
points(edata[1,] ~ jitter(as.numeric(pdata$gender)),col=as.numeric(pdata$gender))

#“dummy” variables for categorical variables
dummy_m = pdata$gender=="M"
dummy_m <- dummy_m*1

lm2 = lm(edata[1,] ~ pdata$gender)
tidy(lm2)
tidy(lm(edata[1,] ~ pdata$tissue.type ))
#Adjusting for variables,consider both age and gender
lm3 = lm(edata[1,] ~ pdata$age + pdata$gender)
tidy(lm3)
#You can add interaction terms, but there are a couple of things to keep in mind:
#Interpretation can be challenging
#You should always also include main effects
#Interactions in genomics data are in general hard to detect, and can be fraught for example in gene by environment interactions, the interaction term only tells you that there may be an association, but there are lots of reasons for this statistical association that aren’t biological interaction.
lm4 = lm(edata[1,] ~ pdata$age * pdata$gender)
tidy(lm4)

index = 1:19
lm5 = lm(edata[6,] ~ index)
plot(index,edata[6,],col=2)
abline(lm5,col=1,lwd=3)
lm6 = lm(edata[6,-19] ~ index[-19])
abline(lm6,col=3,lwd=3)
legend(3,1000,c("With outlier","Without outlier"),col=c(1,3),lwd=3)

#In general you’d like your residuals to looks symmetrically (e.g. approximately Normally) distributed but they aren’t here. Outliers in the residuals aren’t great either.
par(mfrow=c(1,2))
hist(lm6$residuals,col=2)
hist(lm5$residuals,col=3)
log_data <- log2(edata[1,]+1)
lm7 <-lm(log_data ~ index)
plot(index,log_data)
abline(lm7,col=3,lwd=5)
hist(lm7$residuals,col=1)
#Be careful when two variables in your regression model are very highly correlated. This is called co-linearity and can lead to highly variable and uninterpretable results. R fails “gracefully” (i.e. doesn’t tell you) when this happens so you have to check by hand. This is also a problem if you fit too many variables for example.
lm8 = lm(log_data ~ pdata$tissue.type + pdata$age)
tidy(lm8)
#Another good idea is to look for “patterns” in the residuals
colramp = colorRampPalette(trop)(17)
lm9 = lm(edata[2,] ~ pdata$age)
plot(lm9$residuals,col=colramp[as.numeric(pdata$tissue.type)])
library(limma)
library(edgeR)
library(edge)
library(devtools)
library(Biobase)
library(dplyr)
library(gplots)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))

#Here we will transform the data and remove lowly expressed genes.
edata = log2(edata + 1)
edata = edata[rowMeans(edata) > 10, ]
fdata = fData(bot)

#Fit many regression models at once.
mod = model.matrix(~ pdata$strain)
fit = lm.fit(mod,t(edata))#very fast for the multiple regression model fitting
names(fit)
fit$coefficients[,1]#coefficient from the first model fit and equals to the first gene LM model
lm <- lm(edata[1,] ~ pdata$strain)
tidy(lm)
plot(edata[1,], col=as.numeric(pdata$strain))
plot(pdata$strain,edata[1,])
plot(edata[1,] ~ pdata$strain,border=c(1,2))
points(edata[1,] ~ jitter(as.numeric(pdata$strain)),col=as.numeric(pdata$strain))
abline(lm,col=3,lwd=3)
#Look at the coefficients across genes
par(mfrow=c(1,2))
hist(fit$coefficients[1,],breaks=100,col=2,xlab="Intercept")
hist(fit$coefficients[2,],breaks=100,col=2,xlab="Strain")
abline(v=0,lwd=3,col=1)
#Look at the residuals for a couple of genes
par(mfrow=c(1,2))
plot(fit$residuals[,1],col=2)
plot(fit$residuals[,2],col=2)
#Fit many regressions with an adjustment
mod_adj = model.matrix(~ pdata$strain + as.factor(pdata$lane.number))
fit_adj = lm.fit(mod_adj,t(edata))
fit_adj$coefficients[,1]
#Fit many regressions with the limma package
fit_limma = lmFit(edata,mod_adj)
names(fit_limma)
fit_limma$coefficients[1,]
#Fit many regressions with the edge package
edge_study = build_study(data=edata,grp=pdata$strain,adj.var=as.factor(pdata$lane.number))
fit_edge = fit_models(edge_study)
summary(fit_edge)
fit_edge@beta.coef[1,]
#install.packages("~/Desktop/jackstraw_1.3.tar.gz", repos = NULL, type = "source")


library(sva)
library(bladderbatch)
library(snpStats)
data("bladderdata")
pheno = pData(bladderEset)
edata = exprs(bladderEset)
#Adjusting for batch effects with a linear model
#We will use two models. One with the variable we care about (cancer status) and the other that is just the known adjustment variables (in this case we will assume none)
mod = model.matrix(~pheno$cancer+as.factor(pheno$batch))
fit <- lm.fit(mod,t(edata))#coefficients for the variable across all the samples
hist(fit$coefficients[2,],breaks=100,col=2)

table(pheno$cancer,pheno$batch)

#combat is another way to adjust for batch effects
#Adjusting for batch effects with Combat
#Another approach is to use Combat. Combat returns a “cleaned” data matrix after batch effects have been removed. Here we pass a model matrix with any known adjustment variables and a second parameter that is the batch variable.

batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)#only intercept
modcancer = model.matrix(~cancer, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)#clean the dataset with combat
combat_fit = lm.fit(modcancer,t(combat_edata))
hist(combat_fit$coefficients[2,],col=2,breaks=100)
#Comparing Combat and linear adjustment
#We can compare the estimated coefficients from Combat and linear adjustment by looking at the right coefficients for each model.

plot(fit$coefficients[2,],combat_fit$coefficients[2,],col=2,
     xlab="Linear Model",ylab="Combat",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)

#Adjusting for batch effects with sva
mod = model.matrix(~cancer,data=pheno)#with the variable I care about
mod0 = model.matrix(~1, data=pheno)#the null model
sva1 = sva(edata,mod,mod0,n.sv=2)#compare two models, sva1 is the new covariates created
#See if any of the variables correlate with batch
summary(lm(sva1$sv ~ pheno$batch))#correlate the batch effects with the actual observed batch effect
# the second surragate variable is highly correlated
boxplot(sva1$sv[,2] ~ pheno$batch)
points(sva1$sv[,2] ~ jitter(as.numeric(pheno$batch)),col=as.numeric(pheno$batch))
#Add the surrogate variables to the model matrix and perform the model fit
modsv = cbind(mod,sva1$sv)
fitsv = lm.fit(modsv,t(edata))
#Compare the fit from surrogate variable analysis to the other two.

par(mfrow=c(1,2))
plot(fitsv$coefficients[2,],combat_fit$coefficients[2,],col=2,
     xlab="SVA",ylab="Combat",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)
plot(fitsv$coefficients[2,], fit$coefficients[2,],col=2,
     xlab="SVA",ylab="linear model",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)


#Principal components for populatin structure
#Load an example data set and take a smaller subset of samples for computational efficiency

data(for.exercise)
controls <- rownames(subject.support)[subject.support$cc==0]
use <- seq(1, ncol(snps.10), 10)
ctl.10 <- snps.10[controls,use]
#Calculate the PCs

xxmat <- xxt(ctl.10, correct.for.missing=FALSE)
evv <- eigen(xxmat, symmetric=TRUE)
pcs <- evv$vectors[,1:5]
#Let’s compare the PCs to the population labels and see that PC1 captures the population variable very well

pop <- subject.support[controls,"stratum"]
plot(pcs[,1],pcs[,2],col=as.numeric(pop),
     xlab="PC1",ylab="PC2")
legend(0,0.15,legend=levels(pop),pch=19,col=1:2)

