## R COMMANDS FOR STACKR 
## ANALYSIS OF VCF FILE

## AUTHOR : CHARLES PERRIER

#####################################
##	SET WORKING DIRECTORY

setwd("~/Documents")
#####################################




#####################################
##	INSTALL MISSING LIBRARIES

biocLite()
# press "a" and then "y" to the questions


library("BiocInstaller", lib.loc="/usr/local/lib/R/site-library")

biocLite("BiocInstaller")
# press "a" and then "y" to the questions


biocLite("gdsfmt")
#####################################




#####################################
##	LOAD LIBRARIES

library("stackr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(reshape)
library(gplots)
#####################################




#####################################
## DATA WORK

vcf.fn <-  "batch_1.vcf"

## REFORMAT DATA

snpgdsVCF2GDS(vcf.fn, "data.gds", method="biallelic.only")
snpgdsSummary("data.gds")

## RELOAD FORMATED DATA

genofile.ad <- snpgdsOpen("data.gds")
samp.id=read.gdsn(index.gdsn(genofile.ad, "sample.id"))
#####################################




#####################################
##	IBS 1

IBSad<-snpgdsIBS(genofile.ad, sample.id=NULL, snp.id=NULL, verbose=TRUE, autosome.only =FALSE)

#IBSad

## WRITE DATA ON MATRIX FORM, RE-LOAD IT AND SAVE THE FILE

write.table(IBSad$ibs, "IBSad.txt", col.names = F, row.names = F)
IBSadhm= read.table("IBSad.txt", header= F, dec=".", sep=" ")
IBSadhm=data.matrix(IBSadhm)
diag(IBSadhm)=NA
write.table(IBSadhm, "IBSad2.txt")

## MAKE HEAT MAP

#heatmap.2(as.matrix(IBSadhm), dendrogram = "none", scale="none", key=T, keysize=1,density.info="none",trace="none",cexCol=0.5,Colv=FALSE,Rowv=FALSE, cexRow = 0.5, col="rainbow")

heatmap.2(as.matrix(IBSadhm), dendrogram = "both", scale="none", key=T, keysize=1,density.info="none",trace="none",cexCol=0.5,Colv=FALSE,Rowv=FALSE, cexRow = 0.5, col="rainbow")
#####################################




#####################################
##	IBS 2
ibs <- snpgdsIBS(genofile.ad, num.thread=2,autosome.only =FALSE)
image(ibs$ibs, col=terrain.colors(16))

loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[, 1]
y <- loc[, 2]

colvec<-c(rep("forestgreen",10), rep("black",10), rep("orange",11))

plot(x, y, xlab = "", ylab = "",col=colvec,main = "Multidimensional Scaling Analysis (IBS)")
#####################################






#####################################
##	IBS 3
set.seed(100)

ibs <- snpgdsHCluster(snpgdsIBS(genofile.ad, num.thread=2,autosome.only =FALSE))

rv <- snpgdsCutTree(ibs)

snpgdsDrawTree(rv, type="z-score", main="bt")

snpgdsDrawTree(rv, main="bt", edgePar=list(col=rgb(0.5,0.5,0.5, 0.75), t.col="black"))

plot(rv$dendrogram, leaflab="none"
)
table(rv$samp.group)


##	ADD PATH TO POP_MAP
pop_cod=read.table("")

pop_code=pop_cod$V1

rv2 <- snpgdsCutTree(ibs, samp.group=as.factor(pop_code))

plot(rv2$dendrogram, leaflab="none")
#####################################





#####################################
##	PCA

pca <- snpgdsPCA(genofile.ad, missing.rate=0.25,maf=0.001,snp.id = NULL, autosome.only = FALSE, num.thread=2)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = pca$sample.id, EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2],stringsAsFactors = FALSE)

colvec<-c(rep("forestgreen",10), rep("black",10), rep("orange",11))

plot(tab$EV1, tab$EV2, col=colvec, xlab="EV1 4.98%", ylab="EV2 0.55%", cex=0.3)

#plot(tab$EV2, tab$EV4, col=colvec, xlab="EV3 4.98%", ylab="EV4 0.55%", cex=0.3)
#####################################


#####################################
##	CALCULATE FST FOR THE WHOLE DATASET

##	ADD PATH TO POP_MAP
popmap <- read.table("", header = FALSE)

popmapf <- as.factor(popmap$V2)

snpgdsFst(genofile.ad, popmapf, method="W&C84", autosome.only = FALSE )
#####################################
