rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#setwd("../Google Drive/projects/LCR.QTL_2014-01")
library(qtl)
library(lhs)
source('fncs.R')
################################################
###                SINGLE QTL               ###
################################################
## Simulate data: even-spaced markers, all global vars
length.chr <- 5000
number.marker <- 2000
n.ind <- 1000
mymodel <<- rbind(c(1, 10, 0.7), c(1, 300, 0.7),c(1, 2000, 0.7),c(1, 3000, 0.7),c(1, 4000, 0.7))
sim(len =length.chr, n.mar = number.marker, n.ind)
#simA$geno[[1]]$data
#View(simA$geno[[1]]$map)
## MR method
out.mr <- scanone(simA, method="mr")
#plot(out.mr, chr=c(1))
#max(out.mr)
## Compute LOD for each marker
#evalLOD(simA=simA, chr=1, marker='D1M256', n=n.ind)
## LHD Sampling
# sample size
n <-20
sample.mk <- LHD(size=n, chr=1, length.chr = length.chr)
sample.mk
QTL <- findQTL(sample.mk, quantile=.80, chr=1)
# plot
plot(out.mr, chr=c(1))
points(mymodel[,2],c(rep(0,5)), pch = 6, col='blue', )
points(QTL$cM,QTL$LOD, pch = 8, col='red', )

plot(x=mymodel[,2], y=mymodel[,2], pch = 6, col='blue', xlab = 'Chromosome', ylab ='cM', 
     xlim=c(0,length.chr), ylim=c(0,length.chr) )
points(x = QTL$cM,y = QTL$cM, pch = 8, col='red')
