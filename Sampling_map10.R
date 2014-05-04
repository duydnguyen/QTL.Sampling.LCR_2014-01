rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
dev.off()
#setwd("../Google Drive/projects/LCR.QTL_2014-01")
library(qtl)
library(lhs)
source('fncs.R')
## Simulate data: 20 chromosomes from map10
#data(map10)
# map: 20 chrs, each length 250cM and 1000 markers/chr
length.chr <-250
num.chr <- 20
num.mk <- 1000
map <- sim.map(rep(length.chr, num.chr), num.mk)
#map <- map10
n.ind <- 1000
# QTL effect
a <- 0.7
mymodel <- rbind(c(1, 50, a), c(1, 100, a-0.3), c(14, 100, a), c(8, 200, a))
simA <- sim.cross(map, type="bc", n.ind=n.ind, model=mymodel)
## MR method
out.mr <- scanone(simA, method="mr")
plot(out.mr, chr=c(1:19))
## Sampling: n = samples per chr
n <-50
# sample.mk <- LHD(size=n, chr=1, length.chr = length.chr)
# sample.mk
#QTL <- findQTL(sample.mk, quantile=.90)
QTL <- matrix(c(rep(0,num.chr * n)), nrow = num.chr, ncol = n)
for (i in 1:num.chr){
  sample.mk <- LHD(size=n, chr=i, length.chr = length.chr)
  res <- findQTL(sample.mk, quantile=.90, chr=i)
  QTL[i,] <- c(res$cM, c(rep(0,n-length(res$cM))))
}

# plot: panel of 19 chrs
source('fncs.R')
layout(matrix(c(1:20),nrow=4,ncol=5))
for (i in 1:19) plotChr(out.mr=out.mr, chr=i, n=n)

layout(matrix(c(1:3),nrow=1,ncol=3))
for (i in c(1,8,14)) plotChr(out.mr=out.mr, chr=i, n=n)
