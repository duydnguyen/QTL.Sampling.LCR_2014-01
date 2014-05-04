rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#setwd("../Google Drive/projects/LCR.QTL_2014-01")
library(lhs)
source('fncs.R')
## Simulate data
length.chr <-500
num.chr <- 2
num.mk <-500
map <- sim.map(rep(length.chr, num.chr), num.mk, include.x=FALSE)
n.ind <- 5000
a <- 0.7
mymodel <- rbind(c(1, 50, a), c(2, 150, a))
simA <- sim.cross(map, type="bc", n.ind=n.ind, model=mymodel)
#out2.sim.mr <- scantwo(simA, method='mr',chr=c(1,2),verbose=FALSE)
#summary(out2.sim.mr)
#### LHD
# sampling size
n <-20
## create a LHD design: n= sample size, k = number of chrs
design <- floor(maximinLHS(n=n , k=2)*n)
## Sample markers from LHD design
sam <- LHD2(size=n,design=design,chr=c(1,2),length.chr=c(length.chr,length.chr))
## Compute lod.full (from evalLODf) 
lod.full <- c(rep(0,dim(sam)[1]))
lod.fv1 <- c(rep(0,dim(sam)[1]))
for (i in 1:dim(sam)[1])
  if ((!is.na(sam$chr1[i]) & (!is.na(sam$chr2[i])))){
    lod.full[i] <- evalLODf(chr1=1, chr2=2,marker1=as.character(sam$markers1[[i]]), 
           marker2=as.character(sam$markers2[[i]]))$lod.full
    lod.fv1[i] <- evalLODfv1(chr1=1, chr2=2,marker1=as.character(sam$markers1[[i]]), 
                             marker2=as.character(sam$markers2[[i]]),
                             sample=sam, fulldata=FALSE)$lod.fv1
  }
sam <- data.frame(sam, lod.full=lod.full, lod.fv1=lod.fv1)  
## plot
dev.off()
plotQTL2.full()





