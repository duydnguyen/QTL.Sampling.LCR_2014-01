rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#setwd("../Google Drive/projects/LCR.QTL_2014-01")
library(qtl)
library(lhs)
library('DoE.base')
source('fncs.R')
## Simulate data
length.chr <-500
num.chr <- 4
num.mk <- 500
map <- sim.map(rep(length.chr, num.chr), num.mk, include.x=FALSE)
n.ind <- 2000
a <- 0.7
mymodel <- rbind(c(1, 50, a), c(2, 150, a))
simA <- sim.cross(map, type="bc", n.ind=n.ind, model=mymodel)

### Try real data
# data(listeria)
# simA <-listeria
# simA <- sim.geno(simA, n.draws=128, step=1)
## Sampling
s <- 5
n <- s^2 #run size 
chr1 <- 4
chr2 <- 4
## create a U design
design <- genU(nfactors=2, nlevels=s)
## Sample markers from U design
sam <- LHD2(size=n,design=design,chr=c(chr1,chr2),length.chr=c(length.chr,length.chr))
## Compute lod.full (from evalLODf) 
lod.full <- c(rep(0,dim(sam)[1]))
lod.fv1 <- c(rep(0,dim(sam)[1]))
for (i in 1:dim(sam)[1])
  if ((!is.na(sam$chr1[i]) & (!is.na(sam$chr2[i])))){
    lod.full[i] <- evalLODf(chr1=chr1, chr2=chr2,marker1=as.character(sam$markers1[[i]]), 
                            marker2=as.character(sam$markers2[[i]]))$lod.full
    lod.fv1[i] <- evalLODfv1(chr1=chr1, chr2=chr2,marker1=as.character(sam$markers1[[i]]), 
                             marker2=as.character(sam$markers2[[i]]),
                             sample=sam, fulldata=FALSE)$lod.fv1
  }
sam <- data.frame(sam, lod.full=lod.full, lod.fv1=lod.fv1) 
#sam[which(sam$lod.full>0.3),]
## plot
dev.off()
cutoff <- quantile(sam$lod.full, probs=0.9)
index <-which(sam$lod.full>cutoff)
sam1 <- sam
sam1$lod.full <- c(rep(0, length(sam$lod.full)))
sam1$lod.full[index] <-sam$lod.full[index]
plotQTL2.full(chr1=chr1, chr2=chr2, s=s, sam1)

#evalLODf(chr1,chr2, marker1='D1M473', marker2='D2M485')
#plot(out2.sim.mr)


sam.chr3 <- sam$lod.full
sam.chr4 <- sam$lod.full
lod <-t(rbind(t(sam.chr3), t(sam.chr4)))
lod  <- as.vector(lod)
df <- data.frame(lod.full=lod, 
                 groups = factor(c(rep(3,25),rep(4,25))) )
bwplot(lod~groups, data=df, ylab="lod.full", xlab="Chromosomes", 
       main="LOD_f of chromosomes 3 and 4"  )



       
