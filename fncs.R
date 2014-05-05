## Simulate data 
sim <- function(len, n.mar, n.ind){
  mapA <<- sim.map(len = len, n.mar= n.mar, include.x= FALSE, eq.spacing=TRUE )
  plot(mapA)
  simA <<- sim.cross(mapA, type='bc', n.ind=n.ind, model=mymodel)
}

## Compute LOD for each marker: simA =simulated data, chr = chromosome,
#marker = 'D1M256', n=individuals
evalLOD <- function(simA, chr, marker, n){
  index <- which(names(simA$geno[[chr]]$data[1,])==marker)
  x <- factor(simA$geno[[chr]]$data[,index])
  y <- simA$pheno[,1]
  F <- summary(aov(y ~ x))[[1]]$F[1]
  return((n/2)*log10(F/(n-2) +1))
}
#which(names(simA$geno[[1]]$data[1,])=='D1M256')
## Compute lod.full: 
evalLODf <- function(chr1,chr2, marker1, marker2){
  index1 <- which(names(simA$geno[[chr1]]$data[1,])==marker1)
  index2 <- which(names(simA$geno[[chr2]]$data[1,])==marker2)
  q1 <- factor(simA$geno[[chr1]]$data[,index1])
  q2 <- factor(simA$geno[[chr2]]$data[,index2])
  y <- simA$pheno[,1]
  l.0 <- logLik(lm(y~1))[1]
  l.f <- logLik(lm(y~q1+q2))[1]
  return(data.frame(lod.full = log10(exp(l.f - l.0))))
}
## Compute lod.1
evalLOD1 <- function(chr, marker){
  index1 <- which(names(simA$geno[[chr]]$data[1,])==marker)
  q1 <- factor(simA$geno[[chr]]$data[,index1])
  y <- simA$pheno[,1]
  l.0 <- logLik(lm(y~1))[1]
  l.1 <- logLik(lm(y~q1))[1]
  return(data.frame(lod.1 = log10(exp(l.1 - l.0))))
}
## Compute lod.fv1(s,t): markers(s,t)
evalLODfv1 <- function(chr1,chr2, marker1, marker2, sample, fulldata){
  log.full <- evalLODf(chr1,chr2, marker1, marker2)[1,]
  M.1 <- evalM1(chr1, chr2,sample,fulldata)
  if ((log.full - M.1)>0)
    lod.fv1 <- log.full - M.1
  else lod.fv1 = 0
  return(data.frame(lod.fv1 = lod.fv1))
}
## Compute M.f(j,k): chrs=j,k
evalMf <- function(chr1, chr2){
  marker1 <- names(simA$geno[[chr1]]$map)
  marker2 <- names(simA$geno[[chr2]]$map)
  max <- -Inf
  for (i in 1:length(marker1))
    for (j in 1:length(marker2)){
      log.full <- evalLODf(chr1=chr1,chr2=chr2, marker1=marker1[i], 
                           marker2=marker2[j])[1,]  
      if (log.full > max) max <- log.full
    }
  return(max)
}
## Compute M.1(j,k): chrs=j,k
evalM1 <- function(chr1, chr2, sample, fulldata){
  if (fulldata){
    marker1 <- names(simA$geno[[chr1]]$map)
    marker2 <- names(simA$geno[[chr2]]$map)
  } else{
    marker1 <- as.vector(sample$markers1)
    marker2 <- as.vector(sample$markers2)
  } 
  max1 <- -Inf
  for (i in 1:length(marker1)){
    lod1 <- evalLOD1(chr=chr1,marker=marker1[i])[1,]
    if (lod1 > max1) max1 <- lod1
  }
  max2 <- -Inf
  for (i in 1:length(marker2)){
    lod1 <- evalLOD1(chr=chr2,marker=marker2[i])[1,]
    if (lod1 > max2) max2 <- lod1
  }
  return(max(max1,max2))  
}
## Compute M.fv1(j,k): chrs=j,k
evalMfv1 <- function(chr1,chr2,sample,fulldata){
  return(evalMf(chr1,chr2) - evalM1(chr1,chr2,sample,fulldata))
}
## Sample from LHD: size = sample size, order(design)
LHD <- function(size, chr, length.chr){
  design <- sort(maximinLHS(n=size , k=1) * length.chr)
  foo <- simA$geno[[chr]]$map
  sample.mk <- c(rep(0,size))
  names <- c(rep('',size))
  for (i in 1:(size)){
    if (i==1) temp <-foo[ (foo <= design[1]) & (foo >= 0) ]
      else temp <-foo[ (foo <= design[i]) & (foo >= design[i-1]) ]
    names.temp <- names(temp)
    temp <- as.vector(temp)
    index <- sample(x=c(1:length(temp)), size=1)
    if (length(names.temp[index]) > 0){
      sample.mk[i] <- temp[index]
      names[i]<- names.temp[index]
    } else{
      sample.mk[i] <- temp[1]
      names[i] <- names.temp[1]
    }  
  }
  return(data.frame(markers = names,cM = sample.mk))
}
## Find the 95th quantile of LODs from sampled markers
findQTL <-function(sample.mk, quantile, chr){
  # sample size
  n <- dim(sample.mk)[1]
  LOD <- c(rep(0,n))
  for (i in 1:n){
    marker <- sample.mk$markers[i]
    if (!is.na(marker))
      LOD[i] <- evalLOD(simA, chr=chr, marker=marker, n =n.ind)
  }
  index <- which(LOD>quantile(LOD, probs=c(quantile)))
  return(data.frame(sample.mk[index,], LOD = LOD[index]))
}
## Plot true QTL and sampled QTL given a chr
# out.mr = mr model; chr = chr to plot; n= samples per chr
plotChr <- function(out.mr, chr, n){
  pos.init <- (chr-1)*num.mk + 1
  pos.end <- pos.init + num.mk -1
  chr.pos <- out.mr$pos[pos.init:pos.end]
  chr.lod <- out.mr$lod[pos.init:pos.end]
  xlab <- paste('Chr',chr, sep=' ')
  plot(chr.pos, chr.lod, type='l', ylab="LOD", xlab=xlab )
  len <- n - length(which(QTL[chr,]==0))
  points(x =QTL[chr,1:len], y= c(rep(0,len)), pch = 8, col='red')
  index <- which(mymodel[,1]==chr)
  points(x =mymodel[index,2], y= c(rep(0.5,length(mymodel[index, 2]))), 
         pch = 6, col='blue')
}
## 2-dim LHD: size= sample size; design = LHD design on total chrs; 
#  chr = vector of chrs to be sampled from; length.chr = vector of chr lengths
LHD2 <- function(size, design, chr, length.chr){
  pos.chr1 <- LHD(size=size, chr=chr[1],length.chr=length.chr[1])$cM
  mk.chr1 <- LHD(size=size, chr=chr[1],length.chr=length.chr[1])$markers
  pos.chr2 <- LHD(size=size, chr=chr[2],length.chr=length.chr[2])$cM
  mk.chr2 <- LHD(size=size, chr=chr[2],length.chr=length.chr[2])$markers
  # fix case when having unbalanced samples
  len1 <- length(pos.chr1[design[,1]])
  len2 <- length(pos.chr2[design[,2]])
  len = min(len1,len2)
  return(data.frame(markers1=mk.chr1[1:len],chr1=pos.chr1[design[,1]][1:len], 
                    markers2=mk.chr2[1:len],chr2=pos.chr2[design[,2]][1:len] ))
}
## plot() for 2 qtls
plotQTL2.full <-function(){
  chr1<-1
  chr2<-2
  # cutoff <- quantile(sam$lod.full,probs=0.8)
  # index <- which(sam$lod.full>cutoff)
  plot(x = sam$chr1, y = sam$chr2, pch = 1, col='black', 
       xlab='chromosome', ylab='chromosome',  
       xlim=c(0,length.chr), ylim=c(0,length.chr), main='Sampling Scheme')
  # plot lod.full
  symbols(x = sam$chr1, y = sam$chr2, circles=sam$lod.full,
          inches=1/3, ann=F, bg="steelblue2", fg=NULL,
          xlab='chromosome', ylab='chromosome',
          xlim=c(0,length.chr), ylim=c(0,length.chr), add=TRUE)
  # plot lod.fv1
  symbols(x = sam$chr1, y = sam$chr2, circles=sam$lod.fv1,
          inches=1/3, ann=F, bg="red", fg=NULL, add=TRUE)
  # points(x = sam$chr1[index], y = sam$chr2[index], pch = 8, col='red', 
  #      xlab='chromosome', ylab='chromosome')
  grid(5,5, lwd =1, col="black")
  # plot true QTLs
  index.x <- which(mymodel[,1]==chr1)
  index.y <- which(mymodel[,1]==chr2)
  abline(v=mymodel[index.x[1],2], col='blue')
  abline(h=mymodel[index.y[1],2], col='blue')
}
## Plot for comparing U design and LHS
plotU_LHD <- function(U){
  #x1=c(1,2,3,6,5,4,8,9,7)
  #x2=c(3,4,8,1,5,7,2,6,9)
  x1 <- U[,1]
  x2 <- U[,2]
  lhd <- ceiling(maximinLHS(n=9 , k=2)*9)
  plot(x1,x2, xlab='chromosome', ylab='chromosome', main='U design vs. LHD')
  points(lhd, pch=3, col = 'blue')
  grid(3,3, lwd=2, col = 'black')
  legend('topleft', c('U design', 'LHD') , pch=c(1,3),
         col=c('black','blue'), bty='n', cex=.75)
}
## generate U design from an Orthogonal Array OA of s symbols
genU <- function(OA, s){
  nrow <- s^2 #also run size
  ncol <- 2
  ## randomized OA
  OA <- OA[sample(c(1:nrow), size = nrow),]
  OA <- OA[,sample(c(1:ncol), size = ncol)]
  rsymbols <- sample(c(1:s)) 
  OA.rand <- matrix(c(rep(0,nrow*ncol)),nrow=nrow,ncol=ncol)
  for (i in 1:ncol)
    for (j in 1:s ){
      index <- which(OA[,i]==j)
      OA.rand[index, i] <- rsymbols[j]
    }
  ## for each column, replace s positions with entry k by a random permutation 
  #  on (k-1)s+1,..., ks
  U <- matrix(c(rep(0,nrow*ncol)),nrow=nrow,ncol=ncol) 
  for (i in 1:ncol)
    for (k in 1:s){
      index <- which(OA.rand[,i]==k)
      val <- sample(c((k-1)*s+1):(k*s))
      U[index,i] <- val
    }
  return(U)
}