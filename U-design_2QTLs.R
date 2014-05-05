library('DoE.base')
#oa.design(nfactors=2, nlevels=3,nruns=9)
#### Orthogonal Array
s <- 3
OA <- matrix(c(1,1,1,2,2,2,3,3,3,
               1,2,3,1,2,3,1,2,3),
             nrow=nrow, ncol=ncol, byrow=FALSE)
genU(OA, s=s)


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


