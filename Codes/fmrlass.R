#file <- '../LCR.QTL_2014-01/fmrlasso_1.0.tar.gz'
#install.packages(file, repos = NULL, type="source")
########## fmrlassopath and cvfmrlassopath ##############
sim <- function(x,prob,beta,ssd){
  ## Purpose: generates data from FMR model
  ## ----------------------------------------------------------------------
  ## Arguments: x (design-matrix);prob (component probabilities)
  ## beta ((p+1)*k Matix of regression coef);ssd (standard dev.)
  ## ----------------------------------------------------------------------
  k <- length(prob)
  n <- dim(x)[1]
  y <- numeric(n)
  s <- sample(1:k,n,replace=TRUE,p=prob)
  for (i in 1:n){
    y[i] <- rnorm(1,mean=beta[,s[i]]%*%x[i,],sd=ssd[s[i]])
  }
  y
}
set.seed(1)
#model specifications
n <- 150
p <- 50
x <- cbind(rep(1,n),matrix(rnorm(n*p),n,p))
beta <- cbind(c(0,3,3,0,0,0,rep(0,p-5)),c(0,-4,-2,0,0,0,rep(0,p-5)),
              c(0,-3,0,0,5,0,rep(0,p-5)))
prob <- c(0.3,0.3,0.4)
ssd <- c(0.5,0.5,0.5)
#generate data
y <- sim(x=x,ssd=ssd,beta=beta,prob=prob)
#initialisation of E-step
ex.ini <- ini.ex(k=3,n=n)
#set lambda grid
la <- seq(5,30,length=15)   
#fit fmrlassopath
fit <- fmrlassopath(x,y,k=3,lambda=la,ssd.ini=0.5,ex.ini=ex.ini)
View(fit$coef[1:51,1:3,7])
plotfmrlassopath(fit)
which.min(fit$bic)
#perform crossvalidation
fit.cv <- cvfmrlassopath(x,y,k=3,lambda=la,ssd.ini=0.5,ex.ini=ex.ini)
fit.cv$cv
which.min(fit.cv$cv)
