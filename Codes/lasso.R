## Simulate data: n= 100 indivuduals, 1 QTL on chr 1 at 48.8 cM (marker 6)
setwd("../LCR.QTL_2014-01/Codes/")
source("functions.R")
data(map10)
plot(map10)
n <- 100
#a <- 2 * sqrt(0.08/ (1 - 2*0.08))
a <- 0.8
mymodel <- rbind(c(1, 48.84, a))
sim <- sim.cross(map10, type='bc', n.ind=n, model=mymodel)
sim$geno[[1]]$map
y <- sim$pheno[1:n,]
X <- sim$geno[[1]]$data
X <- factorX(X)
df <- data.frame(y, X)


## fit a simple linear model
out <- lm(y ~ ., data = df)
summary(out)
step(out, direction="both")
## Lasso
library(glmnet)
cv=cv.glmnet(X,y, alpha = 1)
model=glmnet(X,y,type.gaussian="covariance",lambda=cv$lambda.min, alpha=1)
predict(model, type="coefficients")
coef(model)
