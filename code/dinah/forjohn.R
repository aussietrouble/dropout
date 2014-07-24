

set.seed(64)
library(glmnet)
testdat = gen.banded.data(n=20000, p=250, rho=.7, s=20, sigma=1) #cov is rho^|i-j|
testfit = cv.glmnet(testdat$x, testdat$y, intercept=F)
lmin = testfit$lambda.min

#note that I'm not doing any clarkson-woodruff stuff here
out2=lasso(n=20000, p=250, s=20, rho=.7, lambda=lmin, sigma=1, alpha=1, alpha0=1, t=4096, ose=F, silent=T, seed=64)

#these two should be the same, but they're not :(
relative.risk(testdat$beta0, coef(testfit, s=lmin)[2:251], testdat$Sigma)
out2$risk


# This is what I get on my computer:
# > relative.risk(testdat$beta0, coef(testfit, s=lmin)[2:251], testdat$Sigma)
# [,1]
# [1,] 0.01897039
# > out2$risk
# 1 x 1 Matrix of class "dgeMatrix"
# [,1]
# [1,] 0.05856597