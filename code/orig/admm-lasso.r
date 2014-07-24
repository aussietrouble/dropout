library(mnormt)

mbinom = function(n, theta) {
  p = length(theta)
  B = t(matrix(rbinom(n*p, size=1, prob=theta), p, n))
  return(B)
}

banded.Sigma = function(p, rho) {
  S = diag(1,p)
  for (j in 1:p) {
    for (k in 1:p) {
      S[j,k] = rho^abs(j-k)
    }
  }
  cat(sprintf("Generated %d by %d covariance [%f^|i-j|]\n", p, p, rho))
  return(S)
}


gen.banded.data = function(n, p, rho, s=p, sigma=1) {
  Sigma = banded.Sigma(p, rho)
  x = rmnorm(n, mean=rep(0,p), varcov=Sigma)
  x = scale(x)
  x = x/sqrt(n-1)
  beta0 = runif(p, min=5, max=20)
  z = 2*rbinom(p, 1, prob=rep(.5, p))-1
  beta0 = beta0*z
  if (s < p) {
    beta0[(s+1):p] = 0
  }
  y = x %*% beta0 + sigma*rnorm(n)
  return(list(x=x, y=y, Sigma=Sigma, beta0=beta0))
}

soft.threshold = function(x,lambda) {
  return(sign(x) * (abs(x)-lambda) * (abs(x) > lambda))
}

admm.lasso = function(X, Y, lambda, beta0, alpha1=1, alpha2=1,T=10, rho=.1) {
  x=0*X[1,]; y=0*X[1,]; z=1*X[1,]
  A = X
  theta = rep(alpha2, ncol(X))
  Z = mbinom(nrow(X), theta)
  A = X * Z
  for (k in 1:T) {
    x = solve(t(A) %*% A + rho * diag(ncol(X))) %*% (t(A) %*% Y + rho*z - y)
    z = soft.threshold(x + y/rho, lambda/rho)
    y = y + rho*(x-z)
    A = X
    theta = rep(alpha1, ncol(X))
    theta[which(z==0)] = alpha2
    Z = mbinom(nrow(X), theta)
    A = X * Z
    cat(sprintf("\niter=%d\n", k))
    dsp = cbind(beta0, z, x, y, theta)
    colnames(dsp) = c("beta", "z", "x", "y", "theta")
    print(dsp)
  }
}

lasso = function(n, p, s=p, lambda, sigma=1, alpha1=1, alpha2=1, rho=0, T=1, seed=NA) {
  if (!is.na(seed)) {
    set.seed(seed)
  }
  dat = gen.banded.data(n, p, rho, s, sigma)
  X = dat$x
  n = nrow(X)
  p = ncol(X)
  Y = dat$y
  out = admm.lasso(X, Y, lambda, beta0=dat$beta0, alpha1=alpha1, alpha2=alpha2)
}

# example call:
# out = lasso(n=10000, p=10, s=3, lambda=.7, sigma=1, alpha=.2, seed=10)
# s is the number of nonzero coefficients; lambda is the lasso regularization;
# alpha is the sampling rate theta_j for the current "irrelevant variables" in admm
# the relevant variables are given theta_j=1



