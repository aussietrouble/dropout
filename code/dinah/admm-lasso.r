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
    #beta0[(s+1):p] = 0
    beta0[sample(p, p-s)] = 0
  }
  y = x %*% beta0 + sigma*rnorm(n)
  return(list(x=x, y=y, Sigma=Sigma, beta0=beta0))
}

soft.threshold = function(x,lambda) {
  return(sign(x) * (abs(x)-lambda) * (abs(x) > lambda))
}

admm.lasso = function(X, Y, lambda, beta0, alpha=1, alpha0=1, T=10, rho=.1, t, ose=F, silent=F) {
  require(Matrix)
  x=0*X[1,]; w=0*X[1,]; z=1*X[1,]
  A = X
  theta = rep(alpha0, ncol(X))
  
  S = 0
  Yproj = Y
  if(ose){
    S = cwsmall(nrow(X), t)
    Yproj = S %*% Y
  }
  
  for (k in 1:T) {
    A = X
    Z = mbinom(nrow(X), theta)
    Z = scale(Z, center=F, scale = theta)
    A = X * Z
    if(ose){
      A = as(A, "sparseMatrix")
      #print(dim(S))
      #print(dim(A))
      A = S %*% A
    }
    
    x = solve(t(A) %*% A + rho * Diagonal(ncol(X)) ) %*% (t(A) %*% Yproj + rho*z - w)
    z = soft.threshold(x + w/rho, lambda/rho)
    w = w + rho*(x-z)
    
    theta = rep(1, ncol(X))
    theta[which(z==0)] = alpha

    
    if(!silent){
      cat(sprintf("\niter=%d\n", k))
      dsp = cBind(beta0, z, x, w, theta)
      colnames(dsp) = c("beta", "z", "x", "w", "theta")
      print(dsp) 
    }
    
  }
  return(list(beta = beta0, z =z))
}

risk = function(a,b,C) {
  return(t(a-b) %*% C %*% (a-b))
}


relative.risk = function(a,b,C) {
  return(risk(a,b,C)/risk(a,0,C))
}

lasso = function(n, p, s=p, lambda, sigma=1, alpha=1, alpha0=1, t=t, rho=0, T=10, seed=NA, ose=F, silent=F) {
  if (!is.na(seed)) {
    set.seed(seed)
  }
  dat = gen.banded.data(n, p, rho, s, sigma)
  X = dat$x
  n = nrow(X)
  p = ncol(X)
  Y = dat$y
  start = proc.time()
  out = admm.lasso(X, Y, lambda, beta0=dat$beta0, alpha=alpha, alpha0=alpha0, t=t, ose=ose, silent=silent)
  end = proc.time()
  out$risk = relative.risk(dat$beta0, out$z, dat$Sigma)
  out$relerr = relative.risk(dat$beta0, out$z, diag(p))
  out$time = end-start
  return(out)
}

# example call:
# out = lasso(n=10000, p=10, s=3, lambda=.7, sigma=1, alpha=.2, seed=10)
# s is the number of nonzero coefficients; lambda is the lasso regularization;
# alpha is the sampling rate theta_j for the current "irrelevant variables" in admm
# the relevant variables are given theta_j=1


