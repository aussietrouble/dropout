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
  x = x/sqrt(n-1)  #why? so that the noise scales with n?
  beta0 = runif(p, min=5, max=20)
  z = 2*rbinom(p, 1, prob=rep(.5, p))-1 #randomly generates sign
  beta0 = beta0*z
  if (s < p) {
    beta0[(s+1):p] = 0
  }
  y = x %*% beta0 + sigma*rnorm(n)
  return(list(x=x, y=y, Sigma=Sigma, beta0=beta0))
}

risk = function(a,b,C) {
  return(t(a-b) %*% C %*% (a-b))
}


relative.risk = function(a,b,C) {
  return(risk(a,b,C)/risk(a,0,C))
}


get.marginals = function(X, Y){
  p = dim(X)[2]
  out = rep(0, p)
  for (i in 1:p){
    out[i] = lm(Y~0+X[,i])$coefficients
  }
  return(out)
}

trim = function(v, minval, maxval){
  out = pmin(v, maxval)
  out = pmax(out, minval)
  return(out)
}

grad.drop = function(n, p, s=p, sigma=1, rho=0, theta=1, T=1, adapt=F, init=F, seed=NA, silent=F) {
  if (!is.na(seed)) {
    set.seed(seed)
  }
  dat = gen.banded.data(n, p, rho, s, sigma)
  X = dat$x
  n = nrow(X)
  p = ncol(X)
  Y = dat$y
  Y.hat = 0*Y
  beta.hat = rep(0,p)
  thetat = theta * rep(1,p)
  if (init){
    thetat = abs(get.marginals(X, Y))
    thetat = p*theta * thetat/sum(thetat)
    thetat = trim(thetat, 10/n, 1)
  }
  nnz = 0
  
  dat$beta.ridge = solve(t(X)%*%X + diag((1-thetat)/thetat)) %*% t(X) %*% Y
  dat$risk.ridge = relative.risk(dat$beta0, dat$beta.ridge, dat$Sigma)
  
  for (t in 1:T) {
    R = Y - Y.hat
    #thetat = theta * rep(1,p)
    Z = mbinom(n, thetat)
    X.sparse = X * Z
    beta = solve(t(X.sparse)%*%X.sparse)%*%t(X.sparse)%*%R
    Y.hat = Y.hat + X.sparse %*% beta
    beta.hat = beta.hat + thetat * beta
    nnz = nnz + n*sum(thetat)
    if(!silent){
      cat(sprintf("\nt=%d, mean(|grad|)=%.2f, nnz=%.2f\n", t, mean(abs(thetat*beta)), n*sum(thetat)))
      print(round(rbind(dat$beta0, as.vector(beta.hat), as.vector(thetat*beta), as.vector(thetat)), digits=3))
    }
    # set next theta values
    if (adapt){
      thetat = abs(thetat * beta)
      thetat = p*theta * thetat/sum(thetat)
      thetat = trim(thetat, 10/n, 1)
    }
    
  }
  dat$beta = as.vector(beta.hat)
  dat$risk = relative.risk(dat$beta0, dat$beta, dat$Sigma)
  dat$nnz = nnz/T
  return(dat)
}

# example call:
# out = grad.drop(n=1000,  p=10, s=3, sigma=.5, rho=.2, theta=.5, T=5, seed=11)
# s is the number of nonzero coefficients; rho is the correlation in the cov = [\rho^{|i-j|]
