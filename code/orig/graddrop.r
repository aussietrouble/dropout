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

risk = function(a,b,C) {
  return(t(a-b) %*% C %*% (a-b))
}


relative.risk = function(a,b,C) {
  return(risk(a,b,C)/risk(a,0,C))
}

update.theta = function(alpha, thetat, beta, adaptive=FALSE) {
  theta = rep(alpha, length(thetat))
  return(theta)
}

grad.drop = function(n, p, s=p, sigma=1, rho=0, alpha=1, T=1, seed=NA) {
  if (!is.na(seed)) {
    set.seed(seed)
  }
  dat = gen.banded.data(n, p, rho, s, sigma)
  X = dat$x
  n = nrow(X)
  p = ncol(X)
  Y = dat$y
  Y.hat = 0*Y
  thetat = rep(alpha, p)
  beta = rep(1, p)
  beta.hat = rep(0,p)
  for (t in 1:T) {
    R = Y - Y.hat
    thetat = update.theta(alpha, thetat, beta)
    Z = mbinom(n, thetat)
    X.sparse = X * Z
    beta = solve(t(X.sparse)%*%X.sparse)%*%t(X.sparse)%*%R
    Y.hat = Y.hat + X.sparse %*% beta
    beta.hat = beta.hat + thetat * beta
    cat(sprintf("\nt=%d, mean(|grad|)=%.2f, nnz=%.2f\n", t, mean(abs(thetat*beta)), n*sum(thetat)))
    print(rbind(dat$beta0, as.vector(beta.hat), as.vector(thetat*beta), as.vector(thetat)))
  }
  dat$beta = as.vector(beta.hat)
  dat$risk = relative.risk(dat$beta0, dat$beta, dat$Sigma)
  return(dat)
}

# example call:
# out = grad.drop(n=1000,  p=10, s=3, sigma=.5, rho=.2, alpha=.5, T=5, seed=11)
# s is the number of nonzero coefficients; rho is the correlation in the cov = [\rho^{|i-j|]
