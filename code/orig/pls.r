
library(mnormt)

# Generate a random covariance matrix by
# adding random edges.  Like an Erdos-Renyi graph
# but the diagonal is incremented in each step
# to keep the spectrum non-negative.

rand.Sigma = function(p, m) {
  S = diag(1,p)
  for (i in 1:m) {
    edge = sample(p,2,replace=F)
    v = runif(1,min=0,max=2)
    s = sign(rnorm(1))
    S[edge[1],edge[2]] = S[edge[1],edge[2]] + s*v
    S[edge[2],edge[1]] = S[edge[2],edge[1]] + s*v
    S[edge[1],edge[1]] = S[edge[1],edge[1]] + v
    S[edge[2],edge[2]] = S[edge[2],edge[2]] + v
  }
  nonzeros = length(which(abs(S) > 0))
  zeros = length(which(abs(S) == 0))
  sparsity = (nonzeros-p)/(p*(p-1))
  L = diag(1/sqrt(diag(S)))
  S = L %*% S %*% L
  cat(sprintf("Generated %d by %d covariance with %d edges, sparsity level %.2f%%\n", p, p, m, 100-100*sparsity))
  return(S)
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

rand.model = function(Sigma) {
  p = nrow(Sigma)
  beta = (1:p)^{-.7}
  return(beta)
}

unif.model = function(Sigma) {
  p = nrow(Sigma)
  beta = 2*runif(p)
  return(beta)
}

gen.data = function(n, p, m) {
  Sigma = rand.Sigma(p,m)
  x = rmnorm(n, mean=rep(0,p), varcov=Sigma)
  x = scale(x)
  beta0 = unif.model(Sigma)
  sigma = 1.0
  y = x %*% beta0 + sigma*rnorm(n)
  return(list(x=x, y=y, Sigma=Sigma, beta0=beta0))
}

gen.banded.data = function(n, p, rho, s=p) {
  Sigma = banded.Sigma(p, rho)
  x = rmnorm(n, mean=rep(0,p), varcov=Sigma)
  x = scale(x)
  x = x/sqrt(n-1)
  beta0 = 10*runif(p)
  if (s < p) {
    beta0[(s+1):p] = 0
  }
  sigma = 1.0
  y = x %*% beta0 + sigma*rnorm(n)
  return(list(x=x, y=y, Sigma=Sigma, beta0=beta0))
}

sparsify = function(X, alpha=1) {
  n = nrow(X)
  p = ncol(X)
  B = matrix(rbinom(n*p, size=1, prob=alpha), n, p)
  Xs = X * B
  Xs = scale(Xs, center=F)/sqrt(n-1)
  return(Xs)
}

pls = function(n, p, rho=0, s=p, sprob=.75) {
  dat = gen.banded.data(n, p, rho, s)
  X = dat$x
  Y = dat$y
  dat$betahat = solve(t(X)%*%X) %*% t(X) %*% Y 
  Yhat = 0*Y + mean(Y)
  Xs = sparsify(X,alpha=sprob)
  for (j in 1:p) {
    alpha = t(Xs) %*% Y
    cat(sprintf("|alpha|=%g\n", sqrt(sum(alpha*alpha))))
    Z = Xs %*% alpha
    Z = Z / sqrt(sum(Z*Z))
    Yhat = Yhat + sum(Z*Y) * Z
    Xs = sparsify(Xs - Z %*% t(Z) %*% Xs, alpha=sprob)
  }
  dat$yhat = Yhat
  return(dat)
}

gds = function(n, p, s=p, etac=100, T=1000) {
  dat = gen.banded.data(n, p, .7, s)
  X = dat$x
  n = nrow(X)
  Y = dat$y
  bhat = 0*t(t((X[1,])))
  eta = etac
  for (t in 1:T) {
    bhat = bhat + eta*t(X)%*%(Y - X%*%bhat)/n
#    eta = eta*sqrt(t)/sqrt(t+1)
    eta = etac/sqrt((t+1))
  }
  dat$bhat = bhat
  return(dat)
}

