
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

gen.banded.data = function(n, p, rho) {
  Sigma = banded.Sigma(p, rho)
  x = rmnorm(n, mean=rep(0,p), varcov=Sigma)
  x = scale(x)
  beta0 = unif.model(Sigma)
  sigma = 1.0
  y = x %*% beta0 + sigma*rnorm(n)
  return(list(x=x, y=y, Sigma=Sigma, beta0=beta0))
}

risk.path = function(x, y, Sigma, beta0, lambda, thresh) {
  n = nrow(x)
  p = ncol(x)
  Shat = (1/n) * t(x) %*% x
  inds = which(abs(Shat) < thresh)
  active = which(abs(Shat) >= thresh)
  Shat[inds] = 0
  edges = (length(active)-p)/2
  cat(sprintf("threshold=%f, edges=%d\n", thresh, edges))
  z = (1/n) * t(x) %*% y
  risk = rep(0, length(lambda))
  for (i in 1:length(lambda)) {
    A = solve(Shat + lambda[i]*diag(p))
    betahat = A %*% z
    risk[i] = t(beta0 - betahat) %*% Sigma %*% (beta0 - betahat)
    #cat(sprintf("i=%d, lambda=%f, risk=%f\n", i, lambda[i], risk[i]))
  }
  return(list(risk=risk, edges=edges))
}


test.sparsify = function() {
  p = 10
  m = 50
  Sigma = rand.Sigma(p,m)
  beta0 = unif.model(Sigma)
  sigma = 1.0
  firstN = 20
  lastN = 1000
  byn=10
  ns = seq(firstN, lastN, by=byn)
  trials = 500
  alpha = c(1.0, 0.9, 0.8, 0.7, 0.6, 0.5)
  meanrisk = matrix(0, length(ns), length(alpha))
  for (i in 1:length(ns)) {
    n = ns[i]
    cat(sprintf("n=%d\n", n))
    for (j in 1:length(alpha)) {
      risk = rep(0, trials)
      for (t in 1:trials) {
        x = rmnorm(n, mean=rep(0,p), varcov=Sigma)
        x = scale(x)
        y = x %*% beta0 + sigma*rnorm(n)
        w = matrix(rbinom(n*p, size=1, prob=alpha[j]), n, p)
        x = x * w
        Shat = (1/n) * t(x) %*% x
        z = (1/n) * t(x) %*% y
        A = solve(Shat)
        betahat = A %*% z
        risk[t] = t(beta0 - betahat) %*% Sigma %*% (beta0 - betahat)
      }
      meanrisk[i,j] = mean(risk)
      cat(sprintf("    alpha=%f, risk=%f\n", alpha[j], meanrisk[i,j]))
    }
  }
  matplot(ns, meanrisk, 'l', lwd=2, ylim=c(0,median(meanrisk[1,])))
}


test.banded = function(n, p, rho) {
  cat(sprintf("test.banded: n=%d, p=%d, rho=%f\n", n, p, rho))
  out = gen.banded.data(n, p, rho)
  x = out$x; y=out$y; beta0=out$beta0; Sigma=out$Sigma
  p = ncol(x)
  lambda = seq(0.001, .75, .01)
  edges = p*(p-1)/2
  dev.set(3)
  plot(0:1, rep(edges, 2), 'l', ylab='edges (computation)', xlab='', ylim=c(0, edges))
  out = risk.path(x, y, Sigma, beta0, lambda, 0)
  out.max = risk.path(x, y, Sigma, beta0, lambda, rho/2)
  dev.set(2)
  plot(lambda, out$risk, 'l', col=1, xlab='lambda', ylab='risk', ylim=c(0,min(out.max$risk)))
  thresh = seq(.01, rho/2, .002)
  for (i in length(thresh):1) {
    out = risk.path(x, y, Sigma, beta0, lambda, thresh[i])
    dev.set(2)
    lines(lambda, out$risk, 'l', col=i)
    dev.set(3)
    lines(0:1, rep(out$edges, 2), 'l', col=i)
  }
}

