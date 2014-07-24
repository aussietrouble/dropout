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

dropout = function(X, Y, theta0) {
  A = X
  theta = rep(theta0, ncol(X))
  Z = mbinom(nrow(X), theta)
  A = X * Z
  betahat = solve(t(A) %*% A) %*% t(A) %*% Y
  return(betahat)
}

subsample = function(X, Y, m) {
  n = nrow(X)
  inds = sample(1:n, m)
  A = X[inds, ]
  W = Y[inds, ]
  betahat = solve(t(A) %*% A) %*% t(A) %*% W
  return(betahat)
}

doss = function(n, p, thetas, sigma=1, rho=0, seed=NA) {
  if (!is.na(seed)) {
    set.seed(seed)
  }
  dat = gen.banded.data(n, p, rho, p, sigma)
  beta0 = dat$beta0
  do.error = rep(0, length(thetas))
  do.comp = rep(0, length(thetas))
  for (i in 1:length(thetas)) {
    betahat = dropout(dat$x, dat$y, theta0=thetas[i])
    error = sqrt(sum((betahat-beta0)^2))
    cat(sprintf("theta=%.3f  nnz/p=%.0f  error=%.3f\n", thetas[i], n*thetas[i], error))
    do.error[i] = error
  }

  ms = round(thetas*n/p)
  ss.error = rep(0, length(thetas))
  for (i in 1:length(ms)) {
    betahat = subsample(dat$x, dat$y, ms[i])
    error = sqrt(sum((betahat-beta0)^2))
    cat(sprintf("m=%d  m*p=%.0f  error=%.3f\n", ms[i], ms[i]*p, error))
    ss.error[i] = error
  }
  do.comp = thetas * n * p
  ss.comp = ms * (n*(p + p*(p-1)/2))
  plot(do.comp, do.error, 'l', lwd=2)
  lines(ss.comp, ss.error, 'l', lwd=2)
  return(data.frame(do.comp, do.error, ss.comp, ss.error))
}

plot.mse = function(n, p, thetas, do, ss, pdf=F) {
  do.mean = apply(do, 1, mean)
  do.error = apply(do, 1, sd)
  do.comp = thetas * n * p
  ss.mean = apply(ss, 1, mean)
  ss.error = apply(ss, 1, sd)
  ms = round(thetas*n/p)
  ss.comp = ms * p^2
  if (pdf) {
    pdf(sprintf("n%dp%d.pdf", n, p), width=10, height=8)
  }
  titlez = sprintf("n=%d   p=%d", n, p)
  plot(ss.comp, ss.mean, 'l', lwd=2, ylab='MSE', col='lightsalmon2', xlab='computation',ylim=c(min(do.mean),ss.mean[2]),main=titlez)
  lines(do.comp, do.mean, col='lightblue2', 'l', lwd=2)
  epsilon = (do.comp[2]-do.comp[1])*.1
  segments(do.comp, do.mean-do.error, do.comp, do.mean+do.error)
  segments(do.comp-epsilon, do.mean-do.error, do.comp+epsilon, do.mean-do.error)
  segments(do.comp-epsilon, do.mean+do.error, do.comp+epsilon, do.mean+do.error)
  segments(ss.comp, ss.mean-ss.error, ss.comp, ss.mean+ss.error)
  segments(ss.comp-epsilon, ss.mean-ss.error, ss.comp+epsilon, ss.mean-ss.error)
  segments(ss.comp-epsilon, ss.mean+ss.error, ss.comp+epsilon, ss.mean+ss.error)
  if (pdf) {
    dev.off()
  }
}

doit = function(n, p, T) {
  do = c()
  ss = c()
  thetas = seq(.02, .5, .02)
  for (t in 1:T) {
    cat(sprintf("\n\nSample %d of %d\n========================\n\n", t, T))
    dat = doss(n=n, p=p, thetas=thetas, rho=.2)
    do = cbind(do, dat$do.error)
    ss = cbind(ss, dat$ss.error)
    plot.mse(n, p, thetas, do, ss)
  }
  plot.mse(n, p, thetas, do, ss, pdf=T)
}

  




