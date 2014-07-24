theta = seq(.01,.99,.01)
risk = rep(0, length(theta))
n=10000
p=10
s=3
for (i in 1:length(theta)) {
  out = grad.drop(n,  p, s, sigma=.5, rho=.2, theta=theta[i], T=5, seed=11)
  risk[i] = out$risk
}




