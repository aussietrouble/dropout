#an attempt to translate Saurabh Paul's implementation of the CW sparse embeddings into R

cw = function(n, t){
  require(Matrix)
  k = t/4
  v=t/16
  
  q=t/k
  w=k/v
  
  #hash h: [n] -> [q]
  idx = sample(q, size=n, replace=TRUE)
  a = rep(0, q)
  #he initializes S here
  S = 0
  
  for (i in 1:q) {
    a[i] = length(which(idx==i))
  }
  
  #there are q of the B(i) matrices each of size k x a[i]
  
  for (j in 1:q){
    aj = a[j]
    #he initializes B here
    B = 0
    
    for (j2 in 1:w){
      #form phi*D_i
      idx2 = sample(v, size=aj, replace=TRUE)
      
      g = runif(aj)
      ent = rep(1, aj)
      ent[which(g <.5)] = -1
      Di = Diagonal(aj, ent)
      
      phi = sparseMatrix(i = idx2, j = 1:aj, x= 1, dims=c(v, aj))
      
      if(j2==1){
        B = phi %*% Di
      } else {
        B = rBind(B, phi %*% Di)
      }
      
    }
    
    B = B*sqrt(v/k)
    if(j==1){
      S = B
    } else{
      S = bdiag(S, B)
    }

  }
  
  p = sample(n)
  P <- as(p, "pMatrix")
  
  R = S %*% P
  return(R)
}


cwsmall = function (n, t){
  #form phi*D_i
  idx = sample(t, size=n, replace=TRUE)
  
  g = runif(n)
  ent = rep(1, n)
  ent[which(g <.5)] = -1
  D = Diagonal(n, ent)
  
  phi = sparseMatrix(i = idx, j = 1:n, x= 1, dims=c(t, n))
  
  S = phi %*% D
  return(S)    
}