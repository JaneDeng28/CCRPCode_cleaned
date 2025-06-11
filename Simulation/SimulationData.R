SimulationData = function(p,q,weightc1,weightc2,weightc3,eta,c1idx,c2idx,c3idx,seed){
  n = nrow(stateinfo)
  set.seed(seed)
  
  #Compositional data
  X1 = matrix(0, nrow = n, ncol = p)
  for(i in 1:p){
    X1[,i] = runif(n)
  }
  # Normalize so they sum to 1
  X1 = data.frame(X1)
  X1[,p+1] = rowSums(X1[,1:p])
  for(i in 1:p){
    X1[,i] = X1[,i]/X1[,ncol(X1)]
  }
  X1 = log(X1)
  X1 = as.matrix(X1[,1:p])
  z = helm_transoform(X1)
  sumz = rep(0,n)
  sumz[c1idx] = rowSums(z[c1idx,]*(weightc1))
  sumz[c2idx] = rowSums(z[c2idx,]*(weightc2))
  sumz[c3idx] = rowSums(z[c3idx,]*(weightc3))
  
  #Non-compositional data
  x2 = matrix(0, nrow = n, ncol = q)
  for(i in 1:q){
    x2[,i] = runif(n)
  }
  x2 = data.matrix(x2)
  sumx2 = rep(0,n)
  sumx2 = rowSums(x2*eta)
  
  #Error term
  e = rnorm(n,0,0.1)
  
  #Outcome
  y = sumz + sumx2 + e
  
  return(list(z = z, sumz = sumz, y = y, x2 = x2, sumx2 = sumx2, e = e))
}
