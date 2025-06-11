helm = function(n) {
  mat = matrix(0, n - 1, n) 
  i = 2:n
  r = 1 / sqrt( i * (i - 1) )
  for ( j in 1:(n - 1 ) )   mat[j, 1: c(j + 1) ] = c( rep(r[j], j), - j * r[j]) 
  mat
}

helm_transoform = function(x){
  dx = ncol(x) #p
  H = helm(dx)
  M1 = orthogonal_comp(H)
  x = x %*% M1
  return(x)
}

inverse_helm_transform = function(beta_tilde) {
  dx = ncol(beta_tilde) + 1
  H = helm(dx)
  M1 = t(orthogonal_comp(H))  # Compute the same transformation matrix
  M1_inv = ginv(M1)        # Compute the pseudo-inverse of M1
  beta_original = M1_inv %*% t(beta_tilde)  # Transform back
  return(t(beta_original))
}