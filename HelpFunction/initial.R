initial = function(x, y, z, n, lam, varth, gam, lambda_star){
  dx = ncol(x) #p
  dz = ncol(z) #q
  
  Pz = solve(t(z) %*% z) %*% t(z)
  Qz = diag(n) - z%*%Pz #diag() construct a diagonal matrix. f not sure about the dim of I
  D_ini = matrix(0,nrow = n,ncol = (n*(n-1)/2))
  e = diag(n)
  p = dx
  I_p = diag(p) #PXP matrix
  x_block = t(x[1,]) 
  for (i in 1:(n-1)){
    ind = (i-1)*n-i*(i-1)/2 
    D_ini[,(ind+1):(ind+n-i)]=e[,i]-e[,(i+1):n] #D is the difference between ei and ej, is Delta in old code
    x_block = bdiag(x_block, t(x[(i+1),]))
  }
  A = kronecker(t(D_ini), I_p)
  AA =  t(A) %*% A
  
  #Find Beta0 and eta0 
  beta_Rp1 = (t(x_block) %*% Qz %*% x_block + lambda_star * AA) #ATA = np x np, so XTQzX need to be np x np, then need XD
  beta_Rp1 = solve(beta_Rp1)
  beta_R = beta_Rp1 %*% t(x_block) %*% Qz %*% y
  eta_R = solve(t(z) %*% z) %*% t(z) %*% (y - x_block %*% beta_R)
  
  beta_R = matrix(beta_R, nrow = dx, ncol = n)
  

  rank_idx = data.frame(rank(apply(beta_R, 2, median))) #calculate median by row and give rank number

  colnames(rank_idx) = "rankidx"
  rank_idx$origional_idx = row.names(rank_idx) #add the original order
  ranked_idx = rank_idx[order(rank_idx$rankidx),]
  K_star = floor(0.5*n^(1/2)) #number of group, we have n given
  ranked_idx$groupassign_idx = base::cut(ranked_idx$rankidx, breaks = K_star, labels = FALSE)
  
  groupassign_mtx = matrix(0, nrow = n, ncol = K_star) #w
  for (k in 1:K_star) {
    group_idx = as.numeric(ranked_idx$origional_idx[ranked_idx$groupassign_idx == k])
    thecol = groupassign_mtx[,k]
    thecol[group_idx] = 1
    groupassign_mtx[,k] = thecol
  }
  
  beta0_group = matrix(0, nrow = K_star, ncol = dx)
  eta0_group = matrix(0, nrow = K_star, ncol = dz)
  for (k in 1:K_star) {
    #k=1
    group_idx = as.numeric(ranked_idx$origional_idx[ranked_idx$groupassign_idx == k])
    
    #current group data
    x_group = x[group_idx, ]
    y_group = y[group_idx]
    z_group = z[group_idx, ]
    
    #least square
    beta0_group[k,] = solve(t(x_group) %*% x_group, tol = 1e-50) %*% t(x_group) %*% y_group
    eta0_group[k,] = solve(t(z_group) %*% z_group, tol = 1e-50) %*% t(z_group) %*% (y_group - x_group %*% beta_R[,k])
  }
  
  # TODO: Assign beta and eta to observation according to assigned group
  beta0 = groupassign_mtx %*% beta0_group #dx x n
  eta0 = groupassign_mtx %*% eta0_group
  
  #delta initiate2
  delta_ini = t(D_ini) %*% beta0
  
  #initiate ups
  ups_ini = matrix(0, nrow = (n*(n-1)/2), ncol = dx)
  return(list(x_block = x_block, 
              D_ini = D_ini, 
              A = A, 
              AA = AA, 
              beta0 = beta0, 
              eta0 = eta0, 
              delta_ini = delta_ini, 
              ups_ini = ups_ini, 
              Qz = Qz, 
              K_star = K_star))
}
