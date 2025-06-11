estimation = function(a, x, y, z, n, lam, varth, gam, x_block, D_ini, 
                      A, AA, beta0, eta0, delta_ini, ups_ini, Qz, USgdist){
  dx = ncol(x) #p
  dz = ncol(z) #q
  Pz = solve(t(z) %*% z) %*% t(z)
  Qz = diag(n) - z%*%Pz 
  
  D_weight = lam * USgdist
  
  betaold = beta0 #beta0
  etaold = eta0 #eta0
  deltaold = t(delta_ini) #dx x n(n-1)/2
  upsold = t(ups_ini)
  a = a #ep
  step = 0
  r_new = 2
  rr = numeric(300)
  rr[1] = 10
  
  
  beta_newp1 = t(x_block) %*% Qz %*% x_block + varth * AA
  beta_newp1 = solve(beta_newp1)
  
  while(r_new > a & step < 300){ 
    step = step+1

    beta_newp3 = deltaold - (1/varth) * upsold 
    beta_newp3 = beta_newp3 %*% t(D_ini)
    vbeta_newp3 = as.vector(beta_newp3)
    beta_newp2 = t(x_block) %*% Qz %*% y + varth * vbeta_newp3
    beta_new = beta_newp1 %*% beta_newp2
    beta_newmtx = matrix(beta_new,nrow=dx,ncol=n)
    
    #update eta
    eta_new = Pz %*% (y - x_block %*% beta_new)
    
    #group-wise soft threshold operator
    BetaD = t(D_ini) %*% t(beta_newmtx)
    BetaD = t(BetaD)
    zeta = BetaD + (1/varth) * upsold
    zeta = t(zeta)
    zetanorm = apply(zeta^2, 1, sum)
    zetanorm = sqrt(zetanorm)
    lam1 = lam/varth
    deltaS = 1 - lam1/zetanorm
    deltaS = (((deltaS > 0) * deltaS) %*% t(seq(1, 1,length = dx))) * zeta
    
    #update delta
    delta_new = (deltaS/(1-(varth*gam)^(-1))) * (zetanorm <= (gam*D_weight)) + 
      zeta * (zetanorm > (gam*D_weight))
    delta_new = t(delta_new)
    
    #update ups (v)
    ups_new = upsold + varth * (BetaD - delta_new)
    
    #update r (residual)
    r_new = BetaD - delta_new #r is the primal residual
    r_new = as.vector(r_new)
    r_new = sqrt(t(r_new) %*% r_new)
    rr[step] = r_new
    
    betaold = beta_new
    etaold = eta_new
    upsold = ups_new
    deltaold = delta_new
  }
  
  beta_result=matrix(betaold,nrow=dx,ncol=n)
  beta_result=t(beta_result)
  delta_result=t(deltaold)
  
  pair_index <- combn(n, 2)
  delta_result_t=apply(abs(delta_result),1,sum)
  group = rep(0, n)
  used = rep(FALSE, n)
  K = 1
  
  #Clustering
  for (i in 1:n) {
    if (used[i]) next
    group[i] = K
    for (k in 1:ncol(pair_index)) {
      ii = pair_index[1, k]
      jj = pair_index[2, k]
      if (ii == i && delta_result_t[k] == 0 && !used[jj]) {
        group[jj] = K
        used[jj] = TRUE
      }
      if (jj == i && delta_result_t[k] == 0 && !used[ii]) {
        group[ii] = K
        used[ii] = TRUE
      }
    }
    used[i] = TRUE
    K = K + 1
  }
  
  K = max(group)
  Wtilda = matrix(0, nrow=n, ncol=K)
  for (k in 1:K) {
    Wtilda[group == k, k] = 1
  }
  
  alphaold = matrix(0, nrow=K, ncol=dx)
  for (k in 1:K) {
    members = which(group == k)
    if (length(members) > 1) {
      alphaold[k, ] = colMeans(beta_result[members, , drop=FALSE])
    } else {
      alphaold[k, ] = beta_result[members, ]
    }
  }
  
  W = kronecker(Wtilda, diag(dx))
  alphaold = t(alphaold)
  alphaold = as.vector(alphaold)
  
  sig=sqrt(t(y-z%*%etaold-x_block%*%W%*%alphaold)%*%
             (y-z%*%etaold-x_block%*%W%*%alphaold)/(n-dz-K*dx))
  
  return(list(
    K = K,
    beta = betaold,
    beta_result = beta_result,
    eta = etaold,
    ups = upsold,
    delta = deltaold,
    x_block = x_block,
    group = group,
    alpha = alphaold,
    Wtilda = Wtilda, # 10
    W = W,
    sigma = sig,
    final_residual = r_new,
    residual_trace = rr,
    n_iter = step
  ))
}
