estimate_process = function(lamL,lamU,a, x, y, z, n, lam, varth, gam, x_block, D_ini, 
                         A, AA, beta0, eta0, delta_ini, ups_ini, Qz, USgdist){
  betain=beta0
  lams1 = seq(lamL,4,length=10) #you can have a finer tuning on smaller lambda value
  lams2 = seq(4,lamU,length=5)
  lams = c(lams1,lams2)
  
  allBIC = lapply(lams, 
                  function(lam) ModifiedBIC(a, x, y, z, n, lam, varth, gam, x_block, D_ini, 
                                        A, AA, beta0, eta0, delta_ini, ups_ini, Qz, USgdist))

  allBIC = unlist(allBIC)
  BIC = allBIC[-1]
  lambda_star = lams[(max(which(BIC==min(BIC)))+1)]
  
  result = estimation(a, x, y, z, n, lam=lambda_star, 
                      varth,  gam, x_block, D_ini, A, AA, 
                      beta0, eta0, delta_ini, ups_ini,Qz, 
                      USgdist)
  
  Group = result$group
  beta_est = result$beta_result
  eta_est = result$eta
  r_est = result$final_residual
  rr = result$residual_trace
  step = result$n_iter
  x_block = result$x_block
  W = result$W
  alphaold = result$alpha
  
  return(list(lambda_star = lambda_star, 
              Group = Group, 
              beta = beta_est, 
              eta = eta_est, 
              final_residual = r_est, 
              residual_trace = rr, 
              n_iter = step, 
              x_block = x_block, 
              alpha = alphaold,
              W = W))
}

