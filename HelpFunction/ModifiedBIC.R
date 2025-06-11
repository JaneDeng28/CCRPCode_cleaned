ModifiedBIC = function(a, x, y, z, n, lam, varth, gam, x_block, D_ini, 
                   A, AA, beta0, eta0, delta_ini, ups_ini, Qz, USgdist){

  result=estimation(a, x, y, z, n, lam, varth, gam, x_block, D_ini, 
                    A, AA, beta0, eta0, delta_ini, ups_ini, Qz, USgdist)
  
  etahat=as.vector(as.matrix(result$eta))
  alphahat=result$alpha
  W=result$W
  Khat=result$K
  x_block=as.matrix(result$x_block)
  df=Khat*ncol(x)+ncol(z)
  Qn=t(y-z%*%etahat-x_block%*%W%*%alphahat)%*%(y-z%*%etahat-x_block%*%W%*%alphahat)/n
  Qn=Qn[1,1]
  
# Modified BIC
  BIC=log(Qn)+log(n*ncol(x)+ncol(z))*log(n)*df/(n)

  return(BIC)
}