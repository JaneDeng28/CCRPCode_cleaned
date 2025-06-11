library(dplyr)
library(Matrix)
library(MASS)
library(fossil)
library(parallel)

stateinfo = readRDS("~Data/state_fips_name.rds")
load("~Data/USgdist.RData")

#Partition setting--------
design = 4
c1 <- c("CA", "WA", "OR", "NV", "AZ", "UT", "ID", "MT", "WY", "CO", "TX", "OK", "KS", "NM", "AK", "HI", "NE")
c2 <- c("ND", "SD", "MN", "IA", "MO", "AR", "LA", "MS", "AL", "TN", "KY", "IL", "WI", "MI", "IN", "GA", "FL")
c3 <- c("OH", "SC", "NC", "WV", "VA", "MD", "DC", "PA", "DE", "NJ", "NY", "CT", "RI", "MA", "NH", "ME", "VT")
c1idx = which(stateinfo$abbr %in% c1)
c2idx = which(stateinfo$abbr %in% c2)
c3idx = which(stateinfo$abbr %in% c3)

#Coefficient setting--------
p = 3
q = 2

n = nrow(stateinfo)
weightc1 = matrix(rep(c(2, 2), each = length(c1idx)), ncol = 2) + as.vector(rnorm(length(c1idx), 0, 0.1))
weightc2 = matrix(rep(c(1, 1), each = length(c2idx)), ncol = 2) + as.vector(rnorm(length(c2idx), 0, 0.1))
weightc3 = matrix(rep(c(-3, -3), each = length(c3idx)), ncol = 2) + as.vector(rnorm(length(c3idx), 0, 0.1))
eta = c(1,1)

#Simulate data-----

seedset = sample(100:1000, 50)
sim100 = lapply(seedset, function(seed) SimulationData(p,q,weightc1,weightc2,eta,c1idx,c2idx,seed))

#Simulation setup------

simfct = function(sim){

  z = sim$z; sumz = sim$sumz; y = sim$y; x = sim$x2; sumx = sim$sumx2; e = sim$e
  n = nrow(x); a = 10^(-3); lam = 0.2; varth=1; gam=3; lambda_star = 0.001
  
   #Initiate
  initial_result = initial(x=z, y, z=x, n, lam, varth, gam, lambda_star)
  x_block = initial_result$x_block; D_ini = initial_result$D_ini; A = initial_result$A
  AA = initial_result$AA; beta0 = initial_result$beta0; eta0 = initial_result$eta0
  delta_ini = initial_result$delta_ini; ups_ini = initial_result$ups_ini
  
  lamU=30; lamL=0.1
  
  #Tuning Lambda
  estimation_result = 
    estimate_process(lamL,lamU,a, x=z, y, z=x, n, lam, varth, gam, x_block, D_ini, 
                  A, AA, beta0, eta0, delta_ini, ups_ini, Qz, USgdist)
  lambda_star = estimation_result$lambda_star; Group = estimation_result$Group
  beta_est = estimation_result$beta; eta_est = estimation_result$eta
  r_est = estimation_result$final_residual;rr = estimation_result$residual_trace
  
  return(list(lambda_star = lambda_star, 
              Group = Group, 
              beta = beta_est, 
              eta = eta_est, 
              final_residual = r_est, 
              residual_trace = rr))
}


#Estimation Parallelized------

#Choose spatial weight
caselist = c(6,7,20,21,29)

for (i in caselist){
  USgdist = USgdist_mtxs[[i]]
  tbl_name = paste0("simresult", i)
  start_time = Sys.time()
  simresult = mclapply(sim100, simfct, mc.cores = 6);end_time=Sys.time();time_taken=end_time - start_time; print(time_taken)
  assign(tbl_name, simresult)
  save(simresult, file = paste0("~/US",design,"simresult", i, ".RData"))
  }










