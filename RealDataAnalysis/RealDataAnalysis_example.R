library(dplyr)
library(MASS)
library(fossil)
library(Matrix)
library(stats)
library(expm)
library(colorspace)
library(usmap)
library(tidyr)
library(purrr)

source("~/helm_transoform.R")
source("~/Function/initial.R")
source("~/Function/inverse_helm.R")
source("~/Real_data/estimate_full_realus.R")
stateinfo = readRDS("~/Data/state_fips_name.rds")
FIPScode = read.csv("~/Data/us-state-ansi-fips.csv")
FIPScode$stusps = substring(FIPScode$stusps, 2)

load(file = "~/Data/USgdist.RData")
load(file = "~/Data/cleans22_19.RData")

#y = beta*z + eta*x
#z: compositional data - income
#x: other

oriincome22 = as.data.frame(Incomenew22)
oriincome22$State = FIPScode$stusps

#z
z = data.matrix(Incomealr22)
Z = helm_transoform(z)

#x
x_mtx22 = left_join(GDP22clean %>% dplyr::select(GeoName, Manufacturing), 
                    Arthritis2022 %>% dplyr::rename(Arthritis = Data_Value) %>% 
                      dplyr::select('LocationDesc', 'Arthritis'), 
                    by = c(GeoName = 'LocationDesc')) 
x_mtx22 = left_join(x_mtx22, smoke2022 %>% dplyr::rename(Smoke = Data_Value) %>% 
                      dplyr::select('LocationDesc', 'Smoke'), 
                    by = c(GeoName = 'LocationDesc'))
x_mtx22 = left_join(x_mtx22, obese2022 %>% dplyr::rename(Obese = Data_Value) %>% 
                      dplyr::select('LocationDesc', 'Obese'), 
                    by = c(GeoName = 'LocationDesc'))                     
x_mtx22 = left_join(FIPScode,x_mtx22, by=c(stname = 'GeoName'))
X = data.matrix(x_mtx22[,4:7])
X <- scale(X)

#y
y = log(COPD22$DataValue) 

#Estimation------

estimation_real = function(x = Z, y, z = X, USgdist){
n = nrow(x); a = 10^(-3); lam = 0.2; varth=1; gam=3; lambda_star = 0.001
#Initiate
initial_result = initial(x, y, z, n, lam, varth, gam, lambda_star)
x_block = initial_result$x_block; D_ini = initial_result$D_ini; A = initial_result$A
AA = initial_result$AA; beta0 = initial_result$beta0; eta0 = initial_result$eta0
delta_ini = initial_result$delta_ini; ups_ini = initial_result$ups_ini

#setting lambda range
lamU=1; lamL=0.05
  
#Tuning Lambda
estimation_result = 
  estimate_full(lamL,lamU,a, x, y, z, n, lam, varth, gam, x_block, D_ini, 
                A, AA, beta0, eta0, delta_ini, ups_ini, Qz, USgdist)
lambda_star = estimation_result$lambda_star; Group = estimation_result$Group
beta_est = estimation_result$beta; eta_est = estimation_result$eta
r_est = estimation_result$final_residual;rr = estimation_result$residual_trace
min_BIC = estimation_result$min_BIC

return(list(lambda_star = lambda_star, 
            Group = Group, 
            beta = beta_est, 
            eta = eta_est, 
            final_residual = r_est, 
            residual_trace = rr, 
            min_BIC = min_BIC))
}

#parallelized------

caselist = c(4,5,18,19,29)
for (i in caselist){
  USgdist = USgdist_mtxs[[i]]
  tbl_name = paste0("Result", i )
  start_time = Sys.time()
  result = estimation_real(x, y, z, USgdist);end_time=Sys.time();time_taken=end_time - start_time; print(time_taken)
  assign(tbl_name, result)
  save(result, file = paste0("~/RealUSResult", i, ".RData"))
}
