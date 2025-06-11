#read
USgdist = readRDS(file="~Data/USgdist.rds")
#write.csv(USgdist, file="~Data/USgdist.cvs",  row.names = FALSE)
#USgdist = read.csv("~Data/USgdist.cvs")

#deal with inf
AKidx = which(colnames(USgdist) == "AK");AKidx
DCidx = which(colnames(USgdist) == "WA");DCidx
HIidx = which(colnames(USgdist) == "HI");HIidx
CAidx = which(colnames(USgdist) == "CA");CAidx
USgdist[2,] = USgdist[9,] + 1
USgdist[,2] = USgdist[,9] + 1
USgdist[12,] = USgdist[5,] + 1
USgdist[,12] = USgdist[,5] + 1
USgdist[2,2] = 0
USgdist[12,12] = 0

save(USgdist, file = "/Users/jane/GRA/Data/USgdist_beginning.RData")

USgdist_exp = exp(-USgdist)
lowerI = lower.tri(USgdist_exp)
USgdist_expv = c(USgdist_exp[lowerI])

#Exp Matrix - with different wright try 1-10
distmateix = function(weit){
  USgdist_exp = exp(-USgdist/weit)
  lowerI = lower.tri(USgdist_exp)
  USgdist_expv = c(USgdist_exp[lowerI])
  return(USgdist_expv)
}

weitU1=0.8
weitL1=1
weits1=seq(weitU1,weitL1,length=5)
weitU2=1.5
weitL2=4
weits2=seq(weitU2,weitL2,length=6)
weits = c(weits1,weits2,6,8,10)

USgdist_mtxs = lapply(weits, distmateix)

for(w in weits){
  USgdist_nb = USgdist
  USgdist_nb[USgdist_nb != 1] = exp(-USgdist_nb[USgdist_nb != 1]/w)
  lowerI = lower.tri(USgdist_nb)
  USgdist_neb = c(USgdist_nb[lowerI])
  tbl_name = paste0("USgdist_neb", w)
  assign(tbl_name, USgdist_neb)
}

#Matrix - All 1
USgdist_constant = matrix(1, nrow = 51, ncol = 51)
lowerI = lower.tri(USgdist_constant)
USgdist_constantv = c(USgdist_constant[lowerI])

#Matrix only 0-1
USgdist_bi = USgdist
USgdist_bi[USgdist_bi != 1] = 0
lowerI = lower.tri(USgdist_bi)
USgdist_biv = c(USgdist_bi[lowerI])

for (w in weits) {
  USgdist_mtxs[[length(USgdist_mtxs) + 1]] = get(paste0("USgdist_neb", w))
}
USgdist_mtxs[[length(USgdist_mtxs) + 1]] = USgdist_constantv #29
USgdist_mtxs[[length(USgdist_mtxs) + 1]] = USgdist_biv #30

save(USgdist_mtxs, file = "/Users/jane/GRA/Data/USgdist.RData")

