source("aeS.R")
source("sgTs_thre.R")
source("thresh.R")
source("autocovm.R")
source("opbw.R")

library(foreach)
library(expm)
library(MASS)
library(fastclime)

p = 50
n = 300
S1 = diag(1, p, p)
for(ii in c(1:p)){
  for(jj in c(1:p)){
    S1[ii, jj] = 0.995^(abs(ii-jj))
  }
}
S11 = sqrtm(S1)
res_50 <- foreach(ii = c(1:100), .combine = "rbind", .errorhandling = "remove")%dopar%{
  message(ii)
  set.seed(ii+200)
  X = S11 %*% matrix(rt(n*p, df = 8), ncol = n)
  k_max = 10
  kk = seq(2, k_max, 2)
  M = 500
  k0 = 10
  delta = 1.5
  alpha = 0.05
  wntest(X, M, k_max, kk, type = 1, opt = 1)
}
apply(res_50, 2, mean)

