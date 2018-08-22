rm(list = ls())
source("aeS.R")
source("sgTs_thre.R")
source("thresh.R")
source("autocovm.R")
source("opbw.R")
source("wntest.R")
library(foreach)
library(expm)
library(MASS)
library(fastclime)
library(Rmpi)
library(snow)
library(doSNOW)
library(foreach)
library(doRNG)
library(MASS)
library(stats)
library(parallel)
library(doParallel)
cl = getMPIcluster()
# cl = makeCluster(2, type = "SOCK")
ncores = length(cl)
registerDoParallel(cl)
clusterEvalQ(cl, c(library(MASS),
                   library(foreach),
                   library(fastclime),
                   library(expm),
                   library(clime)))
ncores = length(cl)
print(paste0("ncores_", ncores))

NAI = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
m = NAI %% 3
if(ceiling(NAI/3) %in% c(1, 2)){
  p = 15
}else{
  p = 50
}
# p = 50
n = 300

if(m == 0){
  S1 = diag(1, p, p)
  for(ii in c(1:p)){
    for(jj in c(1:p)){
      S1[ii, jj] = 0.995^(abs(ii-jj))
    }
  }
  S11 = sqrtm(S1)
}else if(m == 1){
  S <- diag(1, p, p)
  ll <- ceiling(p/2.5)
  for(k in c(1:floor(p/ll))){
    S[(ll*(k-1)+1):(ll*k),(ll*(k-1)+1):(ll*k)] = 0.8
  }
  S1 <- S + diag(0.2, p, p)
  S11 <- sqrtm(S1)
}else if(m == 2){
  S11 <- matrix(runif(n = p^2, min = -1, max = 1), nrow = p)
  A <- S11
  S1 <- A %*% t(A)
}
iid <- 1 - ceiling(NAI/3) %in%c(1, 3)
# lambda_list <- c(seq(0.01, 0.4, by = 0.03))
lambda_list <- c(seq(0.001, 0.01, by = 0.002))
# for(jj in c(1:20)){
  # message(jj)
  # result_name <- paste0("res_p_", p, "_iid_", iid, "_", m,"_jj_", jj, "_clime.rda")
  res_full <- foreach(ll = lambda_list)%do%{
    res_50 <- foreach(ii = c(1:100), .combine = "rbind")%dopar%{
      message(ii)
      print(ii)
      set.seed(ii+200)
      if(ceiling(NAI/3) %in%c(1, 3)){
        r0 <- matrix(runif(n = p, min = 1/4, max = 1/2), nrow = p)
        r1 <- rep(0, p)
        for(ii in c(1:p)){
          r1[ii] <- runif(n = 1, min = 0, max = 1/2 - r0[ii])
        }
        # r1 <- t(r1)
        Z1 <- rnorm(n = p) * sqrt(r0/(1 - r1))
        Z2 <- rnorm(n = p) * sqrt(r0/(1 - r1))
        Z <- cbind(Z1, Z2)
        e <- matrix(rnorm(n = p * (n - 2)), ncol = n-2)
        for(jj in c(3 : n)){
          Z_pre <- sqrt(r0 + r1 * Z[, jj - 1]^2) * e[, jj-2]
          Z <- cbind(Z, Z_pre)
        }

        X <- S11 %*% Z

      }else{
        X = S11 %*% matrix(rt(n*p, df = 8), ncol = n)
      }

      k_max = 10
      kk = seq(2, k_max, 2)
      M = 2000
      k0 = 10
      delta = 1.5
      alpha = 0.05
      wntest(X, M, k_max, kk, type = 1, opt = 1, lambda = ll)
    }
    list(ll, res_50)
  }
  # save(res_full, file = result_name)
  # rm(res_full)
# res_full <- foreach(ll = lambda_list)%do%{
#   res_50 <- foreach(ii = c(1:100), .combine = "rbind")%dopar%{
#     message(ii)
#     print(ii)
#     set.seed(ii+200)
#     if(ceiling(NAI/3) %in%c(1, 3)){
#       r0 <- matrix(runif(n = p, min = 1/4, max = 1/2), nrow = p)
#       r1 <- rep(0, p)
#       for(ii in c(1:p)){
#         r1[ii] <- runif(n = 1, min = 0, max = 1/2 - r0[ii])
#       }
#       # r1 <- t(r1)
#       Z1 <- rnorm(n = p) * sqrt(r0/(1 - r1))
#       Z2 <- rnorm(n = p) * sqrt(r0/(1 - r1))
#       Z <- cbind(Z1, Z2)
#       e <- matrix(rnorm(n = p * (n - 2)), ncol = n-2)
#       for(jj in c(3 : n)){
#         Z_pre <- sqrt(r0 + r1 * Z[, jj - 1]^2) * e[, jj-2]
#         Z <- cbind(Z, Z_pre)
#       }
#
#       X <- S11 %*% Z
#
#     }else{
#       X = S11 %*% matrix(rt(n*p, df = 8), ncol = n)
#     }
#
#     k_max = 10
#     kk = seq(2, k_max, 2)
#     M = 2000
#     k0 = 10
#     delta = 1.5
#     alpha = 0.05
#     wntest(X, M, k_max, kk, type = 1, opt = 1, lambda = ll)
#   }
#   list(ll, res_50)
# }


# apply(res_15, 2, mean)
iid <- 1 - ceiling(NAI/3) %in%c(1, 3)
result_name <- paste0("res_p_", p, "_iid_", iid, "_", m, "_clime.rda")
.Last = NULL
save.image(result_name)
stopCluster(cl)
mpi.quit()

q(runLast = F)

