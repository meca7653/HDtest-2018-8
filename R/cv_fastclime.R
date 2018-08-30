#cross validation
cv_fastclime <- function(X, fold = 5, lambda = seq(1e-4, 1e-2, length.out = 50)){
  n = dim(X)[1]
  p = dim(X)[2]
  ii = ll = NULL #define the variables
  res_cv <- foreach(ii = c(1:fold), .combine = "rbind")%do%{
    ind_train <- c(1:n)[-c(c(1:floor(n/fold)) + (ii - 1) * (n/fold))]
    train <- X[ind_train,]
    L <- fastclime(train, lambda.min = 1e-7, nlambda = 3000)
    test <- X[c(c(1:floor(n/fold)) + (ii - 1) * (n/fold)),]
    res_pre <- foreach(ll = lambda, .combine = "cbind")%do%{
      out2 <- fastclime.selector(L$lambdamtx, L$icovlist,ll)
      M_inv <- out2$icov
      norm(M_inv %*% cov(test) - diag(1, p, p), type = "F")
    }
    res_pre
  }
  # la <- lambda[which.min(apply(res_cv, 2, mean))]
  la <- lambda[which(order(apply(res_cv, 2, mean)) == 30)]
  L_final <- fastclime(X, lambda.min = 1e-7, nlambda = 3000)
  out <- fastclime.selector(L_final$lambdamtx, L_final$icovlist,la)
  M_inv_final <- out$icov
  return(M_inv_final)
}






