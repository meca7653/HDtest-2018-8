sgTs_thre = function(X, k0, delta, opt=1, lambda = 0.1,
                     lambda_search = seq(1e-4, 1e-2, length.out = 50),
                     fold = 5, cv_opt = cv_opt){
  n = dim(X)[1]
  p = dim(X)[2]
  if (opt == 1){
    L <- fastclime(X, lambda.min = 1e-7, nlambda = 800)
    out2 <- fastclime.selector(L$lambdamtx, L$icovlist,lambda)
    M_inv <- out2$icov
    rr <- eigen(M_inv)
    vv <- sqrt(rr$value)
    vv[rr$value <= 0] = 0
    G <- rr$vector
    M1 <- (G) %*% diag(vv) %*% t(G)
    X1 = M1 %*% t(X)
    Xn = X1
  }else if(opt== 2){
    M = cov(X)
    M1 = sqrtm(solve(M))
    X1 = M1 %*% t(X) #A\B ==> solve(A)%*%B
    Xn = X1
  }else if(opt == 3){
    L <- clime(X, standardize = F, lambda.min = 1e-6)
    re.cv <- cv.clime(L, loss = "tracel2",fold = fold)
    re.clime.opt <- clime(X, standardize=FALSE, re.cv$lambdaopt)
    print(re.cv$lambdaopt)
    message(print(re.cv$lambdaopt))
    M_inv <- re.clime.opt$Omegalist[[1]]
    rr <- eigen(M_inv)
    vv <- sqrt(rr$value)
    vv[rr$value <= 0] = 0
    G <- rr$vector
    M1 <- (G) %*% diag(vv) %*% t(G)
    X1 = M1 %*% t(X)
    Xn = X1
  }else if(opt == 4){
    M_inv <- cv_fastclime(X, fold = fold, lambda = lambda_search, cv_opt = cv_opt)
    rr <- eigen(M_inv)
    vv <- sqrt(rr$value)
    vv[rr$value <= 0] = 0
    G <- rr$vector
    M1 <- (G) %*% diag(vv) %*% t(G)
    X1 = M1 %*% t(X)
    Xn = X1
  }

  mean_Xn = apply(Xn, 1, mean)
  # for(ii in c(1:n)){
  #   Xn[,ii] = Xn[,ii] - mean_Xn
  # }

  Wy = diag(1, p, p)

  if (p<6){
    for (k in (1:k0)){
      S = autocovm(t(Xn), k)
      Wy = Wy + (1 - k/(k0+1)) * S %*% t(S)
    }
  }else{
    for (k in (1:k0)){
      Sigma_y = autocovm(t(Xn), k)
      res = thresh(Sigma_y, Xn, mean_Xn, k, n, p, delta) #need to be done
      Sigma_ynew = t(matrix(res, p,p))
      Wy = Wy + (1 - k/(k0 + 1)) * Sigma_ynew %*% t(Sigma_ynew)

    }
  }
  rr = eigen(Wy)
  G = rr$vectors
  ev = rr$values
  # order(ev, decreasing = T)
  X1 = t(t(G) %*%Xn)
  return(list(X1 = X1, M1 = M1))

}
