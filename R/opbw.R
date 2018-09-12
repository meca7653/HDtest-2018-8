

opbw = function(X){
  p = dim(X)[1]
  n = dim(X)[2]

  rho = rep(0, p)
  s2 = rep(0, p)
  for (i in c(1:p)){
    atemp = X[i, 1:(n-1)]
    rho[i] = t(X[i, 2:n]) %*% (atemp)/(t(atemp) %*% (atemp))
    btemp = X[i, 2:n] - rho[i] * X[i, 1:(n-1)]
    s2[i] = t(btemp) %*% btemp/(n-1)
  }
  num = sum(4*rho^2 * s2^2/((1-rho)^8))
  den = sum(s2^2/(1 - rho)^4)
  op_bw = 1.3221*(num/den*n)^(1/5)

  return(op_bw)

}

