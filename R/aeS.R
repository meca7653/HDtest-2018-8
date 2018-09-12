

aeS = function(ft, Sn, W, M, alpha){
  p = dim(ft)[1]
  n = dim(ft)[2]
  xi = t(mvrnorm(n = M, mu = rep(0, n), Sn))/sqrt(n)
  G1 = ft %*% xi
  G = W*G1
  rrr <- apply(abs(G),2, max)
  cv = quantile(rrr, 1-alpha)
  return(list(cv = cv, stat = rrr))
}

