
autocovm = function(Y, k){
  n = dim(Y)[1]
  p = dim(Y)[2]
  Y = t(Y)
  sm = rep(0, p)
  a1 = t(t(Y[, 1:(n-k)]))
  a2 = t(t(Y[, (1+k):n]))
  a1x = sweep(a1, MARGIN = 1, (apply(X = a1, MARGIN = 1, FUN = mean)), FUN = "-")
  a2x = sweep(a2, MARGIN = 1, (apply(X = a2, MARGIN = 1, FUN = mean)), FUN = "-")
  for (t in c(1:(n-k))){
    sm = sm + a1x[, t] %*% t(a2x[, t])
  }
  sm = sm/(n-k-1)
  return(sm)
}
