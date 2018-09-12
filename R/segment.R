
segment = function(Y, mean_y, k, n, p){
  b = rep(0, p^2)
  for (t in (1:(n-k))){
    s = 0
    C = Y[,(t+k)] - mean_y
    D = Y[, t] - mean_y
    for (i in (1:p)){
      for (j in (1:p)){
        s = s+1
        b[s] = b[s] + C[i]*D[j]
      }
    }
  }
  b = b/n
  return (b)
}
