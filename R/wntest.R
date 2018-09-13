#' @name wntest
#' @title Testing for multivariate or high dimensional white noise
#' @description A variety of methods to test multivariate or high-dimensional white noise,
#' including classical methods such as the multivariate portmanteau tests,
#'  Lagrange multiplier test, as well as the new method proposed by Chang, Yao and Zhou (2017)
#'  based on maximum cross correlations.
#' @details  For  a \eqn{p}-dimensional weakly stationary time
#'series \eqn{\varepsilon_t} with mean zero, denote by
#' \eqn{\Sigma(k)=\textrm{cov}(\varepsilon_{t+k},\varepsilon_t)} and
#'  \eqn{\Gamma(k) = \textrm{diag}\{{\Sigma}(0)\}^{-1/2}
#'  {\Sigma}(k)\textrm{diag}\{{\Sigma}(0)\}^{-1/2}}, respectively,
#'  the autocovariance and the autocorrelation of  at lag \eqn{k}.
#'  With the available observations
#'  \eqn{\varepsilon_1, \ldots, \varepsilon_n}, let
#'  \deqn{
#'  \widehat{\Gamma}(k) \equiv \{\widehat{\rho}_{ij}(k)\}_{1\leq i,j\leq p}
#'  =\textrm{diag}\{\widehat{\Sigma}(0)\}^{-1/2} \widehat{\Sigma}(k)\textrm{diag}\{\widehat{\Sigma}(0)\}^{-1/2}
#'  }
#'  be the sample autocorrelation matrix at lag \eqn{k}, where
#'   \eqn{\widehat{\Sigma}(k)} is the sample autocovariance matrix. Consider the hypothesis testing problem
#' \deqn{
#'  H_0:\{\varepsilon_t\}~\mbox{is white noise}~~~\textrm{versus}
#'  ~~~H_1:\{\varepsilon_t\}~\mbox{is not white noise}.
#'}
#'
#' To test the above hypothesis of multivariate or high dimensional white noise,
#' we include the traditional portmanteau tests with test statistics:
#'  \eqn{Q_{1}  = n \sum_{k=1}^K \textrm{tr}\{\widehat{\Gamma}(k)^{T}
#'     \widehat{\Gamma}(k)\}},  \eqn{Q_{2}  =
#'   n^2 \sum_{k=1}^K\textrm{tr}\{\widehat{\Gamma}(k)^{T} \widehat{\Gamma}(k)\}/(n-k)},
#'   and \eqn{Q_{3}  = n \sum_{k=1}^K
#' \textrm{tr}\{\widehat{\Gamma}(k)^{T}\widehat{\Gamma}(k)\} + p^2K(K+1)/(2n)}.
#'   Also, we include the Lagranage multiplier test  as well as the Tiao-Box likelihood ratio test.
#' For  the portmanteau tests, both \eqn{\chi^2}-approximation and normal approximation are reported.
#'
#' Since \eqn{\Gamma(k) \equiv 0} for any \eqn{k\geq1} under \eqn{H_0},
#' the newly proposed maximum cross-correlation-based test uses statistic
#' \deqn{
#'  T_n=\max_{1\leq k\leq K}T_{n,k},
#' }
#'  where \eqn{T_{n,k}=\max_{1\leq i, j\leq p}{n}^{1/2}|\widehat{\rho}_{ij}(k)|}
#'   and \eqn{K\ge 1} is prescribed. Null is rejected whenever \eqn{T_n>\textrm{cv}_\alpha},
#'   where \eqn{\textrm{cv}_\alpha >0} is the critical value determined by novel bootstrap
#'   method proposed by
#'   Chang, Yao and Zhou (2017) with no further assumptions   on the data structures.
#'
#' @param Y A \eqn{p} by \eqn{n} data matrix with \eqn{p} time series of length \eqn{n}.
#' @param k_max The largest time lag to be tested for white noise (default is \eqn{10}).
#' @param kk A vector of time lags using for test (ex. \eqn{kk = \texttt{seq}(2,10, \texttt{by} = 2)}),
#' scalar is allowed and the largest kk must be less than k_max.
#' @param M Number of bootstrap replicates, ex. \eqn{2000}.
#' @param k0 A parameter in time series PCA for pre-transformation (default is \eqn{10}).
#' @param delta The thresholding parameter in time series PCA for pre-transformation
#' (default is \eqn{1.5}).
#' @param type Tests to be performed:
#' 1 is coded for  the newly proposed maximum cross-correlation-based
#' test for high-dimensional white noise by Chang, Yao and Zhou (2017);
#' 2 is coded for the Lagrange multiplier test;
#' 3 is coded for the three portmanteau tests, where results for both
#' \eqn{\chi^2} and normal approximations are reported; and
#' 4 is coded for the Tiao-Box likelihoood ratio-based test.
#' @param alpha Level of significance (default is \eqn{0.05}).
#' @param opt Options for pre-transformation of time series. That is, one considers a transformation
#' matrix \eqn{A_{n\times p}} and corresponding pre-transformed data \eqn{AY}. For parameter
#' `opt',
#' 1 is coded for performing the transformation
#'  using package `fastclime' and user-specific tuning parameter \eqn{\lambda} for
#'  estimating the contempaneous correlations;
#' 2 is coded for performing the  transformation using the sample covariance;
#' 3 is coded for performing the  transformation using
#' package `clime' with build-in cross validation on the tuning parameter for estimating the contempaneous correlations;
#' 4 is coded for performing the transformation using
#' `fastclime' with cross validation on the tuning parameter for estimating the contempaneous correlations; and
#' else do not perform the transformation.
#' @param lambda The tuning parameter used in package `fastclime',
#' which is required for `opt=1'.
#'  The default value is \eqn{0.1}.
#' @param lambda_search The tuning parameters search for package `fastclime', which is required
#' for `opt=4' (default is \eqn{\texttt{seq}(1e-4, 1e-2, \texttt{length.out} = 50)}).
#' @param fold Number of folds used in cross validations (default is \eqn{5}).
#' @param S1 True contempaneous \eqn{p\times p} covariance matrix of the data
#'  if it is known in advance. If provided, pre-transformation will use S1 instead of options in
#'  `opt'.
#' @param cv_opt Specify which tuning parameter and the corresponding
#' estimated contempenous correlation (and the precision) matrix to be used for the pre-transformation.
#'  For example, `cv_opt = 2' will choose \eqn{\lambda} and the estimated contempenous correlation (and the precision) matrix
#'  with the second
#'   smallest cross validation error (default value is \eqn{1}, the minimun error).
#' @import fastclime
#' @import clime
#' @import foreach
#' @author Meng Cao, Wen Zhou
#' @references Chang, J., Yao, Q. and Zhou, W., 2017. Testing for high-dimensional white noise using maximum cross-correlations. Biometrika, 104(1): 111-127.
#' @references Cai, T.T., Liu, W., and Luo, X., 2011. A constrained L1 minimization approach for sparse precision matrix estimation. Journal of the American Statistical Association 106(494): 594-607.
#' @references Lutkepohl, H., 2005. New introduction to multiple time series analysis. Springer Science & Business Media.
#' @return
#' \item{res}{Test output: fail to reject (coded as \eqn{0}) or reject (coded as \eqn{1}).}
#' \item{p_value}{\eqn{p}-values or approximated \eqn{p}-value.}
#' \item{M1}{Square root  of the estimated contempenous precision matrix
#' if pre-transfermation was applied.}
#' @examples
#' library(expm)
#' p = 15
#' n = 300
#' S1 = diag(1, p, p)
#' for(ii in c(1:p)){
#' for(jj in c(1:p)){
#' S1[ii, jj] = 0.995^(abs(ii-jj))
#' }
#' }
#' S11 = sqrtm(S1)
#' X = S11 %*% matrix(rt(n*p, df = 8), ncol = n)
#' k_max = 10
#' kk = seq(2, k_max, 2)
#' M = 500
#' k0 = 10
#' delta = 1.5
#' alpha = 0.05
#' wntest(X, M, k_max, kk, type = 1, opt = 0)
#' \dontrun{
#' wntest(X, M, k_max, kk, type = 1, opt = 4, cv_opt = 1)
#'}
#' @export
wntest = function(Y, M, k_max = 10, kk, type = 1, alpha = 0.05,
                  k0 = 10, delta = 1.5, opt = 1, lambda = 0.01,
                  lambda_search = seq(1e-4, 1e-2, length.out = 50),
                  fold = 5,
                  S1 = NULL,
                  cv_opt = NULL){
  if (type == 1){
    X = Y

    if(opt %in% c(1:4)){
      X_pre = sgTs_thre(t(X), k0 = k0, delta = delta,
                      opt = opt,  lambda = lambda,
                      lambda_search = lambda_search,
                      fold = fold, cv_opt = cv_opt,
                      S1 = S1)
      X <- t(X_pre$X1)
      M1 <- X_pre$M1
    }
    bw = opbw(X)
    p = dim(X)[1]
    n = dim(X)[2]
    R0 = diag(1/sqrt(diag(X%*%t(X)/n)))
    # G0 = ginv(X%*%t(X)/n)####different than matlab code: change inv to ginv
    Tn = rep(0, k_max)
    for (k in c(1:k_max)){
      sm = matrix(0, p, p)
      for (t in c(1:(n-k))){
        sm = sm + X[, t+k]%*%t(X[, t])
      }
      sm = sm/n
      Rmk = max(max(abs(sqrt(n)*R0 %*% sm %*% R0)))
      Tn[k] = Rmk
    }
    ###initial value of Tnn
    Tnn = rep(0, length(kk))
    for (j in c(1:length(kk))){
      Tnn[j] = max(Tn[1:kk[j]])
    }

    btemp = 1/sqrt(diag((X %*% t(X))/n))
    btemp = t(btemp)
    dm = c(1, rep(0, (p^2 -1)))

    for (j in c(1:p)){
      indp = c(((j-1)*p + 1):((j-1)*p + p))
      dm[indp] = btemp*btemp[j]
    }

    W = rep(dm, k_max)

    Sn = diag(rep(1, n), n, n)
    for (i in c(1:n)){
      for (j in c(1:n)){
        if (i != j){
          xu = (i-j)/bw
          Sn[i,j] = 25/(12*pi^2*xu^2)*(sin(6*pi*xu/5)/
                                         (6*pi*xu/5)-cos(6*pi*xu/5))
        }
      }
    }




    res = p_value = rep(0, length(kk))
    for (jj in c(1:length(kk))){
      K = kk[jj]
      nt = n - K

      ft = diag(1, p^2*K, nt)
      for (tt in c(1:nt)){
        for (j in c(1:K)){
          ind = ((j-1)*p^2+1):((j-1)*p^2+p^2)
          atemp = X[,tt+j] %*% t(X[,tt])
          ft[ind, tt] = rbind(atemp)
        }
      }
      mm = apply(ft, 1, mean)
      ft = ft - mm
      cv = aeS(ft = ft, Sn = Sn[1:nt, 1:nt], W = W[1: (p^2*K)], M, alpha = 0.05)
      res[jj] = (Tnn[jj]>cv$cv)
      p_value[jj] = sum(cv$stat > Tnn[jj])/length(cv$stat)

    }
    if(opt %in% c(1:4)){
      result <- list(res = res, p_value = p_value, M1 = M1)
    }else{
      result <- list(res = res, p_value = p_value)
    }



  }else if(type == 2){
    #test_LM
    k = k_max
    p = dim(Y)[1]
    n = dim(Y)[2]
    Y = sweep(Y, MARGIN = 1, (apply(X = Y, MARGIN = 1, FUN = mean)), FUN = "-")
    YZz = c()
    for (j in (k:1)){

      YZz = rbind(YZz, Y[, j:(n-k+j)])
    }
    YZz = rbind(rep(1, n-k+1), YZz)
    nt = dim(YZz)[2]
    Yz = YZz[, 1:(nt-1)]
    Yy = Y[,(k+1):n]
    A =  (Yy %*% t(Yz)) %*% solve(Yz%*%t(Yz)) # A/B in matlab ==> B%*%solve(A)
    Z = Yy - A%*%Yz
    LM = (n - k)*(p - sum(diag(solve(Yy%*%t(Yy)) %*% (Z %*% t(Z))))) # A\B ==> solve(A)%*%B
    if (p > 25){
      Tstat1 = (LM - k*p^2)/sqrt(2*k*p^2)
      res = (abs(Tstat1) > qnorm(1- alpha/2))
      p_value <- 1 - pnorm(abs(Tstat1))

    }else{
      res = (LM > qchisq(p = 1 - alpha, df = k*p^2))
      p_value <- 1 - pchisq(LM, df = k*p^2)
    }
    result <- list(res = res, p_value = p_value)

  }else if(type == 3){
    #test_pre
    # Y = par[[1]]
    # k_max = par[[2]]
    # k = par[3]
    p = dim(Y)[1]
    n = dim(Y)[2]
    G0 = solve(Y%*%t(Y)/n)
    ll = rep(0, k_max)
    ls = rep(0, k_max)
    for (k in (1:k_max)){
      sm = matrix(0, p, p)
      for (t in (1:(n-k))){
        sm = sm + Y[, t+k] %*% t(Y[, t])
      }
      sm = sm/n
      ll[k] = sum(diag(t(sm) %*% t(G0) %*% sm %*% G0))
      ls[k] = ll[k]/(n-k)
    }
    q1 = q2 = q3 = rep(0, length(kk))
    q1n = q2n = q3n = rep(0, length(kk))
    p_value1 = p_value2 = p_value3 = rep(0, length(kk))
    p_value1n = p_value2n = p_value3n = rep(0, length(kk))

    for (ix in (1: length(kk))){
      j = kk[ix]
      Q1 = sum(ll[1:j])*n
      Q2 = sum(ls[1:j])*n^2
      Q3 = Q1 + p^2 * j * (j+1) / (2*n)
      cvK = qchisq(1-alpha, p^2*j)
      q1[ix] = Q1 > cvK
      q2[ix] = Q2 > cvK
      q3[ix] = Q3 > cvK
      p_value1[ix] = 1 - pchisq(Q1, p^2*j)
      p_value2[ix] = 1 - pchisq(Q2, p^2*j)
      p_value3[ix] = 1 - pchisq(Q3, p^2*j)

      Q1n = (Q1 - p^2*j)/sqrt(2*p^2*j)
      Q2n = (Q2 - p^2*j)/sqrt(2*p^2*j)
      Q3n = (Q3 - p^2*j)/sqrt(2*p^2*j)


      q1n[ix] = Q1n > qnorm(1-alpha);
      q2n[ix] = Q2n > qnorm(1-alpha);
      q3n[ix] = Q3n > qnorm(1-alpha);
      p_value1n[ix] = 1 - pnorm(Q1n)
      p_value2n[ix] = 1 - pnorm(Q2n)
      p_value3n[ix] = 1 - pnorm(Q3n)

    }
    res = cbind(q1,q2,q3, q1n, q2n, q3n)
    p_value <- cbind(p_value1, p_value2, p_value3,
                     p_value1n, p_value2n, p_value3n)
    result <- list(res = res, p_value = p_value)
  }else{
    #test_TB
    # Y = par[[1]]
    p = dim(Y)[1]
    n = dim(Y)[2]
    Y = sweep(Y, MARGIN = 1, apply(Y, 1, mean), FUN = "-")
    Y = rbind(t(rep(1, n)), Y)
    Yz = Y[,1:(n-1)]
    Yy = Y[, 2:n]
    A = (diag(rep(1, (p+1))) %*% Yy) %*% t(solve(Yz%*% t(Yz)) %*% Yz) # A\B ==> solve(A)%*%B

    Ytemp = A %*% Yz
    Z = Yy[2:(p+1), ] - Ytemp[2:(p+1),]

    U = max(det(X%*%t(X)/n) / det(Yy[2:(p+1),] %*% t(Yy[2:(p+1),])/n ), 0)
    Tstat = -log(U) * (n-p-3/2)
    if (p > 25){
      Tstat1 = (Tstat - p^2) / sqrt(2*p^2)
      res = (abs(Tstat1) > qnorm(1-alpha/2))
      p_value <- 1 - pnorm(abs(Tstat1))

    }else{
      res = (Tstat > qchisq(1-alpha, p^2))
      p_value <- 1 - pchisq(Tstat, p^2)
    }
    result <- list(res = res, p_value = p_value)
  }

  return(result)
}

