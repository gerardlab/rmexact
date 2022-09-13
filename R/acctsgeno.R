#'Accounting for Genotype Uncertainty using fractional counts
#'
#'@param Posterior probabilities 
#'@example
#'fraccounts(0.1,0.1,0.3,0.4,0.1, N = 100)

gamfunction <- function(y0hat, y1hat, y2hat, y3hat, y4hat) {
  n <- y0hat + y1hat + y2hat + y3hat + y4hat
result <- lgamma(n + 1) - (lgamma(y0hat + 1) + lgamma(y1hat + 1) + lgamma(y2hat + 1) + lgamma(y3hat + 1) + lgamma(y4hat +1))
return(result)
}

fraccounts <- function(gi0, gi1, gi2, gi3, gi4, N) {
  yhat = c(gi0,gi1,gi2,gi3,gi4) * N
  ykhat = sum(yhat)
  sobserved <- gamfunction(y0hat = yhat[1], 
                   y1hat = yhat[2], 
                   y2hat = yhat[3], 
                   y3hat = yhat[4], 
                   y4hat = yhat[5]) 
  umax = min(N - yhat[3] - yhat[1], yhat[2]/2, (N - yhat[3] - yhat[4])/2, yhat[5])
  umin = max(-yhat[1], (yhat[3] + yhat[2] - N)/2, -yhat[4]/2, yhat[3] + yhat[5] - N)
  m = ceiling(umax - umin) + 1
  yhat1 = c(yhat[1] + umax, 
            yhat[2] - (2 * umax), 
            yhat[3], 
            yhat[4] + (2 * umax), 
            yhat[5] - umax)
  s1 <- rep(NA, length = m)
  s1[1] <- gamfunction(y0hat = yhat1[1], 
                    y1hat = yhat1[2], 
                    y2hat = yhat1[3], 
                    y3hat = yhat1[4], 
                    y4hat = yhat1[5])
  for (i in 2:m) {
    yhat1 <- yhat1 + c(-1,2,0,-2,1)
    s1[i] <- gamfunction(y0hat = yhat1[1],
                        y1hat = yhat1[2], 
                        y2hat = yhat1[3], 
                        y3hat = yhat1[4], 
                        y4hat = yhat1[5])
  }
  sihat <- s1 - matrixStats::logSumExp(s1)
  shat <- sobserved - matrixStats::logSumExp(s1) 
  pvalue <- matrixStats::logSumExp(sihat[sihat <= shat])
  expvalue <- exp(pvalue)
  return(expvalue)
}
