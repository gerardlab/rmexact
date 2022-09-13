#'Accounting for Genotype Uncertainty using fractional counts
#'
#'
#' 
gamfunction <- function(y0hat,y1hat,y2hat,y3hat,y4hat) {
  n <- y0hat + y1hat + y2hat + y3hat + y4hat
result <- lgamma(n + 1) - (lgamma(y0hat + 1) + lgamma(y1hat + 1) + lgmamma(y2hat + 1) + lgamma(y3hat + 1) + lgamma(y4hat +1))
return(result)
}

fraccounts <- function(gi0,gi1,gi2,gi3,gi4) {
  yhat = c(gi0,gi1,gi2,gi3,gi4) * N
  ykhat = sum(yhat)
  sobserved <- gamfunction(y0hat = yhat[1], 
                   y1hat = yhat[2], 
                   y2hat = yhat[3], 
                   y3hat = yhat[4], 
                   y4hat = yhat[5]) 
  umax = min(n - y2hat - y0hat, y1hat/2, (n - y2hat - y3hat)/2, y4hat)
  umin = max(-y0hat, (y2hat + y1hat - n)/2, -y3hat/2, y2hat + y4hat - n)
  m = ceiling(umax - umin) + 1
  yhat1 = c(y0hat + umax, 
            y1hat - (2 * umax), 
            y2hat, 
            y3hat + (2 * umax), 
            y4hat - umax)
  s1[1] <- gamfunction(y0hat = yhat1[1], 
                    y1hat = yhat1[2], 
                    y2hat = yhat1[3], 
                    y3hat = yhat1[4], 
                    y4hat = yhat1[5])
  for (i in 2:m) {
    yhat1 <- yi + c(-1,2,0,-2,1)
    s1[i] <- gamfunction(y0hat = yhat1[1],
                        y1hat = yhat1[2], 
                        y2hat = yhat1[3], 
                        y3hat = yhat1[4], 
                        y4hat = yhat1[5])
  }
  sihat <- s1 - lse(s1)
  shat <- sobserved - lse(s) 
  pvalue <- lse(sihat[sihat <= shat])
  expvalue <- exp(pvalue)
}

