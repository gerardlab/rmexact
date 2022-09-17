#'Integrating over posterior genotypes
#'
#'@param gi0,.. Posterior probabilities
#'@param B Number of samples
#'@return calculate the exact p value from this sample
#'@example
genotypelikelihood <- function(gi0,gi1,gi2,gi3,gi4, B) {
  vecpb <- rep(NA, length = B) #create a vector of length to contain all p values generated from for loop
  yb <- rep(NA, length = B)
for (i in 1:B) {
samp <- sample(x = c(0,1,2,3,4), prob = c(gi0,gi1,gi2,gi3,gi4))
yb[i] <- c(y0hat = samp[1], 
           y1hat = samp[2], 
           y2hat = samp[3], 
           y3hat = samp[4], 
           y4hat = samp[5])
n <- sum(yb)
ybk <- sum(yb)
p_value <- tetraploid(y0 = yb[1],
                      y1 = yb[2], 
                      y2 = yb[3], 
                      y3 = yb[4], 
                      y4 = yb[5])
vecpb[i] <- p_value
}
  return(vecpb)
}
