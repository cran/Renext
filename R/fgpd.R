##==============================================================================
## Author: Yves Deville
##
## Find ML estimate of a one parameter GPD (shape 'xi fixed')
##
## The lilelihood is easily maximised in that case
##==============================================================================

fgpd1 <- function(x, shape = 0.2, plot = FALSE) {

  xi <- shape
  names(xi) <- NULL
  if (xi <= -1.0) stop("parameter 'shape' must be > -1")
  
  n <- length(x)
  
  if (any(is.na(x)) || any(x < 0) ) stop("'x' elements must be > 0 and non NA") 
  
  ## exponential case...
  if (xi == 0.0) {
    sigma.hat <- mean(x)
    cov0 <- sigma.hat^2 / n
    return(list(estimate = c(scale = sigma.hat),
                sd = sigma.hat / sqrt(n),
                loglik = -log(sigma.hat) - n,
                cov  = cov0,
                info = 1 / cov0))
  }
  
  ## not used for now
  dlogL1 <- function (sigma) {
    xmod <- x / sigma
    ( -n + (xi + 1.0) * sum(xmod / (1.0 + xi * xmod) ) ) / sigma 
  }
  
  logL1 <- function (sigma) {
    xmod <- x / sigma
    -n * log(sigma)  - (xi + 1) * sum( log(1 + xi*xmod) ) /xi
  }

  if (xi < 0 ) {
    ## when  xi < 0, mind the support
    interv <- c(-xi*max(x), max(x))
  } else {
    ## when  xi > 0, the logL is increasing at min(x)
    ## and decreasing at max(x)
    interv <- range(x)
  }

  res <- optimize(f = logL1, interval = interv, maximum = TRUE)
  sigma.hat <- res$maximum
  
  loglik <- res$objective

  if (xi > -0.5) {
    info <- n/(1+2*xi)/sigma.hat/sigma.hat
    cov0 <- 1 / info
    sdp <- sqrt(cov0)
  } else {
    warning("'shape' is <= -0.5. No variance provided for the estimator of 'shape'") 
    info <- NULL
    cov0 <- NULL
    sdp <- NULL
  }
     
  if (plot) {
    sigmas <- seq(from = interv[1], to = max(x), length.out = 500)
    lL <- sapply(sigmas, logL1)
    plot(x = sigmas, y = lL, type = "l", cex = 0.8,
         main = sprintf("xi = %4.2f", xi))
    abline(v = res$maximum, col = "red")
    abline(v = res$maximum + c(-1, 1) * sdp, col = "red", lty = "dashed")
    ## abline(v = sigma, col = "SpringGreen3")
  }
  
  list(estimate = c(scale = sigma.hat),
       sd = sdp,
       loglik = loglik,
       cov  = cov0,
       info = info)

}
  
