##======================================================
## Y is logexp if log X - log threshold is
## exponential. This is the shifted Pareto distribution
##======================================================

rspareto <- function(n, delta = 1.0, shape = 4.0) {
  
  if (delta <= 0) stop("'delta' must be > 0")
  if (shape <= 0) stop("'shape' must be > 0")
  X <- (runif(n))^(-1/shape) 
  X <- delta*X
  X - delta

}

dspareto <- function(x, delta = 1.0, shape = 4.0, log = FALSE) {
  
  if (delta <= 0) stop("'delta' must be >0")
  if (shape <= 0) stop("'shape' must be > 0")
  f <- rep(0, length(x))
  ind <- x > 0
  f[ind] <- log(shape) - (shape + 1)*log(delta + x[ind]) + shape*log(delta)
  if (!log) f[ind] <- exp(f[ind])
  f
  
}

pspareto <- function(q, delta = 1.0, shape = 4.0, lower.tail = FALSE) {

  if (delta <= 0) stop("'delta' must be >0")
  if (shape <= 0) stop("'shape' must be  >0")
  S <- rep(0, length(q))
  ind <- q > 0 
  S[ind] <- (delta / (delta + q[ind]) )^shape
  if (lower.tail) S[!ind] <- 1
  else S[ind] <- 1 - S[ind]
  S
  
}

qspareto <- function(p, delta = 1.0, shape = 4.0) {

  if (delta <= 0) stop("'delta' must be >0")
  if (shape <= 0) stop("'shape' must be > 0")
  if ( any(p < 0) || any(p > 1) ) stop("'p' must be >=0 and <=1")
  lr <- -log(1-p) / shape  
  delta* (exp(lr) -1)

}

##=====================================================
## Shifted Left Truncated Weibull or "SLTW"
##=====================================================

rSLTW <- function(n, delta = 1.0, shape = 1.0, scale = 1.0) {

  if (delta <= 0) stop("'delta' must be >0")
  if (shape <= 0) stop("'shape' must be >0")
  if (scale <= 0) stop("'scale' must be >0")

  U <- runif(n)
  X <- ( delta^shape - (scale^shape)*log(U) )^(1/shape) - delta
  
}

dSLTW <- function(x, delta = 1.0, shape = 1.0, scale = 1.0, log = FALSE) {
  
  if (delta <= 0) stop("'delta' must be >0")
  if (shape <= 0) stop("'shape' must be >0")
  if (scale <= 0) stop("'scale' must be >0")
  
  f <- rep(0, length(x))
  ind <- x > 0
  xs <- ( x[ind] + delta ) / scale 
  f[ind] <- exp( - xs^shape + (delta/scale)^shape )*(xs ^(shape - 1)) * shape / scale
  
  if (log) { f[ind] <- log(f[ind]) }
  f
  
}

pSLTW <- function(q, delta = 1.0, shape = 1.0, scale = 1.0,  lower.tail = FALSE) {
  
  if (delta <= 0) stop("'delta' must be >0")
  if (shape <= 0) stop("'shape' must be >0")
  if (scale <= 0) stop("'scale' must be >0")

  S <- rep(0, length(q))
  ind <- q > 0 
  S[ind] <- exp( - ((q[ind] +delta)/scale)^shape + (delta/scale)^shape )
  if (!lower.tail) S[ind] <- 1 -  S[ind]
  else S[!ind] <- 1
  S
  
}

qSLTW <- function(p, delta = 1.0,  shape = 1.0, scale = 1.0) {
  
  if (delta <= 0) stop("'delta' must be >0")
  if (shape <= 0) stop("'shape' must be >0")
  if (scale <= 0) stop("'scale' must be >0")
  
  if ( any(p < 0) || any(p > 1) ) stop("'p' must be >=0 and <=1")

  quant <- delta^shape - (scale^shape) * log(1-p)
  quant <- quant^(1/shape) - delta
  
}
