
##=====================================
## probability (distribution) function
##=====================================

pmixexp2 <- function(q,
                     prob1,
                     rate1 = 1.0,
                     rate2 = rate1 + delta,
                     delta,
                     log = FALSE) {
  
  res <- prob1*pexp(q, rate = rate1, log = FALSE) + (1-prob1)*pexp(q, rate = rate2, log = FALSE)

  if (log) res <- log(res)
  res
  
}

##==================
## density function
##==================

dmixexp2 <- function(x,
                     prob1,
                     rate1 = 1.0,
                     rate2 = rate1 + delta,
                     delta,
                     log = FALSE) {
  
  res <- prob1*dexp(x, rate = rate1, log = FALSE) +
    (1-prob1)*dexp(x, rate = rate2, log = FALSE) 
  
  if (log) res <- log(res)
  res
  
  
}

##============
## simulation
##============

rmixexp2 <- function(n,
                     prob1,
                     rate1 = 1.0,
                     rate2 = rate1 + delta,
                     delta) {
  
  g.true <- rbinom(n, size = 1, prob = prob1)
  
  x <- rep(NA, n)
  ind <- g.true == 1
  n1 <- sum(ind)
  if (n1) x[ind] <- rexp(n1, rate = rate1)
  if (n1 < n) x[!ind] <-  rexp(n-n1, rate = rate2)

  x
  
}

##=============
## hazard rate
##=============

hmixexp2 <- function(x,
                     prob1,
                     rate1 = 1.0,
                     rate2 = rate1 + delta,
                     delta) {
                       
  f <- prob1*dexp(x, rate = rate1, log = FALSE) + (1-prob1)*dexp(x, rate = rate2, log = FALSE)
  F <- prob1*pexp(x, rate = rate1, log = FALSE) + (1-prob1)*pexp(x, rate = rate2, log = FALSE)
  f/(1-F)
  
}

##=======================
## cumulater hazard rate
##=======================

Hmixexp2 <- function(x,
                     prob1,
                     rate1 = 1.0,
                     rate2 = rate1 + delta,
                     delta) {
  
  F <- prob1*pexp(x, rate = rate1, log = FALSE) +
    (1-prob1)*pexp(x, rate = rate2, log = FALSE)
  
  -log(1-F)
  
}

##=======================
## quantile function
##=======================

qmixexp2 <- function(p,
                     prob1,
                     rate1 = 1.0,
                     rate2 = rate1 + delta,
                     delta) {

  tol <- 1e-9
  n <- length(p) 

  if (rate1 > rate2) {
    ## warning("using the 'qmixep2' function with rate2 < rate1") 
    prob1 <- 1-prob1
    ratep <- rate1
    rate1 <- rate2
    rate2 <- ratep
  }
  
  Hs <- - log(1-p)
  lprob1 <- log(prob1)

  xs <- rep(NA, n)
  nit <- rep(NA, n)
  x.L.prec <- -Inf
  rate.bar <- prob1*rate1 + (1-prob1)*rate2
  
  for (i in 1:n) {
    
    H.star <- Hs[i]

    x.L <-  max(c(H.star/rate.bar, (lprob1 + H.star)/rate1))
    x.U <-  H.star/rate1
    
    H.L <- Hmixexp2(x.L, prob1 = prob1, rate1 = rate1, rate2 = rate2) 
    H.U <- Hmixexp2(x.U, prob1 = prob1, rate1 = rate1, rate2 = rate2) 
    h.L <- hmixexp2(x.L, prob1 = prob1, rate1 = rate1, rate2 = rate2)
    
    dx <- x.U - x.L
    
    cvg <- FALSE
    iter <- 1

    NR <- FALSE
    
    while ( !cvg && (iter < 20) ) {
    
      if (!NR) {
        
        a <- ( (H.U - H.L) / dx - h.L ) / dx
        
        if ( is.na(a) || (a >= -1e-6) ) {
          ## cat(sprintf("i = %d  it. = %d x.L = %7.4f
          ## x.U = %7.4f dx = %7.4f H.star = %10.7f a = %7.4f\n",
          ## i, iter, x.L, x.U, dx, H.star, a))
          ## warning("a >= 0")
          NR <- TRUE
          x.prov <- x.L + (H.star - H.L) / h.L
        } else {
          Delta <- h.L*h.L - 4*a*(H.L-H.star)
          x.prov <- x.L + ( -h.L + sqrt(Delta) ) / 2 / a          
        }
          
      } else {
        x.prov <- x.L + (H.star - H.L) / h.L
      }
        
      H.prov <- Hmixexp2(x.prov,
                         prob1 = prob1,
                         rate1 = rate1,
                         rate2 = rate2)
      
      if (abs(H.prov - H.star) < tol) {
        x.star <- x.prov
        cvg <- TRUE
      } else {
    
        if (H.prov < H.star) {
          x.L <- x.prov
          h.L <- hmixexp2(x.prov, prob1 = prob1, rate1 = rate1, rate2 = rate2)
          H.L <- H.prov 
        } else {
          x.U <- x.prov
          H.U <- H.prov       
        }
        
        dx <- x.U - x.L
      }
      
      iter <- iter + 1

    }
    
    if (cvg) xs[i] <- x.star
    else {
      print(NR)
      cat(sprintf("i = %d  it. = %d\n x.L = %e x.U = %e dx = %e\n
                  H.L= %e H.U= %e H.star = %e a = %e\n",
                  i, iter, x.L, x.U, dx, H.L, H.U, H.star, a))
      cat(sprintf("prob1 = %e  rate1 = %e rate2 = %e\n",
             prob1, rate1, rate2))
      cat(sprintf("p = %e\n",
                  p[i]))
      stop("diverged i =", i, "\n")
    }
    nit[i] <- iter
  }

  attr(xs, "nit") <- nit
  ## print(nit)
  xs
  
}
