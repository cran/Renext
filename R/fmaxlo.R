##=======================================================================
## Author: Yves Deville
##
## Find ML estimate of a two parameter 'maxlo' distribution using a
## sample x. The likelihood is concentrated with respect to the shape
## parameter
##
## This distribution is a special case of a beta distribution on an
## interval (0, beta)  where 'beta' is considered as an unknown parameter.
##
## It is a reparametrisation of a GPD with location mu = 0 and negative
## shape xi < 0
##
##=======================================================================

`fmaxlo` <- function(x,
                     info.observed = FALSE,
                     plot = FALSE,
                     scaleData = TRUE) {
  Cvg <- TRUE
  if (any(x <= 0)) stop("all elements in 'x' must be > 0")
  
  parnames <- c("shape", "scale")
  n <- length(x)
  M1 <- mean(x)

  if (scaleData) {
    ## xUnscaled <- x    ## useful ???
    x <- x / M1
    CV <- sqrt(1 - 1/n) * sd(x)
    cLogLik <- - n * log(M1)
    trans <- diag(c(1, 1 / M1))
    colnames(trans) <- rownames(trans) <- parnames
  } else {
    CV <- sqrt(1 - 1/n) * sd(x) / M1
  }
  
  if (CV > 1.00) stop("CV > 1. Estimation impossible for \"maxlo\"")
  if (CV > 0.999) warning("large CV value: estimation may not converge")
  
  ## could be used to find a zero
  dlogLc <- function (beta) {
    xmod <- x / beta
    R <- -mean(log(1.0 - xmod))
    dR <- -mean( xmod / (1.0 - xmod) ) / beta
    -n*( (1.0 - R) * dR / R + 1.0 / beta ) 
  }
  
  logLc <- function (beta) {
    xmod <- x/beta
    R <- -mean(log(1.0 - xmod))
    -n*( log(R) + log(beta) - R + 1.0 )
  }
  
  logL <- function (parm) {
    rho <- parm[1]
    beta <- parm[2]
    xmod <- x / beta
    R <- -mean(log(1.0 - xmod))
    n*( log(rho) - log(beta) -(rho - 1.0) * R )
  }

  ## 'beta' is scaled, as is 'x"
  log2L <- function (alpha, beta) {
    xmod <- x / beta
    xmod1 <- 1 - xmod
    RL <- mean(-log(xmod1))
    R1 <- mean(1 / xmod1)
    R2 <- mean(1 / xmod1 / xmod1)
    
    alpha1 <- alpha - 1.0
    
    logL <- n* (log(alpha/beta)  - alpha1 * RL) 

    dlogL <- c(shape = n * (1 / alpha - RL),
               scale = n * ( -1 + alpha1 * (R1 - 1) ) / beta)

    d2logL <- array(0, dim = c(2, 2), dimnames = list(parnames, parnames))
    d2logL["shape", "shape"] <- -n / alpha / alpha
    d2logL["shape", "scale"] <- n * (R1 - 1) / beta 
    d2logL["scale", "shape"] <- d2logL["shape", "scale"]
    d2logL["scale", "scale"] <-
      n * (alpha - alpha1 * R2 ) / beta / beta 

    if (scaleData) {
      logL <- logL + cLogLik
      dlogL <- trans %*% dlogL 
      d2logL <- trans %*% d2logL %*% trans
    } 
    
    list(logL = logL,
         dlogL = dlogL,
         d2logL = d2logL)    
  }
  
  ## compute an upper bound for beta. M1 <- mean(x) was done before
  M2 <- mean(x^2)
  M3 <- mean(x^3)
  
  if (scaleData) {
    betaRoots <- polyroot(c(-12*M3, 10*M3 - 6*M2, 3*(M2 - 2)))
  } else {
    betaRoots <- polyroot(c(-12*M1*M3, 10*M3 - 6*M1*M2, 3*(M2 - 2*M1^2)))
  }

  betaUpper <- max(Re(betaRoots))
  mind <- max(x)
  if (betaUpper < 2 * mind) betaUpper <- 2 * mind
  
  interv <- c(mind + 1e-7, betaUpper)
  checks <- unlist(sapply(interv, dlogLc))

  if ( (checks[1] < 0) || (checks[2] > 0) ) {
    warning("no interval found to maximise loglik")
    Cvg <- FALSE
  }
  
  res <- optimize(f = logLc, interval = interv, maximum = TRUE)
  beta.hat <- res$maximum
  alpha.hat <- -1/ mean(log(1 - x / beta.hat))
  
  res2 <- log2L(alpha = alpha.hat, beta = beta.hat)

  if (scaleData) beta.hat <- M1 * beta.hat
  
  if (alpha.hat >= 2) {
    
    if (info.observed) {
      info <- -res2$d2logL
    } else {
      info <- n * c(1 / alpha.hat / alpha.hat,
                    rep(-1 / beta.hat / (alpha.hat - 1), 2),
                    alpha.hat / (alpha.hat - 2) / beta.hat / beta.hat)
      info <- matrix(info, nrow = 2L, ncol = 2L)
      colnames(info) <- rownames(info) <- parnames
      ## print(cbind(info, -res2$d2logL)) ## for checks
    }
    
    cov <- solve(info)
    sds <- sqrt(diag(cov))
    
  } else {
    warning("'shape' is < 2 ML inference results not suitable")
    
    info <- matrix(NA, nrow = 2L, ncol = 2L)
    colnames(info) <- rownames(info) <- parnames
    ## print(cbind(info, -res2$d2logL)) ## for checks
    cov <- info
    sds <- rep(NA, 2L)
    names(sds) <- parnames
  }

  if (plot) {
    
    if (scaleData) beta.sol <- beta.hat / M1
    else beta.sol <- beta.hat

    ## prepare grid and compute the limit
    lcInf <- -n * (1 + log(mean(x)))
    betas <- seq(from = interv[1], to = interv[2], length = 200)   
    fs <- sapply(betas, logLc)
    dfs <- sapply(betas, dlogLc)
    
    ind <- 1L:length(betas)
    ind <- dfs < 20
    Stext <- ifelse(scaleData, "(scaled data)", "")

    ## now plot logLik derivative and logLik
    opar <- par(mfrow = c(2L, 1L))
    
    par(mar = c(0, 5, 5, 5))

    ## First plot: logLik derivative. 
    plot(betas[ind], dfs[ind],
         type = "n",
         main = sprintf("'Maxlo' concentrated log-lik CV = %4.2f %s", CV, Stext),
         xlab = " ", ylab = "dlogL",
         xaxt = "n", yaxt = "n",
         xlim = interv)
    axis(side = 4)
    abline(h = 0)
    abline(v = beta.sol, col = "orangered")
    abline(v = interv, col = "darkcyan", lwd = 2)
    ## abline(v = M1, col = "Chartreuse4", lty = "dotted")
    abline(h = 0, col = "gray")
    lines(betas[ind], dfs[ind],
          type = "l", lty = "solid", col = "red3")
    ## mtext(text = "mean", col = "Chartreuse4",
    ##       side = 3, at = M1, line = 0)
    par(mar = c(5, 5, 0, 5))
    
    ## Second plot: logLik function 
    plot(betas[ind], fs[ind], type = "l",
         lty = "solid", col = "red3",
         xlab = "beta (scale param.)", ylab = "logL",
         xlim = interv, ylim = range(fs[ind], lcInf))
    abline(h = lcInf, col = "orchid")
    mtext(text = "lim.", side = 4, at = lcInf,
          col = "orchid")
    abline(v = interv, col = "darkcyan", lwd = 2)
    abline(v = beta.sol, col = "orangered")
    mtext(text = "betaHat", col = "orangered",
          side = 1, at = beta.sol, line = 0.5)
    ## abline(v = M1, col = "Chartreuse4", lty = "dotted")

    par(opar)

  }

  list(estimate = c(shape = alpha.hat, scale = beta.hat),
       sd = sds,
       CV = CV,
       loglik = res2$logL,
       dloglik = res2$dlogL,
       cov = cov,
       info = info,
       cvg = Cvg)
  
}

##=======================================================================
## Author: Yves Deville
##
## Find ML estimate of a two parameter 'maxlo' distribution WITH KNOWN
## SHAPE.
##
## NB. When the SCALE 'beta' is known, -log(1 - x / beta) is exponential
## with rate 'alpha'.
##
## TODO: add the 'scaleData' argument, find bounds for the optim...
##
##=======================================================================

`fmaxlo1` <- function(x,
                      shape = 1.5,
                      plot = FALSE) {
  
  if (any(x <= 0)) stop("all elements in 'x' must be >0")  
  
  n <- length(x)
  parnames <- c("shape", "scale")

  CV <- sqrt(1 - 1/n) * sd(x) / mean(x) 
  if (CV > 0.99) warning("large CV value: estimation may not converge")

  alpha <- shape
  alpha1 <- alpha - 1
  
  ## could be used to find a zero
  dlogL1 <- function (beta) {
    xmod <- x / beta
    R <- -mean(log(1.0 - xmod))
    S1 <- mean( xmod / (1.0 - xmod) )
    n*( -1  + alpha1*S1 ) /beta 
  }
  
  logL1 <- function (beta) {
    xmod <- x/beta
    R <- -mean(log(1.0 - xmod))
    n*(log(alpha/beta) - alpha1*R)
  }
  
  mind <- max(x)
  interv <- c(mind + 1e-6, 10 * mean(x))
  checks <- unlist(sapply(interv, dlogL1))
  if ( (checks[1] < 0) ||  (checks[2] > 0) ) 
    warning("no interval found to maximise loglik")
  res <- optimize(f = logL1, interval = interv, maximum = TRUE)
  beta.hat <- res$maximum
  loglik <- res$objective
  dloglik <- dlogL1(beta.hat)
  
  if (alpha >= 2) {
    info <-  n*alpha/(alpha-2)/beta.hat/beta.hat
    cov <- solve(info)
    sds <- sqrt(diag(cov))
    
  } else {
    warning("'shape' is < 2 ML inference results not suitable")
    info <- NULL
    cov <- NULL
    sds <- NULL
  }
  
  if (plot) {
    betas <- seq(from = mind, to = 10*mean(x), length = 200)
    fs <- sapply(betas, logL1)
    dfs <- c(NA, diff(fs)/diff(betas))
    ind <- 1L:length(betas)
    ind <- dfs < 20
    plot(betas[ind], dfs[ind], type = "l", lty = "dotted")
    lines(betas[ind], sapply(betas[ind], dlogL1), col = "red")
    abline(v= mind, col = "purple")
    abline(h = 0)
    abline(v = beta.hat)
  }
  
  list(estimate = c(scale = beta.hat),
       sd = sds,
       CV = CV,
       loglik = loglik,
       dloglik = dloglik,
       cov = cov,
       info = info)
  
}

