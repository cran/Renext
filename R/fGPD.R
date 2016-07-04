##**********************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
## 
## Maximum Likelihood estimation of a Generalised Pareto Distibution
## by likelihood concentration.
##
## TODO : situations with one known parameter
##**********************************************************************

fGPD <- function(x,
                 info.observed = TRUE,
                 shapeMin = -0.8,
                 dCV = 1e-4,
                 cov = TRUE,
                 trace = 0) {
  
    parnames <- c("shape", "scale")
    n <- length(x)
    M1 <- mean(x)
    CV <- sqrt(1.0 - 1.0 / n) * sd(x) / M1
    
    if (any(x <= 0)) stop("all elements in 'x' must be > 0")
    if (shapeMin >= 0) stop("'shapeMin' must be negative")
    if (shapeMin < -1) warning("'shapeMin < -1': non-identifiable model")
    if (dCV >= 1e-2) warning("'dCV' should be small < 0.01")
    
    if (CV > 1.0 + dCV) {
        ## Lomax case
        res.lomax <- flomax(x,
                            info.observed = info.observed, cov = cov)
        res.gpd <- lomax2gpd(parLomax = res.lomax$estimate,
                             vcovLomax = res.lomax$cov)
        loglik <- res.lomax$loglik
        if (!is.null(vcov <- attr(res.gpd, "vcov"))) {
            dloglik <- attr(res.gpd, "jacobian")[ , c("shape", "scale")] %*%
                res.lomax$dloglik
            sd.gpd <- sqrt(diag(vcov))
        } else {
            dloglik <- NULL
            sd.gpd <- NULL
        }
    } else if (CV < 1.0 - dCV) {
        ## Maxlo case
        res.maxlo <- fmaxlo(x, shapeMin = - 1.0 / shapeMin,
                            info.observed = info.observed, cov = cov)
        res.gpd <- maxlo2gpd(parMaxlo = res.maxlo$estimate,
                             vcovMaxlo = res.maxlo$cov)
        loglik <- res.maxlo$loglik
        if (!is.null(vcov <- attr(res.gpd, "vcov"))) {
            dloglik <- attr(res.gpd, "jacobian")[ , c("shape", "scale")] %*%
                res.maxlo$dloglik
            sd.gpd <- sqrt(diag(vcov))
        } else {
            dloglik <- NULL
            sd.gpd <- NULL
        }
    } else {
        ## exponential case. Note that the variance of the estimated
        ## scale is greater than that of the one-parameter exponential
        ## distribution. 
        if (trace) {
            cat(sprintf("CV = %6.4f close to 1: exponential\n", CV))
        }
        res.gpd <- c("scale" = M1, "shape" = 0.0)
        loglik <- -n * (1 + log(M1)) 
        if (cov) {
            dloglik <- c("scale" = 0.0, "shape" = 0.0)
            if (info.observed) {
                M3u <- mean((x / M1)^3)
                vcov <- c(M1 * M1 * 2 * (M3u - 3.0) / 3.0,
                          -M1, -M1, 1.0) /
                    (M3u * 2.0 / 3.0 - 3.0) / n
                vcov <- matrix(vcov, nrow = 2L, ncol = 2L)
                colnames(vcov) <- rownames(vcov) <- c("shape", "scale")
                sd.gpd <- sqrt(diag(vcov))
            } else {
                vcov <- matrix(c(2.0 * M1^2, -M1, -M1, 1.0) / n,
                               nrow = 2L, ncol = 2L)
                colnames(vcov) <- rownames(vcov) <- c("shape", "scale")
                rn <- sqrt(n)
                sd.gpd <- c("scale" = M1 * sqrt(2.0) / rn,  shape = 1.0 / rn)
            }
        } else {
            return(list(estimate = res.gpd[1L:2L], CV = CV,
                        loglik = loglik, cvg = TRUE))
        }
        
    }
    
    list(estimate = res.gpd[1L:2L],
         CV = CV,
         loglik = loglik,
         dloglik = dloglik,
         sd = sd.gpd,
         cov = vcov,
         cvg = TRUE)
    
}

##==============================================================================
## Author: Yves Deville
##
## Find ML estimate of a one parameter GPD (shape 'xi' fixed)
##
## The lilelihood is easily maximised in that case
##==============================================================================

fgpd1 <- function(x,
                  shape = 0.2,
                  info.observed = TRUE,
                  plot = FALSE) {

    Cvg <- TRUE
    
    xi <- shape
    names(xi) <- NULL
    if (xi <= -1.0) stop("parameter 'shape' must be > -1")
    
    n <- length(x)
    CV <- sqrt(1.0 - 1.0 / n) * sd(x) / mean(x)
    
    if (any(is.na(x)) || any(x < 0)) stop("'x' elements must be > 0 and non NA") 
    
    ## exponential case...
    if (xi == 0.0) {
        sigma.hat <- mean(x)
        cov0 <- sigma.hat^2 / n
        return(list(estimate = c(scale = sigma.hat),
                    sd = sigma.hat / sqrt(n),
                    CV = CV,
                    loglik = -log(sigma.hat) - n,
                    cov  = cov0,
                    info = 1.0 / cov0,
                    cvg = TRUE))
    }
    
    ## not used for now
    dlogL1 <- function (sigma) {
        xmod <- x / sigma
        (-n + (xi + 1.0) * sum(xmod / (1.0 + xi * xmod))) / sigma 
    }
    
    logL1 <- function (sigma) {
        xmod <- x / sigma
        -n * log(sigma)  - (xi + 1.0) * sum(log(1 + xi * xmod)) / xi
    }
     
    log2L1 <- function (sigma) {
        xi1  <- xi + 1.0
        xmod <- x / sigma
        A <- mean(log(1 + xi * xmod))
        B <- mean(xmod / (1.0 + xi * xmod))
        B2 <-  mean((xmod / (1.0 + xi * xmod))^2)
        
        logL <- -n * (log(sigma)  +  xi1 * A / xi)
        dlogL <-  -n * (1.0  -  xi1 * B) / sigma
        d2logL <- n * (1.0 - 2.0 * xi1 * B + xi * xi1 * B2) / sigma / sigma
        
        list(logL = logL,
             dlogL = dlogL,
             d2logL = d2logL)
    }

    if (xi < 0) {
        ## when  xi < 0, mind the support
        interv <- c(-xi * max(x), max(x))
    } else {
        ## when  xi > 0, the logL is increasing at min(x)
        ## and decreasing at max(x)
        interv <- range(x)
    }
    
    res <- optimize(f = logL1, interval = interv, maximum = TRUE)
    sigma.hat <- res$maximum
    
    loglik <- res$objective
    res2 <- log2L1(sigma = sigma.hat)
    
    if (xi > -0.5) {
        if (info.observed) {
            info <- -res2$d2logL
        } else {
            info <- n / (1.0 + 2.0 * xi) / sigma.hat / sigma.hat
        }
        cov0 <- 1.0 / info
        sdp <- sqrt(cov0)
        
    } else {
        warning("'shape' is <= -0.5. No variance provided for the ",
                "estimator of 'shape'") 
        info <- NA
        cov0 <- NA
        sdp <- NA
    }
    
    if (plot) {
        sigmas <- seq(from = interv[1], to = max(x), length.out = 500)
        lL <- sapply(sigmas, logL1)
        plot(x = sigmas, y = lL, type = "l", cex = 0.8,
             main = sprintf("xi = %4.2f", xi))
        abline(v = res$maximum, col = "red")
        abline(v = res$maximum + c(-1.0, 1.0) * sdp, col = "red", lty = "dashed")
        ## abline(v = sigma, col = "SpringGreen3")
    }
    
    list(estimate = c(scale = sigma.hat),
         sd = sdp,
         CV = CV,
         loglik = loglik,
         dloglik = res2$dlogL,
         cov  = cov0,
         info = info,
         cvg = Cvg)
    
}

##==============================================================================
## Author: Yves Deville
##
## Find ML estimate of a one parameter GPD (shape 'xi' fixed)
##
## The lilelihood is easily maximised in that case
##==============================================================================

fGPD1 <- function(x,
                  shape = 0.2,
                  info.observed = TRUE,
                  cov = TRUE,
                  plot = FALSE) {

    Cvg <- TRUE
    
    xi <- unname(shape)
    if (xi <= -1.0) stop("parameter 'shape' must be > -1")

    n <- length(x)
    if (any(is.na(x)) || any(x < 0)) stop("'x' elements must be > 0 and non NA")
    
    if (plot) cov <- TRUE

    if (!cov) {
        cov0 <- sdp <- info <- NULL
    }
    
    if (xi > 0.0) {
        
        res.lomax1 <- flomax1(x, shape = 1.0 / xi,
                              info.observed = info.observed,
                              cov = cov,
                              plot = FALSE)
        
        sigma.hat <- res.lomax1$estimate * xi
        loglik <- res.lomax1$loglik
        CV <- res.lomax1$CV
        Cvg <- res.lomax1$cvg
        dloglik <- res.lomax1$dloglik / xi 
        
        if (cov) {
            cov0 <- res.lomax1$cov * xi^2
            info <- 1.0 / cov0
            sdp <- sqrt(cov0)
        } 
        
     
    } else if (xi < 0.0) {
        
        res.maxlo1 <- fmaxlo1(x, shape = -1.0 / xi,
                              info.observed = info.observed,
                              cov = cov,
                              plot = FALSE)
        
        sigma.hat <- -res.maxlo1$estimate * xi
        loglik <- res.maxlo1$loglik
        CV <- res.maxlo1$CV
        Cvg <- res.maxlo1$cvg
        dloglik <- -res.maxlo1$dloglik / xi 
        
        if (cov) {   
            if (xi <= -0.5) {
                warning("'shape' is <= -0.5. No variance provided for the ",
                        "estimator of 'shape'") 
                info <- NA
                cov0 <- NA
                sdp <- NA
            }  else {
                cov0 <- res.maxlo1$cov * xi^2
                info <- 1.0 / cov0
                sdp <- sqrt(cov0)
            }
        }
        
    } else if (xi == 0.0) {
        
        sigma.hat <- mean(x)
        CV <- sqrt(1.0 - 1.0 / n) * sd(x) / sigma.hat
        
        loglik <- -n * (log(sigma.hat) + 1.0)
        dloglik <- 0.0
        
        if (cov) {
            cov0 <- sigma.hat^2 / n
            info <- cov0
        }
        
    }
    
    if (plot) {
        
        ## not used for now
        dlogL1 <- function(sigma) {
            xmod <- x / sigma
            (-n + (xi + 1.0) * sum(xmod / (1.0 + xi * xmod))) / sigma 
        }
        
        logL1 <- function(sigma) {
            xmod <- x / sigma
            -n * log(sigma)  - (xi + 1.0) * sum(log(1 + xi * xmod)) / xi
        }
        
        ## mind the support
        if (xi < 0) interv <- c(-xi * max(x), max(x))
        else interv <- range(x)
        
        sigmas <- seq(from = interv[1], to = max(x), length.out = 500)
        lL <- sapply(sigmas, logL1)
        plot(x = sigmas, y = lL, type = "l", cex = 0.8,
             col = "orangered", lwd = 2,
             xlab = "sigma", ylab = "logL",
             main = sprintf("GPD log-lik for xi = %4.2f", xi))
        abline(v = sigma.hat, col = "SpringGreen3")
        abline(v = sigma.hat + c(-1.0, 1.0) * sdp,
               col = "SpringGreen1", lty = "dashed")
        ## abline(v = sigma, col = "SpringGreen3")
    }
    
    list(estimate = c(scale = sigma.hat),
         sd = sdp,
         CV = CV,
         loglik = loglik,
         dloglik = dloglik,
         cov  = cov0,
         info = info,
         cvg = Cvg)
    
}

