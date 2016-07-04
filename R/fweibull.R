##==============================================================================
## Author: Yves Deville
##
## Find ML estimate of a two-parameters Weibull distribution using a sample x.
##
## The parameter 'eta' = beta^alpha is used in place of beta (= scale) where
## alpha = shape. Then 'eta' is concentrated out of the likelihood.
##==============================================================================

"fweibull" <- function(x,
                       info.observed = TRUE,
                       scaleData = TRUE,
                       cov = TRUE,
                       check.loglik = FALSE) {
    
    if (any(x <= 0)) stop("all elements in 'x' must be > 0")  
    
    n <- length(x)
    xbar <- mean(x)
    parnames <- c("shape", "scale")
    
    if (scaleData) {
        x <- x / xbar
        cLogLik <- - n * log(xbar)
    } else {
        cLogLik <- 0.0
    }
    
    lx <- log(x)
    mlx <- mean(lx)
    alpha <- 1.2825 / sd(lx)
    
    phi <- function(alpha) {
        if (alpha < 1e-8) stop("'alpha' must be > 0") 
        xa <- x^alpha
        mean(xa * lx) / mean(xa) -mlx - 1.0 / alpha
    }
    
    res <- uniroot(f = phi, interval = c(1e-8, 6 * alpha),
                   tol = .Machine$double.eps^0.3)
    
    if (abs(res$f.root) > 0.0001) {
        cat("fweibull estimation\n")
        print(res)
        stop("root not found in 'fweibull'")
    }
    
    alpha.hat <- res$root
    
    xa <- x^alpha.hat
    eta.hat <- mean(xa)
    leta <- log(eta.hat)
    
    loglik <- n * (log(alpha.hat) - leta + (alpha.hat - 1.0) * mlx - 1.0) +
        cLogLik
    estimate <- c("shape" = alpha.hat, "scale" = eta.hat^(1.0 / alpha.hat))
    
    if (cov) {
    
        if (info.observed) {
            r0 <- eta.hat      ## = mean(xa)
            r1 <- mean(xa * lx)
            r2 <- mean(xa * lx * lx)
            
            I11 <- n * (1.0 / alpha.hat / alpha.hat + r2 / eta.hat)
            I12 <- -n * r1 / eta.hat / eta.hat
            I22 <- n * (2 * r0 / eta.hat -1) / eta.hat / eta.hat
            
            info <- matrix(c(I11, I12, I12, I22), nrow = 2L, ncol = 2L)
        } else {
            Euler <- -digamma(1.0)   ## 0.577216
            lambda1 <- 1.0 - Euler
            lambda2 <- pi * pi / 6 + Euler^2 - 2 * Euler
            I11 <-  n * ( 1.0 + (leta * leta + 2 * lambda1 * leta + lambda2) ) /
                alpha.hat / alpha.hat 
            I12 <- -n * (lambda1 + leta) / eta.hat / alpha.hat
            I22 <- n * (1.0) / eta.hat / eta.hat
            info <- matrix(c(I11, I12, I12, I22), nrow = 2L, ncol = 2L)
        }

        colnames(info) <- rownames(info) <- c("shape", "eta")
        
        ## derivatives of beta = scale with respect to 'alpha' and 'eta'
        dalpha <- -leta * eta.hat^(1.0 / alpha.hat) / alpha.hat/ alpha.hat
        deta <- eta.hat^(1.0 / alpha.hat - 1.0) / alpha.hat
        
        mat <- solve(info)
        sd.eta = sqrt(mat[2L, 2L])
        L <- matrix(c(1, dalpha, 0, deta), nrow = 2L, ncol = 2L)
        cov.beta <- L %*% mat %*% t(L)

        if (scaleData) {
            d <- xbar^alpha.hat
            estimate["scale"] <- estimate["scale"] * xbar
            eta.hat <- eta.hat * d
            info[ , 2L] <- info[ , 2L] / d
            info[2L, ] <- info[2L, ] / d
            sd.eta <- sd.eta * d
            cov.beta[ , 2L] <- cov.beta[ , 2L] * xbar
            cov.beta[2L, ] <- cov.beta[2L, ] * xbar
        }
        
        sdp <- sqrt(diag(cov.beta))
        names(sdp) <- parnames

        
    } else {
        sdp <- NULL
        sd.eta <- NULL
        info <- NULL
        cov.beta <- NULL
    }
        
    if (check.loglik) {
        check.loglik <-
            sum(stats::dweibull(x = x,
                                shape = estimate["shape"],
                                scale = estimate["scale"],
                                log = TRUE))
    } else {
        check.loglik <- FALSE
    }
    
    list(estimate = estimate,
         sd = sdp,
         loglik = loglik,
         check.loglik = check.loglik,
         cov  = cov.beta,
         eta = eta.hat,
         sd.eta = sd.eta,
         info = info)
    
}


