##======================================================================
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
##
## Find ML estimate of a two parameter Lomax distribution using a
## sample 'x'. The likelihood is concentrated with respect to the
## shape parameter 'alpha' and maximised with respect to 'beta'.
##
## When 'scaleData' is TRUE, the sample is divided by its mean in order
## to make the optimisation safer.
##
## The 'flomax' function was changed on 2013-05-18 to implement
## a new optimisation interval determination as well as a data
## scaling. The plot was also improved, yet remaining an utility
## plot.
## 
##=======================================================================

`flomax` <- function(x,
                     info.observed = TRUE,
                     plot = FALSE,
                     scaleData = TRUE,
                     cov = TRUE) {

    Cvg <- TRUE
    if (any(x <= 0.0)) stop("all elements in 'x' must be > 0")  
    
    parnames <- c("shape", "scale")
    n <- length(x)
    M1 <- mean(x)
    
    if (scaleData) {
        x <- x / M1
        CV <- sqrt(1.0 - 1.0 / n) * sd(x)
        cLogLik <- - n * log(M1)
        trans <- diag(c(1.0, 1.0 / M1))
    } else {
        CV <- sqrt(1.0 - 1.0 / n) * sd(x) / M1
    }
    
    if (CV < 1.00) stop("CV < 1. Estimation impossible for \"lomax\"")
    
    ## compute an upper bound for beta. M1 <- mean(x) was done before
    M2 <- mean(x^2)
    M3 <- mean(x^3)
    
    if (scaleData) {
        ## 1st condition ensures ell(beta) > ell(infty) the second
        ## is better
        ## betaRoots <- polyroot(c(2*M3, 2*M3 - 3*M2, -3*(M2 - 2)))
        betaRoots <- polyroot(c(M3, M3 - M2, 1.0 - M2 / 2))
    } else {
        ## 1st condition ensures ell(beta) > ell(infty), the second
        ## is better
        ## betaRoots <- polyroot(c(2*M1*M3, 2*M3 - 3*M2*M1, -3*(M2 - 2*M1^2))) 
        betaRoots <- polyroot(c(M1 * M3, M3 - M1 * M2, M1^2 - M2 / 2))
    }
    
    betaLower <- 0.15 * min(x)
    betaUpper <- max(Re(betaRoots))
    mind <- max(x)
    if (betaUpper < mind) betaUpper <- mind
    
    interv <- c(betaLower, betaUpper)
    
    ## ## could be used to find a zero. Yet we prefer 'optimise' here
    ## dlogLc <- function (beta) {
    ##     xmod <- x / beta
    ##     R <- mean(log(1.0 + xmod))
    ##     dR <- -mean( xmod / (1.0 + xmod) ) / beta
    ##     -n * ((R + 1.0) * dR/R + 1.0 / beta) 
    ## }
    
    logLc <- function (beta) {
        xmod <- x / beta
        R <- mean(log(1.0 + xmod))
        -n * (log(R) + log(beta) + R + 1.0)
    }
    
    if (plot) cov <- TRUE
    
    if (cov) {
        
        log2L <- function (alpha, beta) {
            xmod <- x / beta
            xmod1 <- 1.0 + xmod
            s1 <- sum(log(xmod1))
            s2 <- sum(xmod / xmod1)
            alpha1 <- 1.0 + alpha
            
            logL <- n * log(alpha / beta)  - alpha1 * s1 
            
            dlogL <- c(shape = n / alpha - s1,
                       scale = (-n + alpha1 * s2) / beta)
            
            d2logL <- array(0, dim = c(2L, 2L), dimnames = list(parnames, parnames))
            d2logL["shape", "shape"] <- -n / alpha / alpha
            d2logL["shape", "scale"] <- s2 / beta
            d2logL["scale", "shape"] <- d2logL["shape", "scale"]
            d2logL["scale", "scale"] <-
                (n  - alpha1 * s2  -alpha1 * sum(xmod / xmod1 / xmod1)) / beta / beta
            
            if (scaleData) {
                logL <- logL + cLogLik
                dlogL <- trans %*% dlogL 
                d2logL <- trans %*% d2logL %*% trans
                ## added on 2015-12-17 Yves
                rownames(dlogL) <- parnames
                dimnames(d2logL) <- list(parnames, parnames)
            } 
            
            list(logL = logL,
                 dlogL = dlogL,
                 d2logL = d2logL)
        }

    } 
    
    res <- optimize(f = logLc, interval = interv, maximum = TRUE)
    beta.hatS <- res$maximum
    alpha.hat <- 1.0 / mean(log(1.0 + x / beta.hatS))

    if (scaleData) beta.hat <- M1 * beta.hatS
    else beta.hat <- beta.hatS
    
    if (!cov) {
        loglik <- res$objective
        if (scaleData) loglik <- loglik + cLogLik
        return(list(estimate = c(shape = alpha.hat, scale = beta.hat),
                    CV = CV,
                    loglik = loglik,
                    cvg = Cvg))
    }

    ## Remind that 'log2L' is called with scaled 'beta'. The
    ## derivatives are corrected for the scaling.
    res2 <- log2L(alpha = alpha.hat, beta = beta.hatS)
    
    if (info.observed) {
        info <- -res2$d2logL
        vcov <- try(solve(info))
        ## changed 2015-12-17: the test was wrong!!!
        if (inherits(vcov, "try-error")) {
            vcov <- NULL
            sds <- NULL
            warning("hessian could not be inverted")
        } else {
            sds <- sqrt(diag(vcov))
        }
    } else {
        a1 <- alpha.hat + 1.0
        a2 <- alpha.hat + 2.0
        i11 <- 1.0 / alpha.hat / alpha.hat
        i12 <- -1.0 / beta.hat / a1
        i22 <- alpha.hat / a2 / beta.hat / beta.hat
        info <- matrix(n * c(i11, i12, i12, i22), nrow = 2L, ncol = 2L)
        colnames(info) <- rownames(info) <- parnames
        c11 <- alpha.hat^2 * a1^2
        c12 <- alpha.hat * a1 * a2 * beta.hat
        c22 <- a1^2 * a2 * beta.hat^2 / alpha.hat
        vcov <- matrix(c(c11, c12, c12, c22) / n, nrow = 2L, ncol = 2L) 
        colnames(vcov) <- rownames(vcov) <- parnames
        sds <- sqrt(diag(vcov))
    }
    
    if (plot) {
       
        dlogLc <- function (beta) {
            xmod <- x / beta
            R <- mean(log(1.0 + xmod))
            dR <- -mean(xmod / (1.0 + xmod)) / beta
            -n * ((R + 1.0) * dR / R + 1.0 / beta) 
        }
        
        if (scaleData) beta.sol <- beta.hat / M1
        else beta.sol <- beta.hat
        
        ## prepare grid and compute the limit
        lcInf <- -n * (1 + log(mean(x)))
        betas <- seq(from = interv[1], to = interv[2], length = 200)   
        fs <- sapply(betas, logLc)
        dfs <- sapply(betas, dlogLc)
        
        ind <- 1L:length(betas)
        ind <- (dfs < 20)
        Stext <- ifelse(scaleData, "(scaled data)", "")
        
        ## now plot logLik derivative and logLik
        opar <- par(mfrow = c(2L, 1L))
        
        par(mar = c(0, 5, 5, 5))
        ## First plot: logLik derivative. 
        plot(betas[ind], dfs[ind],
             type = "n",
             main = sprintf("'Lomax' concentrated log-lik CV = %4.2f %s", CV, Stext),
             xlab = " ", ylab = "dlogL",
             xaxt = "n", yaxt = "n",
             xlim = interv)
        axis(side = 4)
        abline(h = 0)
        abline(v = beta.sol, col = "orangered")
        abline(v = interv, col = "darkcyan", lwd = 2)
        abline(v = M1, col = "Chartreuse4", lty = "dotted")
        abline(h = 0, col = "gray")
        lines(betas[ind], dfs[ind],
              type = "l", lty = "solid", col = "red3")
        
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
        abline(v = M1, col = "Chartreuse4", lty = "dotted")
        par(opar)
    }
    
    list(estimate = c(shape = alpha.hat, scale = beta.hat),
         sd = sds,
         loglik = res2$logL,
         dloglik = res2$dlogL,
         cov = vcov,
         info = info,
         cvg = Cvg)
    
}

##=======================================================================
## AUTHOR: Yves Deville
##
## Find ML estimate of a two parameter Lomax distribution WITH KNOWN
## SHAPE.
##
## NB. When the scale 'beta' is known, log(1 + x / beta) is exponential
## with rate 'alpha'.
##
## TODO: add the 'scaleData' argument, find bounds for the optim.
##
##=======================================================================

`flomax1` <- function(x,
                      shape = 4.0,
                      info.observed = TRUE,
                      scaleData = TRUE,
                      cov = TRUE,
                      plot = FALSE) {

    Cvg <- TRUE
  
    if (any(x <= 0)) stop("all elements in 'x' must be >0")  
    
    n <- length(x)
    M1 <- mean(x)

    if (scaleData) {
        x <- x / M1
        CV <- sqrt(1.0 - 1.0 / n) * sd(x)
        cLogLik <- - n * log(M1)
        trans <- 1.0 / M1
    } else {
        CV <- sqrt(1.0 - 1.0 / n) * sd(x) / M1
    }
    N1 <- mean(1.0 / x)
    
    alpha <- unname(shape)
    alpha1 <- shape + 1.0
    
    ## could be used to find a zero
    dlogL1 <- function (beta) {
        xmod <- x / beta
        R <- mean(log(1.0 + xmod))
        S1 <- mean(xmod / (1.0 + xmod))
        n * (-1.0  + alpha1 * S1) / beta 
    }
    
    logL1 <- function (beta) {
        xmod <- x / beta
        R <- mean(log(1.0 + xmod))
        n * (log(alpha) - log(beta) - alpha1 * R)
    }

    if (plot) cov <- TRUE
    
    if (cov) {
        log2L1 <- function (beta) {
            xmod <- x / beta
            xmod1 <- 1.0 + xmod
            R <- mean(log(xmod1))
            S1 <- mean(xmod / xmod1)
            alpha1 <- 1.0 + alpha
            
            logL <- n * (log(alpha) - log(beta)  - alpha1 * R)
            dlogL <- n * (-1.0 + alpha1 * S1) / beta 
            d2logL <- n * (1.0  - alpha1 * S1  - alpha1 * mean(xmod / xmod1 / xmod1)) /
                beta / beta
            
            if (scaleData) {
                logL <- logL + cLogLik
                dlogL <- trans * dlogL 
                d2logL <- trans * d2logL * trans
            }
            
            list(logL = logL,
                 dlogL = dlogL,
                 d2logL = d2logL)
        }
    }
        
    if (scaleData) interv <- c((1.0 - 1.0 / alpha1) / N1, alpha1)
    else interv <- c((1.0 - 1.0 / alpha1) / N1, M1 * alpha1)
    
    checks <- unlist(sapply(interv, dlogL1))
    if ((checks[1] < 0) || (checks[2] > 0)) {
        warning("no interval found to maximise loglik")
        Cvg <- FALSE
    }
    
    res <- optimize(f = logL1, interval = interv, maximum = TRUE,
                    tol = .Machine$double.eps^0.3)
    beta.hatS <- res$maximum
    alpha.hat <- alpha
    
    if (scaleData) beta.hat <- M1 * beta.hatS
    else beta.hat <- beta.hatS

    if (!cov) {
        loglik <- res$objective
        if (scaleData) loglik <- loglik + cLogLik
        return(list(estimate = c(scale = beta.hat),
                    CV = CV,
                    loglik = loglik,
                    cvg = Cvg))
    }
    
    res2 <- log2L1(beta = beta.hatS)
    loglik <- res2$logL
        
    if (info.observed) {
        info <- -res2$d2logL
    } else {
        info <- n * alpha.hat / (alpha.hat + 2.0) / beta.hat / beta.hat
    }
    
    cov <- 1.0 / info
    sds <- sqrt(cov)
    
    if (plot) {
        Stext <- ifelse(scaleData, "(scaled data)", "")
        betas <- seq(from = interv[1], to = interv[2], length = 200)
        fs <- sapply(betas, logL1)
        dfs <- c(NA, diff(fs) / diff(betas))
        ## ind <- 1L:length(betas)
        ind <- (dfs < 20)
        plot(betas[ind], dfs[ind],
             type = "l", lty = "dotted", 
             main = sprintf("'Lomax' log-lik derivative %s", Stext),
             xlab = "beta (scale)", ylab = "dlogL")
        lines(betas[ind], sapply(betas[ind], dlogL1), col = "orangered")
        abline(h = 0, v = beta.hatS, col = "SpringGreen3")
    }
    
    list(estimate = c(scale = beta.hat),
         sd = sds,
         CV = CV,
         loglik = loglik,
         dloglik = res2$dlogL,
         cov = cov,
         info = info,
         cvg = Cvg)
    
}
