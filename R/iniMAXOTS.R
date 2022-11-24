##=======================================================================
## auxiliary variables for a gumbel distribution
## see "Renext Computing Details".
##
##=======================================================================

gumbelAux <- function(r) {
    res1 <- rep(0, length(r))
    res2 <- rep(pi * pi / 6, length(r))
    for (i in 1:length(r)) {
        if (r[i] > 1) {
            temp <- 1 / (1:(r[i]-1)) 
            res1[i] <- sum(temp)
            res2[i] <- res2[i] - sum(temp^2) 
        }
    }
    list(r = r, nu = res1, kappa2 = res2)
}

## ****************************************************************************
##' Initial values for ML estimation from MAX data.
##'
##' Cheap estimates of \code{lambda} and \code{sigma} for the
##' aggregated POT model with exponential exceedances are obtained
##' using MAX historical data only. One possibility is to use the
##' minimal value in each MAX block as a response for a linear
##' regression on the log-block duration. The noise variance is then
##' proportional to the slope 'beta' which needs some care. These
##' regression estimates allows the estimate the two parameters
##' \code{lambda} and \code{sigma} by a simple transformation. Another
##' possible estimate of \code{sigma} comes from the spacings in
##' relation with Renyi's representation. Here we combine the two
##' estimates of \code{sigma} in order to get an improved estimate and
##' then \code{lambda} is computed.
##' 
##' @title Initial Values from MAX Data
##'
##' @param MAX A list as built by using the \code{\link{makeMAXdata}}
##' function.
##'
##' @param threshold The threshold to be used.
##' 
##' @param distname.y Name of the distribution for the excesses, which
##' must fall in the GPD family.
##' 
##' @return A vector of estimated parameters.
##'
##' @references See the \emph{Renext Computing Details} technical
##' report (Chap. 3) for the details.
##' 
##' @examples
##' set.seed(123)
##' lambda <- 1; scale <- 100; shape <- 0.00
##' nBlocks <- 1L + rpois(1, lambda = 10)
##' w <- rgamma(nBlocks, shape = 3, scale = 2)
##' n <- 1L + rpois(nBlocks, lambda = lambda * w)
##' r <- 1L + rpois(nBlocks, lambda = 5)
##' r <- pmin(n, r)
##' block <- rep(1L:nBlocks, n)
##' y <- rGPD(sum(n), scale = scale, shape = shape)
##' y <- tapply(y, block, identity)
##' MAXdata <- list()
##' for (ib in 1L:nBlocks) MAXdata[[ib]] <- sort(y[[ib]], decreasing = TRUE)[1:r[ib]]
##' MAX <- list("effDuration" = w, "data" = MAXdata, "r" = r)
##' parms <- parIni.MAX(MAX, threshold = 0, dist = "GPD")
##' 
parIni.MAX <- function(MAX, threshold, distname.y = "exp") {

    if (!(distname.y %in% c("exp", "exponential", "gpd", "GPD"))) {
        stop("unauthorised value for 'distname.y'")
    }
    r <- MAX$r
    nBloc <- length(r)

    ## if (nBloc < 2L) stop("The number of blocks must be >= 2") 
    
    ## minimal r-largest statistic
    y <- unlist(lapply(MAX$data, min))
    Spac <- unlist(lapply(MAX$data, spacings))
    nSpac <- length(Spac) 
    aux <- gumbelAux(r)
    
    ## compute sigmaHat
    sigmaSpac <- mean(Spac)
    
    ## regression estimates, see "Renext Computing Details"
    x <- log(MAX$effDuration) - aux$nu
    wt <-  1.0 / aux$kappa2
    covmat <- cov.wt(cbind(x, y), wt = wt, method = "ML")$cov
    
    a <- 1.0 / mean(wt)
    betaReg <- polyroot(c(-covmat[2L, 2L], covmat[1L, 2L], a))
    betaReg <- sort(Re(betaReg))[2L]
    
    if (nSpac) {
        if (nBloc > 1L) {
            den <- nSpac + nBloc / 2.0
            wSpac <- nSpac / den
            sigmaHat <- wSpac * sigmaSpac  + (1 - wSpac) * betaReg
        } else {
            sigmaHat <- sigmaSpac
        }
    } else {
        sigmaHat <- betaReg
    }
    
    alphaReg <- weighted.mean(y, w = wt) - betaReg * weighted.mean(x, w = wt)
    lambdaHat <- exp((alphaReg - threshold) / sigmaHat + digamma(1))
    
    if (distname.y %in% c("exp", "exponential")) {
        return(c("lambda" = lambdaHat, "rate" = 1.0 / sigmaHat))
        ## "alphaReg" = alphaReg, "betaReg" = betaReg, "sigmaSpac" = sigmaSpac)
    } else if (distname.y %in% c("gpd", "GPD")) {
        return(c("lambda" = lambdaHat, "scale" = sigmaHat, "shape" = 0.0))
    }   
    
}

## ****************************************************************************
##' Initial values for ML estimation from OTS data.
##'
##' Cheap estimates of \code{lambda} and \code{sigma} for the
##' aggregated POT model with exponential exceedances are obtained
##' using OTS historical data only. An estimate of \code{lambda} can
##' be obtained by using a Poisson regression with the block duration
##' as covariate. An offset must be used and the canonical link is
##' used. While the estimated slope could be used to find an estimate
##' of \code{sigma}, an estimated based on the spacings is used
##' instead.
##' 
##' @title Initial Values from OTS Data
##'
##' @param OTS A list as built by using the \code{\link{makeOTSdata}}
##' function.
##'
##' @param threshold The threshold to be used.
##' 
##' @param distname.y Name of the distribution for the excesses, which
##' must fall in the GPD family.
##'
##' @return A vector of parameters
##' 
##' @references See the \emph{Renext Computing Details} technical
##' report (Chap. 3) for the details.
##'
##' @examples
##' set.seed(123)
##' lambda <- 1; scale <- 100
##' nBlocks <- 1L + rpois(1, lambda = 10)
##' w <- rgamma(nBlocks, shape = 3, scale = 2)
##' threshold <- rexp(nBlocks, rate = 1 / scale)
##' OTSdata <- list()
##' r <- integer(nBlocks)
##' for (ib in 1L:nBlocks) {
##'    S <- pexp(q = threshold[ib], rate = 1. / scale, lower.tail = FALSE)
##'    r[ib] <- rpois(1, lambda = lambda * S)
##'    OTSdata[[ib]] <- threshold[ib] + rexp(r[ib], rate = 1. / scale)
##' }
##' OTS <- list("effDuration" = w, "data" = OTSdata, "r" = r,
##'             "threshold" = threshold)
##' parms <- parIni.OTS(OTS, threshold = 0, dist = "GPD")
##' 
parIni.OTS <- function(OTS, threshold, distname.y = "exp") {

    if (!(distname.y %in% c("exp", "exponential", "gpd", "GPD"))) {
        stop("unauthorised value for 'distname.y'")
    }
    
    dThreshold <- OTS$threshold - threshold
    
    if (any(dThreshold <= 0)) {
        stop("all values in 'OTSthreshold' must be > 'threshold'")
    }
    
    ## estimate the rate using a GLM, see "Renext Computing Details"
    R <-  OTS$r 
    lw <- log(OTS$effDuration)
    fit <- glm(R ~ dThreshold, family = poisson, offset = lw)
    betaReg <- as.numeric(coef(fit))
    
    lambdaReg <- exp(betaReg[1])
    ## sigmaReg <- - 1.0 / betaReg[2] ## not used 
    
    Ybar <- unlist(lapply(OTS$data, mean)) - unlist(OTS$threshold)
    sigmaHat <- weighted.mean(Ybar, w = OTS$r)
    
    if (distname.y %in% c("exp", "exponential")) {
        return(c("lambda" = lambdaReg, "rate" = 1.0 / sigmaHat)) 
    } else if (distname.y %in% c("gpd", "GPD")) {
        return(c("lambda" = lambdaReg, "scale" = sigmaHat, "shape" = 0.0))
    }
    
}


