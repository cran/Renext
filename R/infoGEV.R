## DO NOT ROXYGENIZE!!!
##
##' Expected Information matrix for the GEV distribution.
##'
##' This matrix is also computed in the function \code{egevd}
##' of the \pkg{EnvStats} package, but not returned.
##' 
##' @title Expected Information matrix for the GEV distribution
##'
##' @param param Named numeric vector of parameters.
##'
##' @param n Sample size.
##' 
##' @return The Fisher information matrix.
##'
##' @references
##'
##' Kotz S. and Nadarajah S.(2000)
##' \emph{Extreme Value Distributions Theory and Applications}.
##' Imperial College Press.
##'
##' Note that there is a typo in the 2000 edition of the book p. 63
##' for the \code{["scale", "scale"]} element.
##'  
##' 
infoGEV <- function(param, n) {

    pn <- c("loc", "scale", "shape")

    if (!setequal(x = names(param), y = pn)) {
        stop("'param' must be a named numeric vector ",
             "with names", pn)
    }

    mu <- param["loc"]
    sigma <- param["scale"]
    xi <- param["shape"]
    xi1 <- 1.0 + xi
    gamxi2 <- gamma(2.0 + xi)
    Euler <- -digamma(1)
    aux <- 1.0 - Euler + 1.0 / xi
    
    p <- xi1^2 * gamma(1 + 2 * xi)
    q <- gamxi2 * (digamma(xi1) + xi1 / xi)

    info <- array(NA, dim = c(3L, 3L), dimnames = list(pn, pn))
    
    info["loc", "loc"] <- p / sigma^2

    info["scale", "scale"] <- (1 - 2 * gamxi2 + p) / (sigma * xi)^2 
    
    info["shape", "shape"] <- (pi^2 / 6 + aux^2 -
                                   2.0 * q / xi + p / xi^2) / xi^2  

    info["loc", "scale"] <- info["scale", "loc"] <-
        -(p - gamxi2) / xi / sigma^2
    
    info["loc", "shape"] <- info["shape", "loc"] <-
        (q - p / xi) / sigma / xi
    
    info["scale", "shape"] <- info["shape", "scale"] <-
        (aux - gamxi2 / xi - q + p / xi) / sigma / xi^2
    
    info <- n * info
    info
    
    
}
