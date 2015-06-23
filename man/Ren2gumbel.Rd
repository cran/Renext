\name{Ren2gumbel}
\alias{Ren2gumbel}

\title{
  
  Translate a vector of coefficients from a Renewal-POT model with
  exponential excesses to a vector of Gumbel parameters
  
}
\description{

  Translate a vector of coefficients from a Renewal-POT model with
  exponential excesses to a vector of Gumbel parameters.

}
\usage{
Ren2gumbel(object,
           threshold = NULL,
           w = 1,
           distname.y = c("exponential", "exp"),
           jacobian = (length(w) == 1L),
           vcovRen = NULL)
}
\arguments{
  \item{object}{

    A named vector of parameters or an object of class \code{"Renouv"}.
    In the first case, the names of the vector element must conform to
    the exponential distribution so the vector must be of length 2 with
    names \code{"lambda"} and \code{"rate"}.

  }
  \item{threshold}{
    
    A threshold associated with the parameters. If \code{object} is an
    object with class \code{"Renouv"}, its threshold slot will be used.

  }
  \item{w}{
    
    A block duration or a vector of block durations.

  }
  \item{distname.y}{

    The name of the distribution for the excesses.  Can be either
    \code{"exponential"} or \code{"exp"}. The choice has no impact on
    the computations, but this name will be attached to the result as an
    attribute and may affect later use.

  }
  \item{jacobian}{
    
    Logical. If \code{TRUE} the jacobian matrix of the transformation
    will be computed and attached to the result as an attribute.

  }
  \item{vcovRen}{
    
    A covariance matrix for the Renouv parameters. 

  }
}
%%\details{}
\value{
  
  A vector of GEV parameters if \code{w} has length 1, and a matrix if
  \code{w} has length \code{> 1}. The returned objects has attributes.

}
%%\references{}
\author{
  
  Yves Deville

}
%%\note{}
\seealso{
  
  \code{\link{Ren2gev}} for the translation of Renouv parameters
  corresponding to GPD excesses.

}
\examples{
## Fit a Renouv model with exponential excesses (default)
fit <- Renouv(Garonne)
## Convert to gumbel (usable for one-year block maxima)
parGumbel <- Ren2gumbel(fit)
## Retrieve the 'Renouv' model by giving the right threshold
parRen <- gumbel2Ren(parGumbel,
                     threshold = 2500,
                     vcovGumbel = attr(parGumbel, "vcov"),
                     plot = TRUE)
## Build a compatible model under the assumption of one event by
## year
parRen2 <- gumbel2Ren(parGumbel,
                      lambda = 1.00,
                      vcovGumbel = attr(parGumbel, "vcov"),
                      plot = TRUE)
parRenNames <- c("lambda", "rate")
## Build a 'Renouv' object without estimation
myVcov <- attr(parRen, "vcov")[parRenNames, parRenNames]
fitNew <- RenouvNoEst(threshold = attr(parRen, "threshold"),
                      estimate = parRen,
                      distname.y = "exp",
                      cov = myVcov)
## Compare return levels
cbind(roundPred(fit$pred)[ , -2], roundPred(fitNew$pred)[ , -2])
## idem for the putative 'Renouv' with rate 1
myVcov2 <- attr(parRen2, "vcov")[parRenNames, parRenNames]
fitNew2 <- RenouvNoEst(threshold = attr(parRen2, "threshold"),
                       estimate = parRen2,
                       distname.y = "exp",
                       cov = myVcov2)
cbind(roundPred(fit$pred)[ , -2], roundPred(fitNew2$pred)[ , -2])
}

%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
