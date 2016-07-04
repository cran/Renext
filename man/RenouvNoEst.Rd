\name{RenouvNoEst}
\alias{RenouvNoEst}
\title{
Define a 'renouvellement' model without estimation
}
\description{
  Build a 'renouvellement' model using parameters
  given by the user.
}
\usage{
RenouvNoEst(threshold,
            estimate = NULL,
            distname.y = "exponential",
            fixed.par.y = NULL,
            trans.y = NULL,
            pct.conf = c(95, 70),
            rl.prob = NULL,
            prob.max = 1 - 1e-04,
            pred.period = NULL,
            cov = NULL,
            nb.OT = NULL,
            infer.method = NULL)
}
\arguments{
  \item{threshold}{

    The threshold.

  }
  \item{estimate}{

    Numeric named vector containing the estimates for the parameters. It
    must be compatible with the distribution chosen, and must contain in
    first position an element named \code{"lambda"} representing an
    estimated event rate in \emph{events by year}.

  }
  \item{distname.y}{

    Character giving the name of the distribution.

  }
  \item{fixed.par.y}{

    Numeric named vector containing values for vectors which are
    considered as fixed (and not estimated).

  }
  \item{trans.y}{

    Transformation as in \code{\link{Renouv}}. Used only when
    \code{distname.y} is equal to \code{"exponential"}.

  }
  \item{pct.conf}{

    Vector of percents for confidence limits.

  }
  \item{rl.prob}{

    Probability used in the return level computations. These values are
    used for instance in return level plots produced with the
    \code{\link{plot.Renouv}} method. When \code{NULL} a default vector
    is used.

  }
  \item{prob.max}{

    Maximal probability for which computations are done.

  }
  \item{pred.period}{

   Vector of periods for which predicted return levels will be computed.

  }
  \item{cov}{

    Covariance matrix for the provided estimated parameters.  Must have
    rownames and colnames in accordance with those of
    \code{estimate}. This covariance matrix is used to build confidence
    limits on parameters and on return levels using the \emph{delta
    method}.

  }
  \item{nb.OT}{

    Number of data over the threshold used for estimation.  This will be
    used only when \code{distname.y} is equal to \code{"exponential"}.

  }
  \item{infer.method}{

  Inference method. Will normally be the \emph{delta method}.

  }
}
\details{

  This function is used for plotting or comparing models with known
  parameter estimates but with no data available.
  
  The parameters estimates should be accompanied with a covariance
  matrix assuming an approximately normal joint distribution of
  these. This matrix is usually obtained by computing the numerical
  derivatives of the log-likelihood at the second order at the
  estimates. This covariance is used to compute approximate confidence
  limits for the return levels of the unknown true distribution that was
  estimated.
  
}
\value{

  An object of class \code{"Renouv"} representing a 'renouvellement'
  model similar to those built with \code{\link{Renouv}}. This is mainly
  a list. Note however that some list elements found in \code{Renouv}
  objects built by \code{Renouv} can not be found here.  For instance,
  the returned objects embeds no goodness-of-fit results since the
  object is created without making use of any data.

}
%%\references{}
\author{
  
  Yves Deville

}
\seealso{

  \code{\link{Renouv}} to estimate such models.

}
\examples{
##======================================================================
## Example from S. Coles' book, page 86 'rainfall data'.
## Note that the first parameter is here the rate 'lambda', and no the
## probability of exceedance as in Coles' book.
##======================================================================
estimate <- c(lambda = 152 / 48, scale = 7.44, shape = 0.184)          
cov <- matrix(c(4.9e-7 * (17531 / 48)^2,  0.0000,  0.0000,
                0.0000,  0.9180, -0.0655,
                0.0000, -0.0655,  0.0102),
              nrow = 3)
colnames(cov) <- rownames(cov) <- names(estimate)
renNE <- RenouvNoEst(threshold = 30, distname.y = "gpd",
                     pct.conf = c(95, 70),
                     estimate = estimate,
                     nb.OT = 152, cov = cov)
summary(renNE)
plot(renNE, main = "Daily rainfall data SW England", ylim = c(0, 400))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
