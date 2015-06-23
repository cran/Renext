\name{fgamma}
\alias{fgamma}
\title{ML estimation of the Gamma distribution}
\description{
  
  Fast Maximum Likelihood estimation of the Gamma distribution.

}
\usage{
fgamma(x, check.loglik = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{

    Sample vector to be fitted. Should contain only positive non-NA
    values.

  }
  \item{check.loglik}{

    If \code{TRUE}, the log-likelihood is recomputed using \code{dgamma}
    function with \code{log =} \code{TRUE}. The result is returned as a
    list element.

  }
}
\details{
  
  The likelihood is concentrated with respect to the scale
  parameter. The concentrated log-likelihood is a strictly concave
  function of the shape parameter which can easily maximised
  numerically.

}
\value{
  A list with the following elements
  \item{estimate}{
    
    Parameter ML estimates.
    
  }
  \item{sd}{

    Vector of (asymptotic) standard deviations for the estimates.

  }
  \item{loglik}{

    The maximised log-likelihood.

  }
  \item{check.loglik}{

    The checked log-likelihood.

  }
  \item{cov}{

    The (asymptotic) covariance matrix computed from theoretical or
    observed information matrix.

  }
  \item{info}{

    The information matrix.

  }
}
\author{
  
  Yves Deville

}
\note{

  The distribution is fitted by using the \code{scale} parameter rather
  than \code{rate} (inverse of \code{scale}).

}
\seealso{

  \code{\link{GammaDist}} in the \pkg{stats} package.

}
\examples{
set.seed(9876)
alpha <- 0.06
beta <- rexp(1)
n <- 30
x <- rgamma(n, shape = alpha, scale = beta)
fit <- fgamma(x, check.loglik = TRUE)

## compare with MASS results
if (require(MASS)) {
   fit.MASS <- fitdistr(x, densfun = "gamma")
   rate <- 1 / fit$estimate["scale"]
   est <- c(fit$estimate, rate = rate)
   der <- rate * rate ## derivative of rate w.r.t scale
   sdest <- c(fit$sd, rate = der * fit$sd["scale"])
   tab <- rbind(sprintf(" \%10.8f ", est),
                sprintf("(\%10.8f)", sdest))
   colnames(tab) <- c("shape", "scale", "rate")
   rownames(tab) <- c("est", "sd")
   noquote(tab)
}

}

