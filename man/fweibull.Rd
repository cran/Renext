\name{fweibull}
\alias{fweibull}
\title{

  ML estimation of classical Weibull distribution

}
\description{
  
  Fast Maximum Likelihood estimation of the classical two parameters
  Weibull distribution.

}
\usage{

fweibull(x, info.observed = TRUE, scaleData = TRUE, cov = TRUE,
         check.loglik = FALSE)

}
\arguments{

  \item{x}{

    Sample vector to be fitted. Should contain only positive
    non-NA values.

  }
  \item{info.observed}{

    Should the observed information matrix be used or the expected one
    be used?

  }
  \item{scaleData}{

    Should the data be scaled before estimation? If \code{TRUE},
    the observations in \code{x} (which are positive) are divided by their
    mean value.  The results are in theory not affected by this
    transformation, but scaling the data could improve the estimation in
    some cases. 
    
  }
  \item{cov}{
    
    Should the covariance of estimates be computed?
    
  }  
  \item{check.loglik}{

    If \code{TRUE}, the log-likelihood is recomputed using
    \code{dweibull} function with \code{log = TRUE}. The result is
    returned as a list element.

  }
}
\details{

  The ML estimates are obtained thanks to a reparameterisation with
  \eqn{\eta = scale^{1/shape}}{eta = scale^{1/shape}} in place of
  \code{shape}. This allows the maximisation of a one-dimensional
  likelihood \eqn{L}{L} since the \eqn{\eta}{eta} parameter can be
  concentrated out of \eqn{L}{L}. This also allows the determination of
  the \emph{expected} information matrix for
  \eqn{[shape,\,\eta]}{[\code{shape}, \code{eta}]} rather than the usual
  \emph{observed} information.

}
\value{
  A list
  \item{estimate}{

    Parameter ML estimates.

  }
  \item{sd}{

    The (asymptotic) standard deviation for estimate.

  }
  \item{cov}{

    The (asymptotic) covariance matrix computed from theoretical or
    observed Information matrix.

  }
  \item{eta}{

    The estimated value for eta.

  }
}
\author{

  Yves Deville

}
\note{

  The default value of \code{info.observed} was set to \code{TRUE} from
  version \code{3.0-1} because standard deviations obtained with this
  choice are usually better.
  
  %% A well-accepted verdict is that expected information is better than
  %% observed information for problems with no missing data as it is the
  %% case here. We suspect in the present estimation context that the
  %% expected information matrix has often a better condition number than
  %% its observed version.
 
}
\seealso{

  \code{\link{weibplot}} for Weibull plots.

}
\examples{

n <- 1000
set.seed(1234)
shape <- 2 * runif(1)
x <- 100 * rweibull(n, shape = 0.8, scale = 1)
res <- fweibull(x)

## compare with MASS
if (require(MASS)) {
   res2 <- fitdistr(x , "weibull")
   est <- cbind(res$estimate, res2$estimate)
   colnames(est) <- c("Renext", "MASS")
   loglik <- c(res$loglik, res2$loglik)
   est <- rbind(est, loglik)
   est
}

## Weibull plot
weibplot(x,
         shape = c(res$estimate["shape"], res2$estimate["shape"]),
         scale = c(res$estimate["scale"], res2$estimate["scale"]),
         labels = c("Renext 'fweibull'", "MASS 'fitdistr'"),
         mono = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
