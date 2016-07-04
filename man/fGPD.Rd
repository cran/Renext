\name{fGPD}
\alias{fGPD}
\title{
  
  Fit a two-parameters Generalised Pareto Distribution from a sample

}
\description{
  
    Fit a two-parameters Generalised Pareto Distribution from a sample.

  }
\usage{

fGPD(x,
     info.observed = TRUE,
     shapeMin = -0.8,
     dCV = 1e-04,
     cov = TRUE,
     trace = 0)

}
\arguments{
  \item{x}{

    The sample data vector.

  }
  \item{info.observed}{
    
    Logical. The observed information matrix is used when \code{TRUE}
    and the expected information matrix is used when \code{FALSE}.

  }
  \item{shapeMin}{

    Lower bound on the shape parameter. This must be \code{>= -1.0}
    since otherwise the ML estimate is obtained with the \code{scale}
    parameter equal to \code{max(x)}.
    
  }
  \item{dCV}{
    
    Half-length of a small interval centered on 1.0.  When the
    Coefficient of Variation (CV) falls in this interval, an exponential
    distribution is fitted, see \bold{Details}.
    
  }
  \item{cov}{
    
     Logical. If \code{FALSE}, a minimal estimation is performed with
     no covariance matrix or derivative returned. This can be useful
     when a large number of ML estimations are required, e.g. to sample
     from a likelihood ratio.

   }
   \item{trace}{
    
     Integer level of verbosity. The value \code{0} prints nothing.

  }


   
}
\details{
  
  This function mainly relies on the \code{\link{flomax}} and
  \code{\link{fmaxlo}} functions. When \code{CV} is larger than
  \code{1.0 + dCV}, a Lomax distribution is fitted by Maximum Likelihood
  using a concentrated log-likelihood.  When instead \code{CV} is
  smaller than \code{1.0 - dCV}, a maxlo distribution is
  fitted. Finally, when \code{CV -1.0} has absolute value \code{<= dCV},
  an exponential distribution is fitted.  In all cases, the result is
  translated into a parameter vector for the GPD.

  Note that when \code{CV} is close to \code{1.0}, fitting a Lomax or a
  maxlo distribution can lead to problems because the estimated values
  of the shape and scale parameter are large, and because the
  concentrated log-likelihood is a flat function of the scale parameter.

}
\value{
  A list
  \item{estimate}{

    Vector containing the estimated values of the unknown parameters.

  }
  \item{CV}{

    The coefficient of variation of \code{x} computed using
    \code{length(x)} as denominator in the variance estimation.

  }
  \item{logLik, dlogLik}{

    The maximised value of the log-likelihood and the vector of its
    first order derivatives, which should be close to zero.

  }
  \item{sd, cov}{

    Vector of approximated standard deviations and covariance matrix
    for the estimated parameters. These are based on the inversion of
    expected information matrix.
    
  }
  \item{sd, cov}{
    
    Vector of standard deviations and covariance matrix for the
    estimates if the \code{cov} formal is \code{TRUE}.
    
  }
  \item{cvg}{

    Logical. Was convergence reached? This logical flag is set to
    \code{TRUE} and remains for compatibility reasons.

  }
  
}
\references{

  See the \emph{Renext Computing Details} document.

}
\author{
  
Yves Deville

}
\note{

  It may happen that the estimated shape parameter is \code{< -0.5}, in
  which case the expected information matrix can not be computed, nor
  does the covariance matrix and standard deviations. In this case, the
  \code{cov} and \code{sd} objects contain \code{NA} values.  This
  problem can arise also when the shape parameter is greater than but
  close to the value \code{-0.5}. Even when \code{info.observed} is
  \code{TRUE}, the information matrix, covariance and standard
  deviations are set to \code{NA}.

  When the true (unknown) value is is \code{< -0.5}, the regularity
  conditions required in the ML approximated inference do not hold.

  The default value of \code{info.observed} was set to \code{TRUE} from
  version \code{3.0-1} because standard deviations obtained with this
  choice are usually better.
  
}
\seealso{

  The \code{\link{fmaxlo}} and \code{\link{flomax}} functions.

}
\examples{
\dontrun{
set.seed(123456)
n <- 500
ns <- 1000
xi <- runif(ns, min = -0.5, max = 0.5)
X <- matrix(nrow = n, ncol = ns)

for (i in 1:length(xi)) {
  Xi <- rgpd(n, scale = 1, shape = xi[i])
  X[ , i] <- Xi
  res1 <- fGPD(Xi)
  res2 <- try(fpot(Xi, threshold = 0.0))
  if (inherits(res2, "try-error")) {
    cat(res2, "\n")
    break
  }
  logLik1 <- res1$loglik; logLik2 <- logLik(res2)
  if (abs(logLik1 - logLik2) > 0.001) {
    cat(sprintf("i = \%d, xi = \%7.4f\n", i, xi[i]))
    mat <- rbind(c(res1$estimate[1:2], logLik = logLik1),
                 c(res2$estimate[1:2], logLik = logLik2))
    rownames(mat) <- c("fGPD", "fpot")
    print(mat)
  }
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
