\name{CV2}
\alias{CV2}
\title{Squared Coefficient of Variation}
\usage{
   CV2(x)
}
\arguments{
  \item{x}{

    Numeric vector or matrix.

  }
}
\value{

  Numeric vector of the squared coefficients of variation.

}
\description{

  Squared Coefficient of Variation.

}
\details{
  
  Compute the squared Coefficient of Variation of one or several samples
  provided as a numeric vector or matrix.

}
\note{
  
  The squared coefficient of variation is the ratio
  \eqn{S^2/\bar{X}^2}{S^2/xbar^2} where \eqn{\bar{X}}{xbar} and \eqn{S^2}
  are the sample mean and the sample variance. The variance is computed
  using the sample size \eqn{n} as denominator, rather than the usual
  \eqn{n-1}.

}
\examples{
n <- 30; nSamp <- 500
X <- matrix(rexp(n * nSamp), nrow= nSamp, ncol = n)
W <- CV2(X)
plot(density(W), main = "CV2 of exponential samples")
}
