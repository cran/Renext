\name{pGreenwood1}
\alias{pGreenwood1}
\title{Probability that the Greenwood's statistic is smaller than one}
\usage{ 
   pGreenwood1(n)
}
\arguments{
  \item{n}{Sample size.}
}
\value{
  
  Probability that the Greenwood's statistic is smaller than one. For a
  random sample of an exponential distribution with size \eqn{n}, this
  is the probability that the coefficient of variation is less than one,
  or the probability that the ML estimate of the GPD shape parameter
  \eqn{\xi} is negative.

}
\description{

  Probability that the Greenwood's statistic is smaller than one.

}
\details{

  The probability was computed by using the approximation of the
  quantile function of the Greenwood's statistic returned by
  \code{\link{qStat}}. The result is found by interpolating the
  distribution function for \eqn{x = 1}.

}
\examples{
n <- 8:500
plot(n, pGreenwood1(n), type = "l", col = "orangered", lwd = 2,
     log ="x", ylim =c(0.5, 0.7), main = "slow convergence to 0.5")
grid() ; abline(h = 0.5, col = "SpringGreen")
}
\author{
Yves Deville
}
