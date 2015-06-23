\name{Hpoints}
\alias{Hpoints}
\encoding{UTF-8}
\title{Plotting positions for exponential return levels}
\usage{
  Hpoints(n)
}
\arguments{
  \item{n}{Sample size.}
}
\value{

  Numeric vector of plotting positions with length \code{n}.

}
\description{
  
  Plotting positions for exponential return level plots.

}
\details{
  
  The plotting positions are numeric values to use as the abscissae
  corresponding to the order statistics in an exponential return level
  plot. They range from 1 to about \eqn{\log n}{log(n)}. They can be
  related to the plotting positions given by \code{\link{ppoints}}.

  The returned vector \eqn{\mathbf{H}}{H} has elements \deqn{H_{i} =
  \frac{1}{n} + \frac{1}{n-1} + \dots + \frac{1}{n + 1 -i}}{ H[i] = 1 /
  n + 1 / (n + 1) + ... + 1 / (n + 1 - i)} for \eqn{1 \leq i \leq n}{1
  <= i <= n}. This is the expectation of the \eqn{i}-th order statistic
  for a sample of the standard exponential distribution, see
  e.g. chap. 4 of Embrechts et al.  } \references{
  
  Embrechts P., Klüppelberg C. and Mikosch T. (1997) \emph{Modelling
   Extremal Events for Insurance and Finance}. Springer.

}
\author{
  Yves Deville
}
\note{
  
  For \eqn{n} large enough, the largest value \eqn{H_n}{H[n]} is
  approximately \eqn{\gamma + \log n}{\gamma + log(n)} where
  \eqn{\gamma} is the Euler-Mascheroni constant, and \eqn{\exp
  H_n}{exp(H[n])} is about \eqn{1.78 n}{1.78 * n}.  Thus if the Hpoints
  are used as plotting positions on a return level plot, the largest
  observation has a return period of about \eqn{1.78 n}{1.78 * n} years.

}
\seealso{
  \code{\link{ppoints}}.
}
\examples{
n <- 30
set.seed(1234)
x <- rGPD(n, shape = 0.2)
plot(exp(Hpoints(n)), sort(x), log = "x",
     main = "Basic return level plot")

}