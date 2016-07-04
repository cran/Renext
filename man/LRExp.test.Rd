\name{LRExp.test}
\alias{LRExp.test}
\title{Likelihood Ratio test of exponentiality vs. GPD}
\usage{
   LRExp.test(x,
              alternative = c("lomax", "GPD", "gpd", "maxlo"),
              method = c("num", "sim", "asymp"),
              nSamp = 15000,
              simW = FALSE)
}
\arguments{
  \item{x}{

    Numeric vector of positive sample values. For the POT context this
    should be the vector of excesses over the threshold.

  }
  \item{alternative}{

    Character string describing the alternative distribution.

  }
  \item{method}{

    Method used to compute the \eqn{p}-value.

  }
  \item{nSamp}{

    Number of samples for a simulation, if \code{method} is
    \code{"sim"}.

  }
  \item{simW}{

    Logical. If this is set to \code{TRUE} and \code{method}
    is \code{"sim"}, the simulated values are returned as
    an element \code{W} in the list.

  }
}
\value{
  
  A list of results with elements \code{statistic}, \code{p.value}
  and \code{method}. Other elements are
  \item{alternative}{

    Character describing the alternative hypothesis.

  }
  \item{W}{

    If \code{simW} is \code{TRUE} and \code{method} is \code{"sim"}
    only.  A vector of \code{nSamp} simulated values of the statistic
    \eqn{W := -2 \log \textrm{LR}}{W = -2 log(LR)}.
    
  }

}
\description{
  
  Likelihood Ratio test of exponentiality vs. GPD.

}
\details{

  The Lomax and maxlo alternatives correspond to a GPD alternative with
  positive shape parameter \eqn{\xi > 0} (Lomax) and GPD with
  \eqn{\xi < 0} (maxlo).
  
  The \emph{asymptotic} distribution of the Likelihood-ratio statistic
  is known. For the GPD alternative, this is a chi-square distribution
  with one df.  For the Lomax alternative, this is the distribution of a
  product \eqn{BC} where \eqn{B} and \eqn{C} are two independent random
  variables following a Bernoulli distribution with probability
  parameter \eqn{p = 0.5} and a chi-square distribution with one df.
  
  \itemize{
  
    \item When \code{method} is \code{"num"}, a numerical
    approximation of the distribution is used. This method
    is not unlike that used by Kozubowski et al., but a different
    approximation is used. However, if \code{x} has a length
    \eqn{n > 500}, the method is turned to \code{"asymp"}.
    
    \item When \code{method} is \code{"sim"}, \code{nSamp} samples of the
    exponential distribution with the same size as \code{x} are drawn
    and the LR statistic is computed for each sample. The \eqn{p}-value
    is simply the estimated probability that a simulated LR is greater
    than the observed LR.

    \item Finally when \code{method} is \code{"asymp"}, the asymptotic
    distribution is used.

  } 
}
\note{
  
  For the Lomax alternative, the distribution of the test
  statistic has \emph{mixed type}: it can take any positive value as
  well as the value \eqn{0} with a positive probability mass. The
  probability mass is the probability that the ML estimate of the GPD
  shape parameter is negative, and a good approximation of it is
  provided by the \code{\link{pGreenwood1}} function. Note that this
  probability converges to its limit \eqn{0.5} \emph{very slowly}, which
  suggests that the asymptotic distribution provides poor results for
  medium sample sizes, say \eqn{< 100}.

  Similarly for a maxlo alternative, the distribution of the test
  statistic has mixed type: it can take any positive value as
  well as the value \eqn{0} with a positive probability mass
  approximately given by \code{1 -pGreenwood1(n)} where \eqn{n}
  is the sample size.

}
\author{
  Yves Deville
}
\seealso{
  
  \code{\link{Lomax}}, \code{\link{Maxlo}}, \code{\link{GPD}} for the
  alternatives used here.

}
\references{

  T.J. Kozubowski, A. K. Panorska, F. Qeadan, A. Gershunov and
  D. Rominger (2009)
  "Testing Exponentiality Versus Pareto Distribution via Likelihood Ratio"
  \emph{Comm. Statist. Simulation Comput.} 38(1),
  pp. 118-139.

  The approximation method used is described in the \emph{Renext
    Computing Details} report.
  
  
}
\examples{
set.seed(1234)
x <- rGPD(n = 50, loc = 0, scale = 1, shape = 0.1)
LRExp.test(x, method = "num")$p.value
LRExp.test(x, method = "asymp")$p.value
\dontrun{
## requires much time
LRExp.test(x, method = "asymp")$p.value
}
}
