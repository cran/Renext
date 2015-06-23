\name{LRGumbel.test}
\alias{LRGumbel.test}
\encoding{UTF8}
\title{Likelihood Ratio test for the Gumbel distribution}
\usage{
   LRGumbel.test(x,
                 alternative = c("frechet", "GEV"),
                 method = c("num", "sim", "asymp"),
                 nSamp = 1500,
                 simW = FALSE)
}
\arguments{
  \item{x}{
    
    Numeric vector of sample values.

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

  Likelihood Ratio test of Gumbel vs. GEV

}
\details{

  The asymptotic distribution of the Likelihood-ratio statistic is
  known. For the GEV alternative, this is a chi-square distribution with
  one df.  For the Fréchet alternative, this is the distribution of a
  product \eqn{XY} where \eqn{X} and \eqn{Y} are two independent random
  variables following a Bernoulli distribution with probability
  parameter \eqn{p = 0.5} and a chi-square distribution with one df.

  \itemize{
    
    \item When \code{method} is \code{"num"}, a numerical approximation
    of the distribution is used.

    \item When \code{method} is \code{"sim"}, \code{nSamp} samples of
    the Gumbel distribution with the same size as \code{x} are drawn and
    the LR statistic is computed for each sample. The \eqn{p}-value is
    simply the estimated probability that a simulated LR is greater than
    the observed LR. This method requires more computation time than
    the tow others.
    
    \item Finally when \code{method} is \code{"asymp"}, the asymptotic
    distribution is used.
    
  }
  
}
\note{

  For the Fréchet alternative, the distribution of the test statistic
  has \emph{mixed type}: it can take any positive value as well as the
  value \eqn{0} with a positive probability mass. The probability mass
  is the probability that the ML estimate of the GEV shape parameter is
  negative.
  
  When \code{method} is \code{"sim"}, the computation can be slow
  because each of the \code{nSamp} simulated values requires two
  optimisations. The \code{"asymp"} method provides an acceptable
  precision for \eqn{n \geq 50}{n >= 50}, and may even be used for
  \eqn{n \geq 30}{n >= 30}.

}
\examples{
set.seed(1234)
x <- rgumbel(60)
res <- LRGumbel.test(x)
}
\author{

  Yves Deville

}
