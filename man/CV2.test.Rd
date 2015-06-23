\name{CV2.test}
\alias{CV2.test}
\title{CV2 test of exponentiality}
\usage{
   CV2.test(x, method = c("num", "sim", "asymp"), nSamp = 15000)
}
\arguments{
  \item{x}{

    Numeric vector giving the sample.

  }
  \item{method}{

    Method used to compute the \eqn{p}-value.  Can be \code{"asymp"},
    \code{"num"} or \code{"sim"} as in \code{\link{LRExp.test}}.

  }
  \item{nSamp}{

    Number of samples used to compute the \eqn{p}-value when
    \code{method} is \code{"sim"}.

  }
}
\value{
  
  A list of test results.

  \item{statistic, p.value}{
    
    The test statistic, i.e. the squared coefficient of
    variation \eqn{\textrm{CV}^2}{CV^2} and the \eqn{p}-value.
    
  }
  \item{df}{
    
    The sample size.
    
  }
  \item{method}{
    
    Description of the test method.
    
  }

}
\description{
  
  Test of exponentiality based on the squared coefficient of variation.

}
\note{

  This test is sometimes referred to as \emph{Wilk's exponentiality
  test} or as \emph{WE1 test}.  It works quite well for a Lomax
  alternative (i.e.  GPD with shape \eqn{\xi >0}), and hence can be
  compared to Jackson's test and the Likelihood-Ratio (LR) test of
  exponentiality.  However, this test has lower power that of the two
  others while having a comparable computation cost due to the
  evaluation of the Greenwood's statistic distribution.
  
}
\details{
  
  The distribution of \eqn{\textrm{CV}^2}{CV^2} is that of
  \emph{Greenwood's statistic} up to normalising constants.  It
  approximately normal with expectation \eqn{1} and standard deviation
  \eqn{2/\sqrt{n}}{2/sqrt(n)} for a large sample size \code{n}. Yet the
  convergence to the normal is known to be \emph{very slow}.
  
}

\references{

  S. Ascher (1990) "A Survey of Tests for Exponentiality"
  \emph{Commun. Statist. Theory Methods}, 19(5), pp. 1811-1525. 

}
\examples{
n <- 30; nSamp <- 500
X <- matrix(rexp(n * nSamp), nrow = nSamp, ncol = n)
pVals <- apply(X, 1, function(x) CV2.test(x)$p.value)
plot(pVals)  ## should be uniform on (0, 1)
}
\author{
Yves Deville
}
\seealso{
  
  The function \code{\link{CV2}} that computes the statistic and
  \code{\link{LRExp.test}} or \code{\link{Jackson.test}} for functions
  implementing comparable tests or exponentiality with the same
  arguments.

}
