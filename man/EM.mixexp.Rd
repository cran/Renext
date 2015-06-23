\name{EM.mixexp}
\alias{EM.mixexp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

  Expectation-Maximisation for a mixture of exponential distributions

}
\description{
  
  Experimental function for Expectation-Maximisation (EM) estimation

}
\usage{
EM.mixexp(x, m = 2)
}
\arguments{
  \item{x}{

    Sample vector with values \code{>0}.

  }
  \item{m}{

    Number of mixture components.

  }
}
\details{
  
   The EM algorithm is very simple for exponential mixtures (as well as
   for many other mixture models).

   According to a general feature of EM, this iterative method leads to
   successive estimates with increasing likelihood but which may
   converge to a local maximum of the likelihood.

 }
\value{
  List with
  \item{estimate}{

    Estimated values as a named vector.

  }
  \item{logL}{

    Vector giving the log-likelihood for successive iterations.

  }
  \item{Alpha}{

    Matrix with \code{m} columns giving probability weights for
    successive iterations. Row with number \code{it} contains the
    \code{m} probabilities at iteration \code{it}.

  }

  \item{Theta}{

    Matrix with \code{m} columns giving the estimates of the
    \code{m} expectations for the successive iterations

  }
}
\author{

  Yves Deville

}
\note{
  
  The estimation is done for expectation (inverse rates) but the
  \code{estimate} vector in the result contains rates for compatibility
  reasons (e.g with exponential).

}
\seealso{

  \code{\link{mom.mixexp2}} and \code{\link{ini.mixexp2}} for "cheap"
  estimators when \code{m = 2}.

}
\examples{
set.seed(1234)
x <- rmixexp2(n = 100, prob1 = 0.5, rate2 = 4)
EM.mixexp(x) -> res
res$estimate
matplot(res$Theta, type = "l", lwd = 2,
        xlab = "iteration", ylab = "theta",
        main = "exponential inverse rates")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
