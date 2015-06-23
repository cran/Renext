\name{LRGumbel}
\alias{LRGumbel}
\title{Likelihood Ratio statistic for Gumbel vs. GEV}
\usage{

   LRGumbel(x, alternative = c("frechet", "GEV"))

}
\arguments{
  \item{x}{

    Numeric vector of sample values.

  }
  \item{alternative}{

    Character string describing the alternative.

  }
}
\value{

  The LR statistic value.

}
\description{

  Likelihood Ratio statistic for the Gumbel distribution vs. GEV.

}
\details{
  
  The Likelihood-Ratio statistic is actually \eqn{W:=-2 \log
    \textrm{LR}}{-2 log LR} where LR is the ratio of the likelihoods
    \emph{Gumbel} to \emph{alternative distribution}.

}
\note{
  
  When the alternative is \code{"frechet"}, the statistic has a
  distribution of mixed type under the null hypothesis of a Gumbel
  distribution.

}
\author{
  Yves Deville
}
\seealso{

  \code{\link{LRGumbel.test}} for the related LR test of Gumbelity.
  
}