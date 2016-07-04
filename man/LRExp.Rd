\name{LRExp}
\alias{LRExp}
\title{Likelihood Ratio statistic for exponential vs. GPD}
\usage{
   LRExp(x, alternative = c("lomax", "GPD", "gpd", "maxlo"))
}
\arguments{
  \item{x}{

    Numeric vector of positive sample values. For the POT context, this
    should be the vector of excesses over the threshold.

  }
  \item{alternative}{

    Character string describing the alternative hypothesis

  }
}
\value{
  The LR statistic value.
}
\description{
  
  Likelihood Ratio statistic for the exponential distribution vs. GPD.

}
\details{
  
  The Likelihood-Ratio statistic is actually \eqn{W:=-2 \log
    \textrm{LR}}{-2 log LR} where LR is the ratio of the likelihoods
    \emph{exponential} to \emph{alternative distribution}.

}
\note{
  
  When the alternative is \code{"lomax"} or \code{"maxlo"}, the
  statistic has a distribution of \emph{mixed type} under the null
  hypothesis of exponentiality. This is a mixture of a distribution of
  continuous type (with positive values) and of a Dirac mass at LR =
  0. The probability mass \eqn{\textrm{Pr}\{\textrm{LR} = 0\}}{Pr[LR =
  0]} can be computed using the \code{\link{pGreenwood1}} function. More
  precisely, the probability mass is \code{pGreenwood1(n)} for the Lomax
  alternative and \code{1.0 - pGreenwood1(n)} for the maxlo alternative,
  where \code{n} is the sample size \code{length(x)}.


}
\seealso{

  \code{\link{LRExp.test}} for the related LR test of exponentiality.
  
}
