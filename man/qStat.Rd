\name{qStat}
\alias{qStat}
\encoding{UTF8}
\title{Quantiles of a test statistic}
\usage{
   qStat(p, n,
         type = c("Greenwood", "Jackson", "logLRGPD", "logLRLomax",
                  "logLRGEV", "logLRFrechet"),
          outNorm = FALSE)
}
\arguments{

  \item{p}{

    Numeric vector of probabilities. Very small values (\code{p < 0.01})
    or very large ones (\code{p > 0.99}) will be truncated as
    \code{0.00} or \code{1.00} to maintain a realistic level of
    precision.

  }
  \item{n}{

    Sample size.

  }
  \item{type}{

    The type of statistic, see \bold{Details}.

  }
  \item{outNorm}{

    Logical. If \code{TRUE} the output is normalized in a such fashion
    that its distribution is the asymptotic one (i.e. standard normal in
    practice). When \code{FALSE}, the quantiles are given in the true
    scale of the statistic: \eqn{\textrm{CV}^2}{CV^2}, Jackson. For LR
    statistics this argument has no impact.

  }
}
\value{
  A vector of quantiles.
}
\description{
  Quantile of a test statistic.
}
\details{
  
  The function provides an approximation of the distribution for several
  statistics.
  
  \itemize{
    
    \item For \code{"Greenwood"}, the statistic is \emph{Greenwood's
    statistic}. The distribution is that of the squared coefficient of
    variation \eqn{\textrm{CV}^2}{CV^2} of a sample of size \code{n}
    from the exponential distribution as computed by \code{\link{CV2}}.
    
    \item For \code{"Jackson"}, the statistic is Jackson's
    statistic, see \code{\link{Jackson}}.
    
    \item For \code{"logLRGPD"} and \code{"logLRLomax"}, the
    statistic is the log of the likelihood ratio of a sample
    from the exponential distribution. The log-likelihoods are
    for an exponential distribution compared to a GPD with
    non-zero shape, or to a GPD with \emph{positive shape}
    (equivalently, a Lomax distribution).
    
    \item For \code{"logLRGEV"} and \code{"logLRFrechet"}, the
    statistic is the log of the likelihood ratio of a sample
    from the Gumbel distribution. The log-likelihoods are for a
    Gumbel distribution compared to a GEV with non-zero shape,
    or to a GEV with \emph{positive shape} (equivalently, a
    Fréchet distribution).
    
  }
  
  The log of Likelihood Ratios are multiplied by \code{2}, so that they
  compare to a chi-square statistic with one degree of freedom.
  
}
\note{
  
  The precision of the result given is limited, and is about
  two-digits. This function is not intended to be used as such and is
  only provided for information.

}
\examples{
res <- qStat(n = 40, type = "Greenwood")
plot(res$q, res$p, type = "o")
}
\author{
  Yves Deville
}
