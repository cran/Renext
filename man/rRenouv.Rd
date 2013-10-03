\name{rRenouv}
\alias{rRenouv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of "Renouvellement" data}
\description{
  Simulation of "Renouvellement" data, i.e. observations Over
  a Threshold and their counts on  blocks (e.g. years).
}
\usage{
    rRenouv(densfun.y = "exponential",
            par.y = list(rate = 1),
            densfun.N = "poisson",
            par.N = list(lambda = 6),
            threshold = 0,
            aggreg = TRUE,
            nb = 50,
            labb = seq(to = 2009, by = 1, length = nb),
            w = rep(1, nb))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{densfun.y}{A character string specifying the distribution of the
    exceedances over the threshold. At the time only \code{"exponential"}, \code{"weibull"}
    and \code{"gpd"} are available.}
  \item{par.y}{Named list giving parameter values for the
    \code{"y"}-distribution. See examples below.}
  \item{densfun.N}{Character string specifying the distribution for the
    counts. At the time only "poisson" and "negative binomial" are allowed.}
  \item{par.N}{Named list giving parameter values for the
    \code{N}-distribution for counts. The names are \code{lambda} in the
    Poisson case (\code{"poisson"}), and \code{size} and \code{prob} in
    the negative binomial case (\code{"negative binomial"} case).}
  \item{threshold}{Threshold value.}
  \item{aggreg}{Only \code{TRUE} is possible at the time.}
  \item{nb}{Number of blocks (or time intervals).}
  \item{labb}{Numeric vector of length \code{nb} that will be used in
    replacement of the block numbers vector.  Typically it can contain
    year numbers. Use \code{NULL} to obtain blocks from \code{1} to \code{nb}.}
  \item{w}{Vector of blocks (time interval) length, i.e. duration.}
}
%\details{ }
\value{
  A list with the following objects
  \item{x}{Vector of \code{x} values i.e. \code{threshold} plus exceedances.}
  \item{N}{Vector of counts.}
  \item{block}{Vector of length \code{length(x)} giving the block number
    for the corresponding element in \code{x}. When coerced to a factor
    \code{block} has \code{length(N)} levels.}
}
%% \references{  }
\author{Yves Deville}
\note{
At the present time, the random drawings of the \code{gpd} distribution
are generated with the the \code{rgpd} of the \code{evd} package.

In future versions, this function might return an object with a special
class (with name such as "renouv"). Then classical \code{simulate}, \code{plot},
...  methods could be defined for that class.

}
\seealso{\code{\link{Renouv}} to fit
  Renouvellement models.}
\examples{
\dontrun{
test1 <- rRenouv(nb = 100,
                threshold = 40,
                par.N = list(lambda = 2),
                densfun.y = "gpd",
                par.y = mom2par("gpd", mean = 30, sd = 36))
}

test2 <- rRenouv(nb = 100,
                threshold = 40,
                par.N = list(lambda = 2),
                densfun.y = "weibull",
                par.y = mom2par("weibull", mean = 30, sd = 36))

test3 <- rRenouv(nb = 100,
                threshold = 40,
                densfun.N = "negative binomial",
                par.N = list(gamma = 10, prob = 0.7),
                densfun.y = "weibull",
                par.y = mom2par("weibull", mean = 30, sd = 36))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
