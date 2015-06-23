\name{NBlevy}
\alias{NBlevy}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Negative Binomial Levy process}
\description{
  
   Negative Binomial Lévy process estimation from partial observations (counts)

 }
\usage{
NBlevy(N,
       gamma = NA,
       prob = NA,
       w = rep(1, length(N)),
       sum.w = sum(w),
       interval = c(0.01, 1000),
       optim = TRUE,
       plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{

    Vector of counts, one count by time period.

  }
  \item{gamma}{

    The \code{gamma} parameter if known (NOT IMPLEMENTED YET).

  }
  \item{prob}{

    The \code{prob} parameter if known (NOT IMPLEMENTED YET).

  }
  \item{w}{

    Vector of time length (durations).

  }
  \item{sum.w}{

    NOT IMPLEMENTED YET. The effective duration. If \code{sum.w} is
    strictly inferior to \code{sum(w)}, it is to be understood that
    missing periods occur within the counts period. This can be taken
    into account with a suitable algorithm (Expectation Maximisation,
    etc.)

  }
  \item{interval}{

    Interval giving min and max values for \code{gamma}.

  }
  \item{optim}{

    If \code{TRUE} a one-dimensional optimisation is used. Else the zero
    of the derivative of the (concentrated) log-likelihood is searched
    for.

  }
  \item{plot}{

    Should a plot be drawn? \emph{May be removed in the future}.

  }
  \item{\dots}{Arguments passed to \code{plot}.}
}
\details{

  The vector \eqn{\mathbf{N}}{N} contains counts for events occurring on
  non-overlapping time periods with lengths given in
  \eqn{\mathbf{w}}{w}.  Under the NB \enc{Lévy}{Levy} process
  assumptions, the observed counts (i.e. elements of
  \eqn{\mathbf{N}}{N}) are independent random variables, each following
  a negative binomial distribution.  The size parameter \eqn{r_k}{r[k]}
  for \eqn{N_k}{N[k]} is \eqn{r_k = \gamma w_k}{r[k] = \gamma * w[k]} and
  the probability parameter \eqn{p}{p} is \code{prob}.  The vector
  \eqn{\boldsymbol{\mu}}{\mu} of the expected counts has elements
  \deqn{\mu_k=\mathrm{E}(N_k)=\frac{1-p}{p} \,\gamma \,w_k.}{\mu[k] =
  E(N[k]) = \gamma * w[k] * (1 - prob) / prob}
  
  The parameters \eqn{\gamma}{\gamma} and \eqn{p \:(\code{prob})}{prob}
  are estimated by Maximum Likelihood using the likelihood concentrated
  with respect to the \code{prob} parameter.

}
\value{
  A list with the results
  \item{estimate}{

    Parameter estimates.

  }
  \item{sd}{

    Standard deviation for the estimate.

  }
  \item{score}{

    Score vector at the estimated parameter vector.

  }
  \item{info}{

    Observed information matrix.

  }
  \item{cov}{

    Covariance matrix (approx.).

  }
}
\references{
  Kozubowski T.J. and \enc{Podgórsky}{Podgorsky} K. (2009)
  "Distributional properties of the negative binomial \enc{Lévy}{Levy} process".
  \emph{Probability and Mathematical Statistics} \bold{29}, pp. 43-71.
  Lund University Publications.
}
\author{Yves Deville}
\note{
  
  The Negative Binomial \enc{Lévy}{Levy} process is an alternative to
  the Homogeneous Poisson Process when counts are subject to
  overdispersion. In the NB process, all counts share the same index of
  dispersion (variance/expectation ratio), namely \code{1/prob}. When
  \code{prob} is close to 1, the counts are nearly Poisson-distributed.

}
\seealso{\code{\link{NegBinomial}} for the negative binomial distribution,
  \code{\link{glm.nb}} from the MASS package for fitting Generalised
  Linear Model of the negative binomial family.}
\examples{
## known parameters
nint <- 100
gam <- 6; prob <- 0.20

## draw random w, then the counts N
w <- rgamma(nint, shape = 3, scale = 1/5)
N <- rnbinom(nint, size = w * gam, prob = prob)
mu <- w * gam * (1 - prob) / prob
Res <- NBlevy(N = N, w = w)

## Use example data 'Brest'
## compute the number of event and duration of the non-skipped periods
gof1 <- gof.date(date = Brest$OTdata$date,
                 skip = Brest$OTmissing,
                 start = Brest$OTinfo$start,
                 end = Brest$OTinfo$end,
                 plot.type = "omit")
ns1 <- gof1$noskip
## fit the NBlevy
fit1 <- NBlevy(N = ns1$nevt, w = ns1$duration)

## use a higher threshold
OT2 <- subset(Brest$OTdata, Surge > 50)
gof2 <- gof.date(date = OT2$date,
                 skip = Brest$OTmissing,
                 start = Brest$OTinfo$start,
                 end = Brest$OTinfo$end,
                 plot.type = "omit")
ns2 <- gof2$noskip
## the NBlevy prob is now closer to 1
fit2 <- NBlevy(N = ns2$nevt, w = ns2$duration)

c(fit1$prob, fit2$prob)
}
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
