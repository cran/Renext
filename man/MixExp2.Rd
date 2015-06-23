\name{MixExp2}
\alias{MixExp2}
\alias{dmixexp2}
\alias{rmixexp2}
\alias{pmixexp2}
\alias{qmixexp2}
\alias{hmixexp2}
\alias{Hmixexp2}
\title{Mixture of two exponential distributions}
\description{
  
  Probability functions associated to the mixture of
  two exponential distributions.

}  
\usage{
   dmixexp2(x, prob1,
            rate1 = 1.0, rate2 = rate1 + delta, delta,
            log = FALSE)
   pmixexp2(q, prob1,
            rate1 = 1.0, rate2 = rate1 + delta, delta,
            log = FALSE)
   qmixexp2(p, prob1,
            rate1 = 1.0, rate2 = rate1 + delta, delta)
   rmixexp2(n, prob1,
            rate1 = 1.0, rate2 = rate1 + delta, delta)
   hmixexp2(x, prob1,
            rate1 = 1.0, rate2 = rate1 + delta, delta)
   Hmixexp2(x, prob1,
            rate1 = 1.0, rate2 = rate1 + delta, delta)
}

\arguments{
  \item{x, q}{

    Vector of quantiles.

  }
  \item{p}{

    Vector of probabilities.

  }
  \item{n}{

    Number of observations.

  }
  \item{log}{

    Logical; if \code{TRUE}, the log density is returned.

  }
  \item{prob1}{

    Probability weight for the "number 1" exponential density.

  }
  \item{rate1}{

    Rate (inverse expectation) for the "number 1" exponential density.

  }
  \item{rate2}{

    Rate (inverse expectation) for the "number 2" exponential
    density. Should in most cases be \code{> rate1}. See \emph{Details}.

  }
  \item{delta}{

    Alternative parameterisation \code{delta = rate2 - rate1}.

  }
}

\value{
  
  \code{dmiwexp2}, \code{pmiwexp2}, \code{qmiwexp2}, evaluates the
  density, the distribution and the quantile functions.  \code{dmixexp2}
  generates a vector of \code{n} random draws from the distribution.
  \code{hmixep2} gives hazard rate and \code{Hmixexp2} gives cumulative
  hazard.

}

\details{
  
  The density function is the mixture of two exponential densities
  \deqn{
    f(x) = \alpha_1 \lambda_1 \, e^{-\lambda_1 x} + (1-\alpha_1)
    \lambda_2 \, e^{-\lambda_2x} \qquad x > 0
  }{
    f(x) = \alpha_1 \lambda_1 exp(- \lambda_1 x) +
          (1 - \alpha_1) \lambda_2 exp(- \lambda_2 x)     (x > 0)
  }
  where \eqn{\alpha_1} is the probability given in
  \code{prob1}
  while \eqn{\lambda_1} and\eqn{\lambda_2} are
  the two rates given in \code{rate1} and \code{rate2}.
  %% \code{prob1*depx(x, rate = rate1) + (1-prob1)*dexp(x, rate = rate2)}
  
  A 'naive' identifiability constraint is \deqn{\lambda_1 < \lambda_2}
  i.e.  \code{rate1 < rate2}, corresponding to the simple constraint
  \code{delta > 0}.  The parameter \code{delta} can be given instead of
  \code{rate2}.
  
  The mixture distribution has a decreasing hazard, increasing Mean
  Residual Life (MRL) and has a thicker tail than the usual
  exponential. However the hazard, MRL have a finite non zero limit and
  the distribution behaves as an exponential for large return
  levels/periods.
  
  The quantile function is not available in closed form and is computed
  using a dedicated numerical method.

}

\examples{
rate1 <- 1.0
rate2 <- 4.0
prob1 <- 0.8
qs <- qmixexp2(p = c(0.99, 0.999), prob1 = prob1,
               rate1 = rate1, rate2 = rate2) 
x <- seq(from = 0, to = qs[2], length.out = 200)
F <- pmixexp2(x, prob1 = prob1, rate1 = rate1, rate2 = rate2)
plot(x, F, type = "l", col = "orangered", lwd = 2,
     main = "Mixexp2 distribution and quantile for p = 0.99")
abline(v = qs[1])
abline(h = 0.99)

}

\keyword{distribution}





