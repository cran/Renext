\name{RLplot}
\alias{RLplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Return level plot}
\description{
  Return level plot for "Renouvellement" data.
}
\usage{
RLplot(data,
       x = NULL,
       duration = 1,
       lambda,
       conf.pct = 95,
       mono = TRUE,
       mark.rl = 100,
       mark.labels = mark.rl,
       mark.col = NULL,
       main = NULL,
       ylim = NULL,
          ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    
    A data.frame object with a column named \code{quant}.

  }
  \item{x}{

    Optional vector of observed levels.

  }
  \item{duration}{

    The (effective) duration corresponding to \code{x} if
    this argument is used.

  }
  \item{lambda}{

    Rate, with unit inverse of that used for \code{duration}, e.g. in
    inverse years when \code{duration} is in years.

  }
  \item{conf.pct}{

    Vector (character or integer) giving confidence levels. See
    \bold{Details} below.

  }
  \item{mono}{

    If \code{TRUE} colours are replaced by black.

  }
  \item{mark.rl}{

    Return levels to be marked on the plot.

  }
  \item{mark.labels}{

    Labels shown at positions in \code{mark.rl}.

  }
  \item{mark.col}{

    Colours for marked levels.

  }
  \item{main}{

    Main title for the return level plot (defaults to empty title).

  }
  \item{ylim}{

    Limits for the y axis (defaults to values computed from the data).

  }
  \item{\dots}{

    Further args to be passed to \code{plot}. Should be removed in
    future versions.

  }
}
\details{

  Percents should match column names in the data.frame as follows. The
  upper and lower limits are expected to be \code{U.95} and \code{L.95}
  respectively. For a \code{70\%} confidence percentage, columns should
  have names \code{"U.70"} and \code{"L.70"}.

  The plot is comparable to the return level described in Coles'book and
  related packages, but the return level is here in log-scale while
  Coles uses a loglog-scale. A line corresponds here to a one parameter
  exponential distribution, while Coles'plot corresponds to Gumbel,
  however the two plots differ only for small return periods.  This plot
  is identical to an \code{expplot} but with x and y scales changed:
  only axis tick-marks differ. The convexity of the scatter plot is
  therefore opposed in the two plots.
  
}
%\value{ }
\references{
  
 Coles S. (2001)  \emph{Introduction to Statistical Modelling
   of Extremes Values}, Springer.

}
\author{Yves Deville}
\note{

  Confidence limits correspond to \emph{two-sided symmetrical
  intervals}. This means that the (random) confidence interval may be
  under or above the true unknown value with the same
  probabilities. E.g. the probability that the unknown quantile falls
  above \code{U.95} is \code{2.5\%}. The two bounds are yet generally
  not symmetrical with respect to \code{quant}; such a behaviour follows
  from the use of "delta" method for approximate intervals.

  It is possible to add graphical material (points, lines) to this plot
  using \code{log(returnlev)} and \code{quantile} coordinates. See
  \bold{Examples} section.

}
\seealso{
  
  See \code{\link{expplot}} for a classical exponential plot.  See Also
  as \code{\link{Renouv}} to fit "Renouvellement" models. The
  \code{return.level} function in the \code{extRemes} package.

}
\examples{
## Typical probability vector
prob <- c(0.0001,
  seq(from = 0.01, to = 0.09, by = 0.01),
  seq(from = 0.10, to = 0.80, by = 0.10),
  seq(from = 0.85, to = 0.99, by = 0.01),
  0.995, 0.996, 0.997, 0.998, 0.999, 0.9995)

## Model parameters rate = #evts by year, over nyear
lambda <- 4
nyear <- 30
theta.x <- 4

## draw points
n.x <- rpois(1, lambda = lambda*nyear)
x <- rexp(n.x, rate = 1/theta.x)

## ML estimation (exponential)
lambda.hat <- n.x / nyear
theta.x.hat <- mean(x)
  
## Compute bounds (here exact)
alpha <- 0.05

quant <- qexp(p = prob, rate = 1/theta.x.hat) 

theta.L <- 2*n.x*theta.x.hat / qchisq(1 - alpha/2, df = 2*n.x)
theta.U <- 2*n.x*theta.x.hat / qchisq(alpha/2, df = 2*n.x)

L.95 <- qexp(p = prob, rate = 1/theta.L) 
U.95 <- qexp(p = prob, rate = 1/theta.U) 

## store in data.frame object
data <- data.frame(prob = prob, quant = quant, L.95 = L.95, U.95 = U.95)

RLplot(data = data, x = x, lambda = lambda.hat,
       duration = nyear,
       main = "Poisson-exponential return levels")

RLplot(data = data, x = x, lambda = lambda.hat, duration = nyear,
       mark.rl = 10, mark.labels = "10 ans", mono = FALSE, mark.col = "SeaGreen",
       main = "Poisson-exponential return levels")

points(x = log(50), y = 25, pch = 18, cex = 1.4, col = "purple") 
text(x = log(50), y = 25, col ="purple", pos = 4, labels = "special event") 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%eyword{ ~kwd1 }
%%eyword{ ~kwd2 }% __ONLY ONE__ keyword per line
