\name{expplot}
\alias{expplot}
\title{

  Classical "exponential distribution" plot

}
\description{
  
  Plot a vector using "exponential distribution" scales

}
\usage{
expplot(x,
        plot.pos = "exp",
        rate = NULL,
        labels = NULL,
        mono = TRUE,
        ...)
}
\arguments{
  \item{x}{

    The vector to be plotted.

  }
  \item{plot.pos}{

    Plotting position for points: either "exp" for \emph{expected} ranks
    or "med" for a \emph{median} rank approximation (see \bold{Details}
    below).

  }
  \item{rate}{

    Rate parameter for one or several "exponential distribution" lines
    to be plotted
    
  }
  \item{labels}{
    Text to display in legend when "exponential distribution" lines
    are specified
  }
  \item{mono}{
    
    Monochrome graph?
    
  }
  \item{\dots}{

    Arguments to be passed to \code{plot}.

  }
}
\details{
  
  This plot shows \eqn{-\log[1-F(x)]}{-log[1-F(x)]} against \eqn{x}{x}
  where \eqn{F(x)}{F(x)} at point \eqn{i}{i} is taken as
  \eqn{i/(n+1)}{i/(n+1)} if \code{plot.pos} is \code{"exp"}, or as the
  "median rank" approximation \eqn{(i-0.3)/(n+0.4)}{(i-0.3)/(n+0.4)} if
  \code{plot.pos} is \code{"med"}.

  If the data in \code{x} is a sample from an exponential distribution,
  the points should be roughly aligned. However the largest order
  statistics have high sampling dispersion.

}
\author{

  Yves Deville

}
\note{

  The log scale for y is emulated via the construction of suitable
  graduations. So be careful when adding graphical material (points,
  etc) to this graph with functions of the "add to plot" family
  (\code{points}, \code{lines}, ...).

  The ML estimate of the \code{rate} parameter is the inverse of the
  sample mean.

}
\seealso{
  
  The \code{\link{weibplot}} function for a classical "Weibull" plot.
  The \code{\link{interevt}} is useful to compute interevents (or
  "interarrivals") that should follow an exponential distribution in the
  homogeneous Poisson process context.

}
\examples{
 x <- rexp(200)
 expplot(x, rate = 1/mean(x), labels = "fitted")
}
%% \keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
