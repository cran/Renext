\name{weibplot}
\alias{weibplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Classical Weibull distribution plot}
\description{
  Plots a vector using Weibull distribution scales
}
\usage{
   weibplot(x,
            plot.pos = "exp",
            shape = NULL, scale = NULL,
            labels = NULL,
            mono = TRUE,
            ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The vector to be plotted.}
  \item{plot.pos}{plotting position for points: either "exp" for
    \emph{expected} ranks or "med" for a \emph{median} rank approximation (see \bold{Details} below).}
  \item{shape}{Shape parameter for one or several Weibull lines to be plotted.}
  \item{scale}{Scale parameter for one or several Weibull lines to be plotted.}
  \item{labels}{Text to display in legend when Weibull lines are
    specified.}
  \item{mono}{Monochrome graph.}
  \item{\dots}{Arguments to be passed to \code{plot}.}
}
\details{
  This plot shows \eqn{\log\{-\log[1-F(x)]\}}{log\{-log[1-F(x)]\}} against
  \eqn{\log(x)}{log(x)} where \eqn{F(x)}{F(x)} at point \eqn{i}{i}
  is taken as \eqn{i/(n+1)}{i/(n+1)} if \code{plot.pos} is \code{"exp"}, or as
  the "median rank" approximation \eqn{(i-0.3)/(n+0.4)}{(i-0.3)/(n+0.4)}
  if \code{plot.pos} is \code{"med"}.
}
%%\value{}
%% \references{}
\author{Yves Deville}
\note{
  The graph displayed uses a log scale for x. The log-log scale for y is
  emulated via the construction of suitable graduations. So be careful when
  adding graphical material (points, etc) to this graph with functions of
  the "add to plot" family (\code{points}, \code{lines}, ...).
}
\seealso{The \code{\link{expplot}} function for an "exponential
  distribution" plot (dedicated to the \code{shape = 1} case), and 
  the \code{\link{fweibull}} function for ML estimation of the
  parameters.
}
\examples{
x <- rweibull(200, shape = 1.2, scale = 1)
weibplot(x, main = "Classical Weibull plot")
## Weibull lines
weibplot(x, shape = c(0.9, 1.3), scale = 1)
weibplot(x, shape = c(0.9, 1.3), scale = 1,
         labels = c("before", "after"))
weibplot(x, shape = c(0.9, 1.3), scale = 1,
         labels = c("before", "after"),
         mono = TRUE)
}
%% \keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
