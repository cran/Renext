\name{plot.Renouv}
\alias{plot.Renouv}
\title{
  Plot an object of class "Renouv" 
}
\description{
  Plot an object of class "Renouv". The plot is a return level plot with
  some supplementary elements to display historical data.
}
\usage{
   \method{plot}{Renouv}(x,
     pct.conf = NULL,
     show.MAX = TRUE, show.OTS = TRUE,
     mono = TRUE,
     rl.mark = NULL, labels.mark = rl.mark, col.mark = NULL,
     main = NULL,
     problim = NULL,
     Tlim = NULL,
     xlab = "periods",
     ylab = "level",
     \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{"Renouv"}.}
  \item{pct.conf}{Percents for confidence limits (lower and
    upper). These levels should be found within those computed in the
    object \code{x}. By default, all computed levels will be used.}
  \item{show.MAX}{Logical. If \code{TRUE} and is \code{x} embeds
    historical data of type \code{MAX}, then these will be shown
    with a symbol differing from the one for ordinary points.
  }
  \item{show.OTS}{Logical. If \code{TRUE} and is \code{x} embeds
    historical data of type \code{OTS}, then these will be shown
    with a symbol differing from the one for ordinary points. An
    exception is when one or several OTS block have no data. Then each
    such block is shown as an horizontal segment; its right end-point
    shows the effective duration of the block and the ordinate shows
    the OTS threshold for this block. No data exceeded the threshold
    within the block.
  }
  \item{mono}{Logical. For a monochrome plot.}
  \item{rl.mark}{Each return level is materialised with a vertical line.}
  \item{labels.mark}{Label(s) for marked return levels.}
  \item{col.mark}{Color(s) for marked return levels.}
  \item{main}{Main title.}
  \item{problim}{Limits for the x-axis in probability scale.}
  \item{Tlim}{Limits for the x-axis in return period scale. The values
    are given as a numeric vector of length 2, containing values
    >=0. The first element (minimal return period can be 0 in which case
    it will be replaced by a very small positive value.}
  \item{xlab}{Label of the x-axis (log scale).}
  \item{ylab}{Label of the y-axis.}
  \item{\dots}{Further arguments to be passed to \code{plot}
    e.g. \code{ylim} to adjust the y-axis.}
}
\details{
  The return level plot is of exponential type i.e. uses a log-scale for
  return periods.
}
%%\value{}
%%\references{}
%%\author{}
%%\note{}
%%\seealso{}
\examples{
data(Brest)
fitB <- Renouv(Brest)
plot(fitB, main = "POT fit for Brest data", xlab = "years")
## Garonne example contains historical data.
data(Garonne)
fitG <- Renouv(Garonne)
plot(fitG, main = "POT fit for Garonne data", xlab = "years")
plot(fitG, Tlim = c(1, 2000))
}
% 
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
