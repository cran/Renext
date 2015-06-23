\name{RLlegend}
\alias{RLlegend}
\alias{RLlegend.ini}
\alias{RLlegend.show}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  
  Legend management for return level plots

}
\description{
  
  Legend management for return level plots produced with the \code{plot}
  and \code{lines} method of the \code{"Renouv"} class.

}
\usage{
    RLlegend.ini(x = "topleft", bty = "n", ...)
    RLlegend.show()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    
    A possible value for the \code{x} argument of \code{\link{legend}}.
    This will usually be a character giving the position e.g,
    \code{"topleft"} or \code{"bottomleft"}. See the
    \code{\link{legend}} function help.

  }
  \item{bty}{

    As in \code{\link{legend}}. The default value \code{"n"} differs
    from the default value of \code{legend}.

  }
%%  \item{envir}{
%%    An environment in which a variable \code{.RLlegend} will be
%%    stored. This argument is for compatibility reasons and is normally
%%    not used.
%%  }
  \item{\dots}{

    Other arguments to be kept in the list and passed later to
    \code{\link{legend}}. These arguments should be chosen among those
    of \code{legend} modifying the global legend appearance (e.g.,
    \code{bg}) but not among those modifying the legend content
    (e.g. \code{col} \code{pt.bg}, \code{legend}, ...) since the content
    is here built semi-automatically.

  }
}
\details{
  
  This function is to be used in conjunction with
  \code{\link{plot.Renouv}} and \code{\link{lines.Renouv}} methods. It
  allows the construction of a legend in a semi-automatic fashion, using
  the value of the \code{par} argument of the \code{plot} and
  \code{lines} methods to specify the legend construction.
  
  Each call to the \code{\link{plot.Renouv}} or
  \code{\link{lines.Renouv}} changes the content of a list variable
  named \code{.RLlegend} in a special environment bound to the
  package. This list is re-created when \code{RLlegend.ini} is called,
  and is used later to draw a legend on the active device when
  \code{RLlegend.show} is called.  Between these two calls, the
  \code{plot} and \code{lines} methods should be used with their arg
  \code{legend} set to \code{FALSE}.

  %% The list variable used to store the needed information can
  %% be assigned into the global environment at the user responsibility
  %% by a suitable choice of \code{envir}.

}
\value{
  
  \code{RLlegend.ini} returns a copy of the variable which is set.
  
  \code{RLlegend.show} returns nothing.

}
\author{
  
  Yves Deville

}
\note{
  The size of symbols (i.e, \emph{plotting characters}) can be set by
  using the \code{\link{RLpar}} function and the \code{\link{par}}
  argument of the methods \code{\link{plot.Renouv}} and
  \code{\link{lines.Renouv}}. However it can not be changed in the
  legend.
}
\seealso{
  \code{\link{plot.Renouv}} and \code{\link{lines.Renouv}} for
  and the \code{\link{RLpar}} function to change the graphical
  parameters of the plot and the legend by using the \code{par}
  argument.
}
\examples{
## use Garonne data
xG <- Garonne$OTdata$Flow
## use special "exponential" distribution 
fit1 <- Renouv(x = xG, threshold = 2500, distname.y = "exponential",
               effDuration = 65, plot = FALSE)

## use 'exp' in black box fashion, hence with delta method
fit2 <- Renouv(x = xG, , threshold = 2500, distname.y = "exp",
               effDuration = 65, start.par.y = c(rate = 1), plot = FALSE)
RLlegend.ini() ## initialise legend
## sample points only
plot(fit1, main = "Two types of confidence lims",
     show = list(OT = TRUE, quant = FALSE, conf = FALSE),
     label = "",
     legend = FALSE)
## quant and confidence lims
lines(fit1,
     show = list(OT = FALSE, quant = TRUE, conf = TRUE),
     label = "exact",
     legend = FALSE)
## quant (overplot) and confidence lims
lines(fit2,
      show = list(OT = FALSE, quant = TRUE, conf = TRUE),
      par = RLpar(quant.lty = 2, quant.col = "SpringGreen2",
        conf.conf1.col = "orangered", conf.conf1.lwd = 3,
        conf.conf2.col = "orangered", conf.conf2.lwd = 3),
      label = "delta",
      legend = FALSE)
RLlegend.show() ## now draw legend

}
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }
