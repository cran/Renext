\name{plot.Renouv}
\alias{plot.Renouv}
\alias{lines.Renouv}
\title{
  Plot an object of class "Renouv"
}
\description{
 Plot an object of class "Renouv". The plot is a return level plot with
  some supplementary elements to display historical data.
}
\usage{
\S3method{plot}{Renouv}(x,
    pct.conf = x$pct.conf,
    show = list(OT = TRUE, quant = TRUE, conf = TRUE,
                MAX = TRUE, OTS = TRUE),
    mono = TRUE,
    predict = FALSE,
    par = NULL,
    legend = TRUE,
    legendEnvir = NULL,
    label = NULL,
    problim = NULL,
    Tlim = NULL,
    main = NULL, xlab = "periods", ylab = "level",
    ...)
\S3method{lines}{Renouv}(x,
    pct.conf = x$pct.conf,
    show = NULL,
    mono = TRUE,
    predict = FALSE,
    par = NULL,
    legend = FALSE,
    legendEnvir = NULL,
    label = NULL,
    ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Object of class \code{"Renouv"}.
  }
  \item{pct.conf}{
   Percents for confidence limits (lower and
   upper). These levels should be found within those computed in the
   object \code{x}. By default, all computed levels will be used.
  }
  \item{show}{
    A list with named elements specifying which parts of the return
    level plot must be drawn. Element \code{OT} is for the
    the sample points (Over the Threshold data), \code{quant}
    is for the quantile curve (or Return Level curve), \code{conf}
    is for the confidence limits. Finally \code{MAX} and \code{OTS}
    are for the two possible types of historical data. When
    the element \code{conf} is \code{TRUE}, only the percent
    levels given in \code{pct.conf} are drawn. Moreover, the levels
    not already computed in the object given in \code{x} will be drawn
    only if the predictions are recomputed with \code{predict} set
    to \code{TRUE}.
  }
  \item{mono}{
    Logical, \code{TRUE} for a monochrome plot.
  }
  \item{predict}{Logical. When \code{TRUE}, predictions are re-computed
    from the model before plotting. One effect is that the points used
    to draw the curves are designed to cover the whole range (if
    specified by the user). One other effect is that the confidence
    limits are recomputed in order to include the percent levels
    given on entry.
  }
  \item{par}{
    A list such as returned by the \code{\link{RLpar}} function.
  }
  \item{legend}{
    Logical. If \code{TRUE}, a legend is built and drawn on the graph. 
  }
  \item{legendEnvir}{
    An optional environment in which legend information will be stored
    for later use. This information will be shown by using
    \code{\link{RLlegend.show}}. This environment is used only when
    \code{legend} is \code{FALSE}, and the legend is drawn directly otherwise.
  }
  \item{label}{
    A character label used to build the labels used in the legend. 
  }
  \item{problim}{
    Limits for the x-axis in probability scale. Can be used as an
    alternative to \code{Tlim}.
  }
  \item{Tlim}{
    Limits for the x-axis in return period scale. The values are given
    as a numeric vector of length 2, containing values
    \eqn{\ge 1}{>=1}. 
    The first element (minimal return period can be 0 in which case
    it will be replaced by a very small positive value.
  }
  \item{xlab}{
    Label of the x-axis (time periods, with log scale).
  }
  \item{main}{
    Main title (character).
  }
  \item{ylab}{
    Label of the y-axis (labels).
  }
  \item{\dots}{
    Other arguments passed to the default \code{\link{plot}} function
    e.g., \code{ylim} to adjust the y-axis.
  }
}
\details{
 
  Historical data embeded in the \code{x} object are \code{list(MAX =
    TRUE)}
  
  \itemize{
    \item{If the \code{MAX} element is \code{TRUE} and if \code{x} embeds
      historical data of type \code{MAX}, then these will be shown
      with a symbol differing from the one for ordinary points.}
    
    \item{If \code{OTS} element is \code{TRUE} and is \code{x} embeds
      historical data of type \code{OTS}, then these will be shown
      with a symbol differing from the one for ordinary points. An
      exception is when one or several OTS block have no data. Then each
      such block is shown as an horizontal segment; its right end-point
      shows the effective duration of the block and the ordinate shows
      the OTS threshold for this block. No data exceeded the threshold
      within the block.}
  }

  This function acts on a list variable named \code{.RLlegend} in the
  \code{legendEnvir} environment. This variable is used
  to build legends for plots produced with multiple commands. See
  the \code{\link{RLlegend}} help page. 
  
}
\value{
  No value returned.
}
%%\references{ }
\author{
  Yves Deville
}
\note{
  The return level plot is of exponential type i.e. uses a log-scale for
  return periods. This contrasts with the Gumbel plot which is also used
  in similar contexts.
}

\section{Caution}{
  Remaind that the methods \code{plot} and \code{lines} may change the
  value of the variable \code{.RLlegend} in the environment
  \code{legendEnvir}. This variable describes the material to be used in
  the legend at the next call of \code{RLlegend.show}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  The
  \code{\link{RLlegend}} page for the legend construction and
  \code{\link{RLpar}} to specify the grahical parameters (colors, line
  types, ...) of the elements.
}
\examples{
## two fits for the Garonne data
fit.exp <- Renouv(x = Garonne, plot = FALSE)
fit.gpd <- Renouv(x = Garonne, distname.y = "gpd", plot = FALSE)

## simple plot (legend is TRUE here)
plot(fit.exp,
     main = "Two fits overlapped",
     label = "",
     ## Tlim = c(1, 5000),
     predict = TRUE)

## Now try 'lines' and RLlegend.xxx functions
plot(fit.exp,
     main = "Fancy legend",
     show = list(OT = FALSE, quant = FALSE, conf = FALSE,
                 OTS = FALSE, MAX = FALSE),
     legend = FALSE,
     Tlim = c(1, 5000))

legEnv <- new.env()
RLlegend.ini(x = "bottomright", bg = "lightyellow", envir = legEnv) ## initialise legend

lines(fit.exp,
      show = list(quant = FALSE, conf = FALSE, OT = TRUE, MAX = TRUE),
      label = "expon",
      par = RLpar(quant.col = "orange", 
        OT.pch = 21, OT.cex = 1.2, OT.col = "SeaGreen", OT.bg = "yellow",
        MAX.block1.col = "purple", MAX.block1.bg = "mistyrose",
        MAX.block1.lwd = 1.4),
      legendEnvir = legEnv
      )

lines(fit.gpd,
      pct.conf = c(95, 70),
      show = list(quant = TRUE, conf = TRUE),
      label = "GPD",
      par = RLpar(quant.col = "darkcyan", conf.conf1.col = "red"),
      legendEnvir = legEnv
      )

RLlegend.show(envir = legEnv) ## now draw legend
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
