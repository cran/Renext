\name{OT2MAX}
\alias{OT2MAX}
\title{Temporal aggregation of a Marked Process}
\usage{
  OT2MAX(OTdata,
         OTmissing = NULL,
         start = NULL,
         end = NULL,
         MAX.r = 1L,
         blockDuration = "year",
         monthGapStat = TRUE,
         maxMissingFrac = 0.05,
         dataFrames = FALSE,
         infMAX = FALSE,
         plot = TRUE,
         plotType = c("max", "gaps"),
         jitterSeed = 123,
         trace = 0L,
         ...)
}
\arguments{
  \item{OTdata}{

    Data frame containing a \code{POSIXct}
    column \code{date} and the marks variable.
    
  }
  \item{OTmissing}{

    Optional data frame with columns \code{start} and \code{end}
    (coerced to \code{POSIXct}) giving the beginning and the end of
    gaps.

  }
  \item{start}{
    
    An object coerced to \code{POSIXct} indicating the beginning of
    reliable/usable information. Unless this is a beginning of block
    (1-st of January for years), the 1-st block will begin \emph{after}
    \code{start} in order to use only qualified information.

  }
  \item{end}{

    An object indicating the end of the reliable/usable
    information. Unless this is a end of block (1-st of January for
    years), the last block will end \emph{before} \code{end} in order to
    use only qualified information.

  }
  \item{MAX.r}{

    \emph{Target} number of observations in the blocks. Can be of length
    one (same number of observations for all blocks) or of length equal
    to the number of blocks, the values being then for the blocks in the
    same order. In both cases, the target number may be impossible to
    reach because of a smaller number of events in the block. If
    \code{infMAX} is \code{TRUE}, the target number of observations will
    be reached by filling if needed with \code{-Inf} values. The
    rationale for this is that a non-existing event is assumed to have
    an arbitrarily small mark.
    
  }
  \item{blockDuration}{

    Duration of the blocks. Can only be \code{"year"} for now.
    
  }
  \item{monthGapStat}{

    Logical. Setting it to \code{TRUE} will compute statistics
    concerning the gaps and return them or show them on a plot.

  }
  \item{maxMissingFrac}{
    
    Maximal fraction of a block duration (between 0 and 1) that can be
    missing without leading to a \code{NA} aggregated value.
    
  }

  \item{dataFrames}{
    
    If \code{TRUE}, the result will contain data frames similar
    to those found in an object with class \code{"Rendata"}. If
    \code{FALSE} the result will contain \emph{list} and \emph{vector}
    objects, similar to those used as inputs in the \code{\link{Renouv}}
    function under the names \code{MAX.data} and
    \code{MAX.effDuration}. Note however, that \code{-Inf} values can be
    found in these objects when \code{infMAX} is \code{TRUE}.
  }
  \item{infMAX}{
    
    If \code{FALSE}, the target number of values the blocks will generally
    not be reached, because the total number of events in a block can be
    lower than the target number. Then, the target number value is revised
    to the number of found values in each block. If \code{TRUE}, the
    target number of values is reached by filling the values with
    \code{-Inf} and the datetimes with (\code{POSIXct}) \code{NA}s.
  }
  \item{plot}{
    
    If \code{TRUE} a simple plot is shown.
    
  }
  \item{plotType}{
    
    Character controlling the plot. With \code{"max"}, the block maxima
    are shown. With \code{plotType = "gap"}, the daily and monthly gap
    rates are shown. This is possible when suitable information
    concerning gaps is provided in \code{OTmissing}. The plot then shows
    the probability that a given day of the year falls in a gap, as well
    as monthly gap rates. Most often one wants that the gap rate does
    not show a seasonal behaviour. Note that gap rates for month-year
    combinations are shown as grey segments after jitterizing them since
    the values \code{0} and \code{1} may be observed for several
    years. An alternative way to is using the \code{monthGapTS}
    multivariate time series returned by the function, see
    \bold{Examples}.

  } \item{jitterSeed}{
    
    Random seed for jittering. Used only when \code{plot} is
    \code{TRUE}, \code{plotType} is \code{"gap"} and when suitable
    information is provided in \code{OTmissing}.
    
  }
  \item{trace}{

    Integer level of verbosity.

  }
  \item{...}{
    
    Other arguments to be passed to \code{plot}.

  }
  
}
\value{
  
  A list, the content of which depends on the value of
  \code{dataFrames}. If this value is \code{TRUE}, the
  following elements are returned.
  
  \item{MAXdata}{
    
    A data frame of largest values by block with one row for
    each observation. The largest values are given as columns
    with names equal to those in the \code{OTdata} data
    frame.
    
  }
  \item{MAXinfo}{
    
    A data frame describing the blocks, with one row by
    block. The two (\code{POSIXct}) columns \code{"start"}
    and \code{"end"} provide the beginning and the end of the
    block.  The numeric column \code{duration} gives the
    \emph{effective duration} (in year) within block.
    
  }
  \item{probMissing}{
    
    A vector with values corresponding to the days in a block
    (year). Each value is a estimation of the probability
    that the day falls in a gap.
   
  
    If \code{dataFrames} is \code{FALSE}, the list still
    contains \code{probMissing} as before, as well as other
    lists as used in \code{\link{Renouv}}.

  }
  \item{effDuration, r}{
    
    Vectors containing the effective duration (\emph{in
      years}) and number of value for the blocks.
    
  }
  \item{data}{
    
    List of maxima or \eqn{r}-largest values for the blocks.
    
  }
  \item{monthGapStat, monthGapTS}{

    Summary information concerning gaps, if \code{monthGapStat} is
    \code{TRUE} and if relevant information is provide via the the
    \code{OTmissing} formal. The element \code{monthGapTS} is a
    multivariate time series with yearly observations and one series
    (column) for each of the 12 months. Each series contains the missing
    fraction of the month for the considered year, ranging from
    \code{0.0} (no gap) to \code{1.0} (full gap). This object can be
    dealt with standard methods for time-series, but the \code{plot}
    method will require to select a reduced number of columns first,
    see \bold{Examples}.
    
  }

}
\description{

  Temporal aggregation of a Marked Process, leading to block maxima or
  \eqn{r}-largest observations.

}
\details{

  The data frame given in \code{OTdata} contains the \emph{events} (or
  \emph{arrivals}) given by the \code{date} column, as well as one mark
  column. Depending on the argument \code{MAX.r}, the maxima or the
  \eqn{r}-largest observations of the marks is computed for each time
  block. When known gaps exist in the data and when they are given in
  \code{OTmissing}, a block for which the total duration of gaps is too
  large will be omitted.
  
}

\note{

  Remind that even when \code{maxMissingFrac} is set to its maximum
  value 1.0, there can still be blocks with no data. When the result is
  intended to be used in the \code{\link{Renouv}} function, the formal
  \code{dataFrames} should be \code{FALSE}; the elements \code{data} and
  \code{effDuration} can then be passed as \code{MAX.data} and
  \code{MAX.effDuration}. At the time \code{infMAX} should also then be
  set to \code{FALSE} since \code{-Inf} values are not yet allowed in
  the \eqn{r}-largest values.
  
}
\examples{
## use Dunkerque data
OTdata <- Dunkerque$OTdata; OTmissing <- Dunkerque$OTmissing
## allow up to 50\\\% gap in a block, or only 5\\\%
MAX1 <- OT2MAX(OTdata = OTdata, OTmissing = OTmissing,
               maxMissingFrac = 0.5,
               main = "impact of the 'maxMissingFrac' formal")
MAX2 <- OT2MAX(OTdata = OTdata, OTmissing = OTmissing, dataFrames = TRUE,
               prefix = "Max", maxMissingFrac = 0.05, plot = FALSE)
lines(MAX2$MAXdata$date, MAX2$MAXdata$Surge, type = "h", col = "red", lwd = 3)
legend("topleft", lw = c(1, 3), col = c("black", "orangered"),
       legend = c("50\\\% max", " 5\\\% max"))

## r-largest obs for r = 4
MAX3 <- OT2MAX(OTdata, OTmissing = OTmissing, MAX.r = 4,
               maxMissingFrac = 0.9, 
               dataFrames = FALSE, trace = TRUE,
               main = "r-largest with r = 4")

## restrict the period
MAX4 <- OT2MAX(OTdata, OTmissing = OTmissing, MAX.r = 4,
               start = "1962-01-01",
               end = "1990-01-01",
               maxMissingFrac = 0.9, 
               dataFrames = FALSE, trace = TRUE,
               main = "r-largest with r = 4 with given 'start' and 'end'")
\dontrun{
  ## use in a block maxima analysis, as if there were no gaps.
  fit <- fGEV.MAX(MAX.data = MAX3$data,
                  MAX.effDuration = rep(1, length(MAX3$effDuration)))     
}
## plot the gap rate
MAX5 <- OT2MAX(OTdata = OTdata, OTmissing = OTmissing,
               maxMissingFrac = 0.5,
               main = "probability of being in a  gap",
               plotType = "gap")

## time series plot (only <= 10 months)
plot(MAX5$monthGapTS[ , c(1:4)], main = "gap rate by month")

## much better with lattice.
\dontrun{
    require(lattice)
    xyplot(MAX5$monthGapTS)
}
}