\name{interevt}
\alias{interevt}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Interevents (or interarrivals) from events dates}

\description{
  
  Compute intervent durations from events dates

}
\usage{
   interevt(date,
            skip = NULL, noskip = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{date}{

    A \code{POSIXct} vector containing the date(time) of the events.

  }
  \item{skip}{

    A data.frame containing two \code{POSIXct} columns \code{start} and
    \code{end} describing the periods to skip over.

  }
  \item{noskip}{

    A data.frame like \code{skip} but where the periods define the NON
    skipped part of the events.

  }
}
\details{
  
   Interevents are the time differences between successive dates. When
   the \code{date} argument contains occurrence times \eqn{T_i}{T[i]}
   for successive events of an homogeneous Poisson process, interevents
   \eqn{T_i -T_{i-1}}{T[i] - T[i-1]} are mutually independent with the
   same exponential distribution.

   When some time intervals are skipped independently from the event
   point process, we may consider the interevents
   \eqn{T_i-T_{i-1}}{T[i] - T[i-1]} between two non-skipped events such
   that the time interval \eqn{(T_{i-1},\,T_i)}{(T[i-1], T[i])} does not
   contains any skipped interval.  These interevents still are mutually
   independent with the same exponential distribution.  When \code{skip}
   or \code{noskip} is not \code{NULL} the computation therefore only
   retains couples of two successive datetimes "falling" in the same
   non-skipped period, which number can therefore be associated with the
   interevent.

 }
 \value{
   
  A list mainly containing a \code{interevt} data.frame.

  \item{interevt}{

    Data.frame. Each row describes a retained interevent through a
    \code{period} integer giving the "noskip" period, a \code{start} and
    \code{end} \code{POSIXct} and a \code{duration} in \bold{days}.

  }
  \item{noskip}{

    Only when \code{skip} or \code{noskip} args have been given. A
    data.frame containing broadly the same information as the
    \code{noskip} arg is it was given or the information deduced from
    the \code{skip} arg if given.

  }
  \item{axis}{

    When needed, a list with some material to build an axis with uneven
    ticks as in the \code{gof.date} with \code{skip.action = "omit"}.

  } 
}
%%\references{ }
\author{

  Yves Deville

}
\note{

  Only one of the two arguments \code{skip} and \code{noskip} should be
  given in the call. In each case, the rows of the returned data.frame
  objects describe periods in chronological order. That is: \code{start}
  at row \code{2} must be after the \code{end} value of row \code{1} and
  so on.

  Note that there are usually less interevents than dates since two
  successive dates will be retained for an interevent only when they are
  not separated by missing period. As a limit case, there can be no
  interevents if the \code{noskip} periods contain only one date from
  the \code{date} vector.

}
\seealso{

  \code{\link{gof.date}} for goodness-of-fit diagnostics for dates of
  events \code{\link{expplot}} for diagnostics concerning the
  exponential distribution.

}
\examples{
## Use Brest data
ie <- interevt(date = Brest$OTdata$date, skip = Brest$OTmissing)

expplot(ie$interevt$duration, rate = 1 / mean(ie$interevt$duration),
  main = "No threshold")

## keep only data over a threshold
ind1 <- Brest$OTdata$Surge >= 35
ie1 <- interevt(Brest$OTdata$date[ind1], skip = Brest$OTmissing)
expplot(ie1$interevt$duration, main = "Threshold = 35")

## increase threshold
ind2 <- Brest$OTdata$Surge >= 55
ie2 <- interevt(date = Brest$OTdata$date[ind2], skip = Brest$OTmissing)
expplot(ie2$interevt$duration, main = "Threshold = 55 cm")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
