\name{skip2noskip}
\alias{skip2noskip}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

  Fix non-skipped periods from skipped ones

}
\description{

  Compute non-skipped periods form start and end of skipped periods.  

}
\usage{
   skip2noskip(skip = NULL,
               start = NULL,
               end = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{skip}{

    A data.frame object with \code{start} and \code{end} columns that
    can be coerced to \code{POSIXct}. Other columns can be present (and will be
    ignored). Each row describes a missing period. Rows must be sorted
    in chronological order and periods should not overlap. Validity
    checks are at the time very limited.

  }
  \item{start}{

    Beginning of the whole period, to be used in \code{as.POSIXct}.

  }
  \item{end}{

    End of the whole period to be used in \code{as.POSIXct}.

  }
}
\details{
  
  In a 'normal' use of this function \code{start} and \code{end} are
  given, and are respectively \emph{before the beginning} of the first
  \code{skip} period and \emph{after the end} of the last \code{skip}
  period.  Thus the returned dataframe will have \code{nrow(skip)+1}
  rows. However, \code{start} and \code{end} can be \code{NULL} in which
  case only the \code{nrows(skip)-1} "inner" non-skipped periods will be
  returned. If \code{start} and \code{end} are \code{NULL} and
  \code{skip} has only one row, the returned result is \code{NULL}.

}
\value{

  A data.frame object with two \code{POSIXct} columns named \code{start} and
  \code{end}. Each row corresponds to a non-skipped period

}
\author{

  Yves Deville

}
%\note{ }
\seealso{

  \code{\link{readXML}} for reading data from XML and csv files.

}
\examples{
## Brest data embeds a description of the gaps

ns <- skip2noskip(skip = Brest$OTmissing)

ns2 <- skip2noskip(skip = Brest$OTmissing,
                   start = Brest$OTinfo$start,
                   end = Brest$OTinfo$end)

## check durations. dur2 should be equal to the effective
## duration (with an error of a fraction of day)
dur <- as.numeric(sum(ns$end-ns$start))/365.25
dur2 <- as.numeric(sum(ns2$end-ns2$start))/365.25


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
