\name{Brest}
\alias{Brest}
\encoding{UTF-8}
\docType{data}
\title{Surge heights at Brest}
\description{
  
  Surge heights near high tide at Brest tide gauge station (France),
  detailed version

}
\usage{Brest}
\format{
  The format is:
  List of 5
  \itemize{
    \item \code{$info} : List of 6
    \itemize{
      \item \code{$name}        : chr \code{"Brest"} 
      \item \code{$shortLab}    : chr \code{"Surge Heights at Brest (France)"} 
      \item \code{$longLab}     : chr \code{"Surge Heights near high tide,  Brest (France)"}  
      \item \code{$varName}     : chr \code{"Surge"}
      \item \code{$varShortLab} : chr \code{"Surge"}
      \item \code{$varUnit} : chr \code{"cm"}
    }
    \item \code{$describe} : chr
    "High tide sea surge over 30 cm at Brest (France)..."
    \item \code{$OTinfo} : List of 4
    \itemize{ 
      \item \code{$start}       : chr POSIXct[1:1], format: \code{"1846-01-01"}
      \item \code{$end}         : chr POSIXct[1:1], format: \code{"2009-01-01"}
      \item \code{$effDuration} : num \code{148}
      \item \code{$threshold}   : num \code{30}
    }
    \item \code{$OTdata} : 'data.frame':	1289 obs. of  2 variables:
    \itemize{
      \item \code{$date} : POSIXct[1:1289], format: \code{"1846-01-14" "1846-01-21"} ...
      \item \code{$Surge} : num [1:1289] \code{36 60 46 40 33} ...
    }
    \item \code{$OTmissing} : 'data.frame':	43 obs. of  3 variables:
    \itemize{
      \item \code{$start}   : POSIXct[1:43], format: \code{"1846-01-01" "1847-01-01"} ...
      \item \code{$end}     : POSIXct[1:43], format: \code{"1846-01-04" "1847-01-21"} ...
      \item \code{$comment} : chr [1:43] \code{"" "" "" ""} ...
    }
    - \code{attr(*, "class")= chr "Rendata"}
  }
}
\details{
  
  Data are provided as a list.
  
  \itemize{
    
    \item \code{info} gives general information about the data
    
    \item \code{OTinfo} gives general information about the Over the
    Threshold part of data. The effective duration (\code{effDuration}
    element) is the total duration for the periods with effective
    measurements.
    
    \item \code{OTdata} give OT measurements

    \item \code{OTmissing} gives start and end of the missing periods
    for OT measurements.
    
  }

  Data come from hourly sea levels measured and predicted by the french
  \emph{Service \enc{Hydrogéographique}{Hydrogeographique} et
  \enc{Océanographique}{Oceanographique} de la Marine} (SHOM).  Observed
  sea levels are available as \emph{REFMAR} data at the url
  \url{http://data.shom.fr}.  Data were processed (declustered) by IRSN
  in order to provide a series of independent surge heights at high
  tide. Surge height at high tide is defined as the difference between
  the observed and the predicted maximal sea levels near high tide.  A
  correction was applied to account for trend in the sea-level over the
  observation period.
  
  The effective duration given in years is defined up to a small
  fraction of year due to leap years and leap seconds.
  
}
\source{
  
  \url{http://data.shom.fr}

}
%%\references{~~ possibly secondary sources and usages ~~}
\examples{
str(Brest)
Brest$OTinfo$start
plot(Brest)
}
\keyword{datasets}
