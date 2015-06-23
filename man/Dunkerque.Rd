\name{Dunkerque}
\alias{Dunkerque}
\docType{data}
\title{Surge heights at Dunkerque}
\description{
  
  Surge heights near high tide at Dunkerque tide gauge station (France)

}
\usage{Dunkerque}
\format{
  The format is:
  List of 7
  \itemize{
    \item \code{$info} : List of 6
    \itemize{
      \item \code{$name}    : chr \code{"Dunkerque"}
      \item \code{$shortLab}   : chr \code{"Surge Heights at Dunkerque (France)"}
      \item \code{$longLab}    : chr \code{"Surge Heights near high tide,  Dunkerque (France)"}
      \item \code{$varName}    : chr \code{"Surge"}
      \item \code{$varShortLab} : chr \code{"Surge"}
      \item \code{$varUnit}    : chr \code{"cm"}
    }
    \item \code{$describe} : chr \code{"High tide sea surge over 30 cm at Dunkerque... "}
    \item \code{$OTinfo} : List of 4
    \itemize{
      \item \code{$start}     : POSIXct[1:1], format: \code{"1956-01-01"}
      \item \code{$end}       : POSIXct[1:1], format: \code{"2009-01-01"}
      \item \code{$effDuration: num 38.8}
      \item \code{$threshold  : num "30"}
    }
    \item \code{$OTdata} : 'data.frame':	740 obs. of 3 variables:
    \itemize{
      \item \code{$date}  : POSIXct[1:740], format: \code{"1956-11-27" "1956-12-03"} ...
      \item \code{$Surge} : num [1:740] \code{67.9 30.9 51.8 30.8 39.8} ...
      \item \code{$comment} : chr [1:740] \code{"" "" "" ""} ...
    }
    \item  \code{$OTmissing} : 'data.frame':	83 obs. of  3 variables:
    \itemize{
      \item \code{$start}  : POSIXct[1:83], format: \code{"1956-01-01" "1956-08-08"} ...
      \item \code{$end}    : POSIXct[1:83], format: \code{"1956-06-07" "1956-11-03"} ...
      \item \code{$comment}: chr [1:83] \code{"" "" "" ""} ...
    }
    \item \code{$MAXinfo}  : 'data.frame' : 1 obs. of  3 variables:
    \itemize{
      \item \code{$start}  : POSIXct[1:1], format: \code{"1706-01-01"}
      \item \code{$end}    : POSIXct[1:1], format: \code{"1956-01-01"}
      \item \code{$duration} : num \code{250}
    }
    \item \code{$MAXdata}  :'data.frame':	1 obs. of  4 variables:
    \itemize{
      \item \code{$block}  : int \code{1}
      \item \code{$date}   : POSIXct[1:1], format: \code{"1953-02-01"}
      \item \code{$Surge}  : num \code{213}
      \item \code{$comment} : chr \code{"1"}}
    - \code{attr(*, "class")= chr "Rendata"}
  }
}
\details{
  See \code{\link{Brest}} and \code{\link{Garonne}} datasets with the
  same list structure.

  An 'historical' surge of 213 cm was observed on 1953-02-01 and is
  considered by experts as having a return period of 250 years. 
}
%\source{ }
%\references{ }
\examples{
Dunkerque$info
plot(Dunkerque)
}
\keyword{datasets}
