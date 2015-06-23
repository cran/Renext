\name{Garonne}
\Rdversion{1.1}
\alias{Garonne}
\encoding{UTF-8}
\docType{data}
\title{
  Flow of the french river La Garonne
}
\description{
  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
  Flow of the french river La Garonne at le Mas d'Agenais
}
\usage{Garonne}
\format{
  The format is:
  List of 7
  \itemize{
    \item \code{$info} : List of 6
    \itemize{
      \item \code{$name}       : chr \code{"Garonne"}
      \item \code{$shortLab}   : chr \code{"La Garonne at Le Mas d'Agenais"}
      \item \code{$longLab}    : chr \code{"River flow of La Garonne at Le Mas d'Agenais"}
      \item \code{$varName}    : chr \code{"Flow"}
      \item \code{$varShortLab} : chr \code{"Flow"}
      \item \code{$varUnit}    : chr \code{"m3/s"}
    }
    \item \code{$describe} : chr \code{"Flow of the french river La Garonne ..."}
    \item \code{$OTinfo} :List of 4
    \itemize{
      \item \code{$start}       : POSIXct[1:1], format: \code{"1913-01-01"}
      \item \code{$end}         : POSIXct[1:1], format: \code{"1978-01-01"}
      \item \code{$effduration} : num \code{65}
      \item \code{$threshold}   : num \code{2500}
    }
    \item \code{$OTdata} : 'data.frame':	151 obs. of  3 variables:
    \itemize{
      \item \code{$date} : POSIXct[1:151], format: \code{"1913-04-08" "1913-04-25"}  ...
      \item \code{$Flow}  : num [1:151] \code{2600 2800 2700 4579 3400}
   ...
      \item \code{comment} : chr [1:151] "" "" "" "" ...
    }
    \item \code{$OTmissing} : \code{NULL}
    \item \code{$MAXinfo} :'data.frame':	1 obs. of  3 variables:
    \itemize{
      \item \code{$start}    : POSIXct[1:1], format: \code{"1770-01-01"}
      \item \code{$end}      : POSIXct[1:1], format: \code{"1913-01-01"}
      \item \code{$duration} : num \code{143}
    }
    \item \code{$MAXdata} :'data.frame':	12 obs. of 4 variables:
    \itemize{
      \item \code{$block}  : num [1:12] \code{1 1 1 1 1 1 1 1 1 1} ...
      \item \code{date} : POSIXct[1:12], format: \code{NA NA} ...
      \item \code{$Flow}      : num [1:12]
      \code{7500 7400 7000 7000 7000 6600 6500 6500 6400 6300} ...
      \item \code{$comment} : chr [1:12] "1 (1875)" "2 (1770)" "3 (1783)" "4 (1855)" ...
    }
    - \code{attr(*, "class")= chr "Rendata"}
    }
    }
\details{
   The data concern the french river \emph{La Garonne} at the gauging
   station named \emph{Le Mas d'Agenais} where many floods occurred during 
   the past centuries.
   
   The data consist in OT data and historical data. The variable is
   the river flow in cube meter per second \eqn{(\textrm{m}^3/\textrm{s})}{(m3/s)} as estimated from 
   the river level using a rating curve. The precision is limited and many ties are present
   among the flow values.

   The OT data or "OTdata" contain flows values over the threshold
   \eqn{u = 2500\,\mathrm{m}}{u = 2500 m} for the \eqn{65}{65} years
   period 1913-1977.
   The historical data or "MAXdata"  is simply the \eqn{r=12}{r=12} largest flows for the period 
   of \eqn{143}{143} years 1770-1912. The exact dates of these events are not known with
   precision but the years are known and given as comments.
  
}
\source{
   The data were taken from the book by Miquel.
}
\references{
  Miquel J. (1984) \emph{Guide pratique d'estimation des
    \enc{probabilités}{probabilites} de crues}, Eyrolles (coll. EDF DER).

  Parent E. and Bernier J. (2003) Bayesian POT modeling for Historical
  data. \emph{Journal of Hydrology} vol. 274, pp. 95-108.
  

   
}
\examples{
plot(Garonne) 
}
\keyword{datasets}
