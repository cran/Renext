\name{Renext-package}
\alias{Renext-package}
\alias{Renext}
\docType{package}
\title{
Renewal method for extreme values extrapolation
}
\description{
   This package proposes fits and diagnostics for the so-called "renewal method", an
   alternative to other "Peaks Over Threshold" (POT) methods.
   The renewal method generalises the classical POT by allowing the exceedances
   over the threshold to follow a probability distribution which can differ from the
   Generalised Pareto Distribution (GPD). Weibull or gamma exceedances
   are sometimes preferred to GPD exceedances. The special case of
   exponential exceedances (which falls in the three families: GPD,
   Weibull and gamma) has a special interest since it allows exact
   inference for the (one dimensional) parameter and for the quantiles
   form OT data (only).
   
   The package allows the joint use of possibly three kinds of data or
   information. The first kind is \emph{classical exceedances}, or
   \emph{"OT data"} and will always be required. It can be completed
   with a second kind of data called \emph{historical data} consisting
   in the largest observations over one or several periods called
   historical blocks. A third usable kind of information called
   \emph{unobserved levels} occurs when a fixed level is known not to
   have been exceeded during a period with known duration. Historical
   data are often available in hydrology (e.g. river flows) for large
   periods such as past centuries. An unobserved level can typically be
   related to a material benchmark.

   Maximum likelihood estimation is made possible in this context of
   heterogeneous data. Inference is based on the asymptotic normality of
   parameter vector estimate and on linearisation ("delta method") for
   quantiles or parameter functions.

   The package allows the use of "marked-process observations" data
   (datetime of event and level) where an interevent analysis can be
   useful. It also allows that the event dates are unknown and replaced
   by a much broader \emph{block} indication, e.g. a year number. The
   key point is then that the "effective duration" (total duration of
   observation periods) is known. Event counts for blocks can be used to
   check the assumption of Poisson-distributed events.
   
   The package development was initiated, directed and financed by the french \emph{Institut de
   radioprotection et de surete nucleaire} (IRSN). The package is a
   non-academic tool designed for applied analysis on case studies and
   investigations or comparisons on classical probabilistic models.
   
   Additional information and material related to this package can be
   found at the URL
   \url{https://gforge.irsn.fr/gf/project/renext}.
   
}
\details{
  \tabular{ll}{
    Package: \tab Renext\cr
    Type: \tab Package\cr
    Version: \tab 1.0-0\cr
    Date: \tab 2010-06-10\cr
    License: \tab GPL (>=2)\cr
    LazyLoad: \tab yes\cr
  }
  This package contains functions \code{\link{fRenouv}} and \code{\link{rRenouv}} to fit
  or simulate from "renouvellement" models.
}
\author{
  Yves Deville <deville.yves@alpestat.com>, IRSN <renext@irsn.fr>
  
Maintainer: IRSN <renext@irsn.fr>
}
\references{
  \itemize{
    \item Miquel J.(1984) \emph{Guide pratique d'estimation des
      probabilites de crues}, Eyrolles (coll. EDF DER)
    \item Coles S. (2001)  \emph{Introduction to Statistical Modelling
      of Extremes Values}, Springer.
  }
}
\keyword{datagen}
\seealso{
  The packages \code{evd},
  \code{ismev},
  \code{extRemes},
  \code{bayesevd},
  \code{POT}.
}
\examples{
## artificial data
test <- rRenouv(nb = 100, threshold = 40,
                par.N = list(lambda = 2),
                densfun.y = "weibull", par.y = mom2par("weibull", mean = 30, sd = 36))  

## exponential exceedances
fit.expon <- fRenouv(x.OT = test$x,
                     sumw.BOT = 100,
                     distname.y = "exponential",
                     threshold = 50,
                     conf.pct = c(70, 95, 99),
                     main = "\"exponential\" distr. (special)")
}
