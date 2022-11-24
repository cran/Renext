\name{Renext-package}
\alias{Renext-package}
\alias{Renext}
\docType{package}
\encoding{UTF-8}
\title{
  
  Renewal Method for Extreme Values Extrapolation

}
\description{

  This package proposes fits and diagnostics for the so-called
  \emph{méthode du renouvellement}, an alternative to other "Peaks Over
  Threshold" (POT) methods.  The \emph{méthode du renouvellement}
  generalises the classical POT by allowing the excesses over the
  threshold to follow a probability distribution which can differ from
  the Generalised Pareto Distribution (GPD). Weibull or gamma
  excesses are sometimes preferred to GPD excesses. The special
  case of exponential excesses (which falls in the three families:
  GPD, Weibull and gamma) has a special interest since it allows exact
  inference for the (scalar) parameter and for the quantiles form OT
  data (only).
   
  The package allows the joint use of possibly three kinds of data or
  information. The first kind is \emph{classical excesses}, or
  \emph{"OT data"}. It can be completed with two kinds of data
  resulting from a temporal aggregation, as is often the case for
  \emph{historical data}. Both types are optional, and concern periods
  or \emph{blocks} that must not overlap nor cross the OT period.
  
  \itemize{
    
    \item \emph{MAX data} correspond to the case where one knows the
    \eqn{r}{r} largest observations over each block. The number
    \eqn{r}{r} may vary across blocks. This kind of data is often
    called '\eqn{r}{r} largest', or "\eqn{r}{r} Largest Order
    Statistics" (\eqn{r}{r} LOS).
    
    \item \emph{OTS data} (for OT Supplementary data) correspond to the
    case where one knows for each block \eqn{b}{b} all the observations
    that exceeded a threshold \eqn{u_b}{u[b]} which is greater (usually
    much greater) than the main threshold \eqn{u}{u}. The number
    \eqn{r_b}{r[b]} of such observations can be zero, in which case we
    may say that \eqn{u_b}{u[b]} is an \emph{unobserved level}. A
    threshold \eqn{u_b}{u[b]} is sometimes called a \emph{perception
      threshold}.
    
  }
  
  Historical data are often available in hydrology (e.g. for river flood
  discharges, for sea-levels or sea surges) and can concern large
  periods such as past centuries. An unobserved level can typically be
  related to a material benchmark.
  
  Maximum likelihood estimation is made possible in this context of
  heterogeneous data. Inference is based on the asymptotic normality of
  parameter vector estimate and on linearisation ("delta method") for
  quantiles or parameter functions.
  
  The package allows the use of "marked-process observations" data
  (datetime of event and level) where an interevent analysis can be
  useful. It also allows the event dates to be unknown and replaced
  by a much broader \emph{block} indication, e.g. a year number. The
  key point is then that the "effective duration" (total duration of
  observation periods) is known. Event counts for blocks can be used to
  check the assumption of Poisson-distributed events.
  
  The package development was initiated, directed and financed by the
  french \emph{Institut de Radioprotection et de \enc{Sûreté}{Surete}
    \enc{Nucléaire}{Nucleaire}} (IRSN). The package is a non-academic
  tool designed for applied analysis on case studies and investigations
  or comparisons on classical probabilistic models.
  
  %% Additional information and material related to this package can be
  %% found at the URL \url{https://gforge.irsn.fr/gf/project/renext}.
  
}
\details{
  
  %% \tabular{ll}{
  %%   Package: \tab Renext\cr
  %%   Type: \tab Package\cr
  %%   Version: \tab 3.1-0\cr
  %%   Date: \tab 2016-07-01\cr
  %%   License: \tab GPL (>= 2)\cr
  %%   LazyLoad: \tab yes\cr

  %% }
  The DESCRIPTION file:
  \packageDESCRIPTION{Renext}
  \packageIndices{Renext}

  This package contains a function \code{\link{Renouv}} to fit
  "renouvellement" models.

}
\author{
  
  \packageAuthor{Renext}
  %% Yves Deville <deville.yves@alpestat.com>, IRSN <renext@irsn.fr>

  Maintainer: \packageMaintainer{Renext}
  %% Maintainer: Lise Bardet <lise.bardet@irsn.fr>

}
\references{

  \itemize{

    \item Miquel J. (1984) \emph{Guide pratique d'estimation des
    \enc{probabilités}{probabilites} de crues}, Eyrolles (coll. EDF DER).

    \item Coles S. (2001) \emph{Introduction to Statistical Modelling of
      Extremes Values}, Springer.

    \item Embrechts P., Klüppelberg C. and Mikosch T. (1997)
    \emph{Modelling Extremal Events for Insurance and
    Finance}. Springer.

  }
}
\keyword{datagen}
\seealso{
  The CRAN packages \CRANpkg{evd}, \CRANpkg{ismev},
  \CRANpkg{extRemes}, \CRANpkg{POT}.
}
\examples{
## 'Garonne' data set
summary(Garonne)
plot(Garonne)

## Weibull excesses
fG <- Renouv(x = Garonne,
             threshold = 3000,
             distname.y = "weibull",
             main = "Weibull fit for 'Garonne'")

coef(fG)
vcov(fG)
summary(fG)
logLik(fG)
## Re-plot if needed
plot(fG)

## Classical 'predict' method with usual formal args 
predict(fG, newdata = c(100, 150, 200), level = c(0.8, 0.9))
}
