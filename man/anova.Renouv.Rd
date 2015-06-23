\name{anova.Renouv}
\alias{anova.Renouv}
\title{
  Compute an analysis of deviance table for two nested Renouv
  objects
}
\usage{
   \method{anova}{Renouv}(object, object1, trace = 1L, ...)
}
\arguments{
  
  \item{object}{

    A \code{Renouv} model as fitted with \code{\link{Renouv}}.

  }
  \item{object1}{

    A \code{Renouv} object such that \code{object} is nested in
    \code{object1}.

  }
  \item{trace}{

    Level of verbosity. The value \code{0} prints nothing.

  }
  \item{...}{

    Not used yet.

  }
}
\value{

  An object of class \code{"anova"} inheriting from class \code{"data.frame"}.

}
\description{
  
Compute an analysis of deviance table for two nested Renouv objects

}
\seealso{
  \code{\link{anova}}, \code{\link{LRExp.test}}.
}
\details{
  
Of special interest is the case when the distribution of the excesses
used in \code{object} is exponential while \code{object1} uses a
two-parameters alternative in the GPD family. We know then that the
convergence to the asymptotic distribution is slow, and a numerical
approximation of the true distribution of the test statistic is used
when possible, i.e. when the objects do not use MAX or OTS data and the
number of exceedances is between 8 and 500.

}
\note{
  
The deviance of the models can not be interpreted: only the
difference of the deviance is meaningful.

}
\examples{
## test using historical data
fit1Exp <- Renouv(Garonne,  distname.y = "exponential", plot = FALSE)
fit1GPD <- Renouv(Garonne, distname.y = "GPD", plot = FALSE)
anova(fit1Exp, fit1GPD)

## test without using historical data
x <- Garonne$OTdata$Flow
dur <- Garonne$OTinfo$effDuration

fit2Exp <- Renouv(x,  threshold = 2700,  effDuration = dur,
                  distname.y = "exponential", plot = FALSE)
fit2GPD <- Renouv(x, threshold = 2700, effDuration = dur,
                  distname.y = "GPD", plot = FALSE)
anova(fit2Exp, fit2GPD)
}
