\name{predict.Renouv}
\alias{predict.Renouv}
\title{
  Compute return levels and confidence limits for a "Renouv" object.
}
\description{
  Compute return levels and confidence limits for an object of class "Renouv".
 }
\usage{
   \method{predict}{Renouv}(object,
        newdata = c(10, 20, 50, 100, 200, 500, 1000),
        cov.rate = FALSE,
        level = c(0.95, 0.7),
        trace = 1, eps = 1e-06,
        \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{"Renouv"} typically created
    by using the \code{Renouv} function.
  }
  \item{newdata}{The return period at which return levels and confidence
    bounds are wanted.
  }
  \item{cov.rate}{If \code{FALSE}, the delta method will not take
    into account the uncertainty on the event rate \code{lambda}
    of the Poisson process. NOT IMPLEMENTED YET. At the time only
    \code{FALSE} is possible.
  }
  \item{level}{Confidence levels as in other 'predict' methods (not
    percentages).
  }
  \item{trace}{Some details are printed when \code{trace} is not zero.
  }
  \item{eps}{Level of perturbation used to compute the numerical
    derivatives in the delta method.
  }
  \item{\dots}{Further arguments passed to or from other methods.
  }
}
\details{
  Unless in some very special cases, the confidence limits are
  approximated ones computed by using the delta method with numerical
  derivatives.
   
 }
\value{
  A data.frame with the expected return levels (col. named
  \code{"quant"}) at the given return
  periods, and confidence limits. The returned object has an
  \code{infer.method} attribute describing the method used to compute
  the confidence limits.
  
}
\references{
   Coles S. (2001)  \emph{Introduction to Statistical Modelling
     of Extremes Values}, Springer.
}
\author{
  Yves Deville
}
\note{
  Despite of its name, this method does not compute true
  predictions. A return period is to be interpreted as
  an average interevent time rather than the duration of a specific
  period of time. For instance, the expected return level for a given
  return period with length 100 years is the level
  that would be on average exceeded once every 100 years
  (assuming that the model description in \code{object} is correct).
}
\seealso{
  \code{\link{Renouv}} to fit \code{Renouv} model.
}
\examples{
library(Renext)
data(Brest)
fit <- Renouv(Brest)
pred <- predict(fit, newdata = c(100, 125, 150, 175, 200),
                level = c(0.99, 0.95))


}
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
