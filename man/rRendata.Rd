\name{rRendata}
\alias{rRendata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

  Simulate a random RenData object

}
\description{
  
  Simulate a random \code{RenData} object that can be used within the
  \code{Renouv} function for tests.

}
\usage{
rRendata(lambda = 1,
         threshold = 0,
         effDuration = 100,
         distname.y = "exp",
         par.y = c(rate = 1),
         start = "1913-01-01",
         name = NULL,
         varName = "X", varUnit = "?",
         simDate = TRUE, roundDate = FALSE,
         MAX.effDuration = NULL,
         MAX.r = rep(1L, length(MAX.effDuration)),
         OTS.effDuration = NULL,
         OTS.threshold = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{

    The rate of the Homogeneous Poisson Process.

  }
  \item{threshold}{

    The threshold for the exceedances.

  }
  \item{effDuration}{

    The effective duration of the main Over Threshold (OT) period.
    This must be a positive value.

  }
  \item{distname.y}{

    Name of the distribution for the excesses to be simulated.  See
    \bold{Details}.

  }
  \item{par.y}{

    A named vector or list giving the parameters values for the
    distribution. The name must conform to the chosen distribution.

  }
  \item{start}{

    A \code{POSIXct} object, or character that can be coerced to
    \code{POSIXct} (e.g. a date given as a character in the
    \code{"YYYY-MM-DD"} format) giving the start of the main OT sample.

  }
  \item{name}{
    
    A name for the dataset which will be attached to it and be used by
    some methods for \code{"Rendata"}.

  }
  \item{varName}{

    Name of the simulated variable.

  }
  \item{varUnit}{

    Unit for the simulated variable (is used by plot). 

  }
  \item{simDate}{

    Logical. If \code{TRUE} the dates will be reported for the
    historical data (MAX and OTS).

  }
  \item{roundDate}{

    Logical. If \code{TRUE} the time part ot the \code{date} column will
    be rounded. Not implemented yet.

  }
  \item{MAX.effDuration}{
    
    Vector of the durations for the \code{MAX} historical blocks.

  }
  \item{MAX.r}{

    Vector of the (positive) numbers of observations for \code{MAX}
    historical blocks. Must have the same length as
    \code{MAX.effDuration}. See \bold{Caution} below for the effect of
    selection large values.

  }
  \item{OTS.effDuration}{
    
    Vector of durations for the \code{OTS} historical blocks.

  }
  \item{OTS.threshold}{

    Vector of numerical thresholds for the observations in \code{OTS}
    historical blocks. Must have the same length as
    \code{OTS.effDuration}. All values must be \code{>= threshold}.

  }
}
\details{
  
  The distribution of the excesses named in \code{distname.y} can be
  any known distribution, provided that when prefixed with the usual
  \code{"r"} letter, the name gives the wanted simulation function. For
  example, with \code{distname.y = "exp"}, the \code{rexp} function is
  used and \code{par.y} must thus contain an element with name
  \code{"rate"}.
  
  When a suitable numeric threshold is given, the simulated marks of the
  marked process are the sum of the threshold and of a random excess
  drawn from \code{distname.y}.  When the threshold is not a finite
  numeric value, the observed marks are the simulated values themselves.

  The main OT sample is assumed to begin at \code{start}. Historical MAX
  blocks (if any) are assumed to be just before \code{start}, and OTS
  are just before \code{start} or just before the beginning of the MAX
  blocks when there are some.  The dates are computed without taking
  into consideration the problems of leap years or leap seconds.

}
\value{

  An object with S3 class \code{"Rendata"}. This class currently has
  \code{plot} and \code{summary} methods.

}
%%\references{}
\author{

  Yves Deville

}
\note{

  When \code{effDuration} is small relative to the inverse of
  \code{lambda} the number of simulated marks in the OT sample may be
  \eqn{0} which can cause problems for some uses of the created data.

}
\section{Caution}{
  
  By construction, each \code{MAX} block contains at least one
  observation, while a random period of the same duration might have
  none. The simulated number of events on a MAX block is generated using
  a censored Poisson distribution.  Care must be taken when estimations
  are made on such data, since creating \code{MAX} blocks obviously
  create a positive bias on \code{lambda}. Such bias then also affects
  the other parameters concerning the excesses, because these
  parameters are no longer orthogonal to the rate parameter
  \code{lambda} when historical data are used. The bias can be severe if
  \code{MAX} blocks with small durations are used, or if large number of
  events are chosen in \code{MAX.r}.

}
\seealso{
  
  \code{\link{plot.Rendata}}, \code{\link{summary.Rendata}}.

}
\examples{
set.seed(1234)
rd <- rRendata(effDuration = 60,
               MAX.effDuration = rep(3, 6),
               MAX.r = rep(4, 6),
               distname.y = "exp", par.y = c(rate = 1/100))
plot(rd)
summary(rd)
rd2 <- rRendata(effDuration = 10,
                MAX.effDuration = rep(60, 2),
                MAX.r = rep(3, 2),
                simDate = FALSE,
                distname.y = "gpd", par.y = c(scale = 20, shape = 0.16))
plot(rd2)
rd3 <- rRendata(effDuration = 10,
                OTS.effDuration = rep(60, 2),
                OTS.threshold = rep(80, 2),
                simDate = FALSE,
                distname.y = "gpd", par.y = c(scale = 20, shape = 0.16))
plot(rd3)
## Renouv fit with historical data
fit <- Renouv(rd)
summary(fit)
}

%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
