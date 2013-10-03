\name{roundPred}
\alias{roundPred}
\title{
  Round quantiles in a pseudo-prediction table
}
\description{
  Round the quantiles of a pseudo prediction table such
  that computed by \code{predict.Renouv}.
}
\usage{
   roundPred(pred, dig.quant = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pred}{The data.frame containing the predicted quantiles and
    return levels.
  }
  \item{dig.quant}{Number of digits. Guessed if not provided.
  }
}
\details{
  Only the columns that can be considered as quantiles are
  rounded. These are assumed to have names \code{"quant"} for
  the expected return level and \code{"L."} or \code{"U."} followed by a
  percentage for lower and upper confidence limits (e.g. \code{"L.95"}
  and \code{"U.95"} for 95\% percent confidence limits.
  The number of digits guessed is experimental.
}
\value{
  A data.frame with the same structure as that given, but with some
  columns rounded.
}
%%\references{}
%%\author{}
%%\note{}
%%\seealso{}
%%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
