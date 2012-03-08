\name{summary.Rendata}
\alias{summary.Rendata}
\alias{print.summary.Rendata}
%%
\title{
  Summary method for "Rendata" objects
}
\description{
  Summary method for "Rendata" objects representing data
  to be used in renouvellement models.
}
\usage{
  \method{summary}{Rendata}(object, \dots)

  \method{print}{summary.Rendata}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object with class \code{"Rendata"}.}
  \item{x}{An object of class \code{"summary.Rendata"}, i.e.
    a result of a call to \code{summary.Rendata}.}
  \item{\dots}{Further arguments passed to or from other methods.}
}
%%\details{}
%%\value{}
%%\references{}
\author{Yves Deville}
%%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{}
\examples{
## no historical data
data(Brest)
summary(Brest)

## example with historical data
data(Garonne)
summary(Garonne)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
