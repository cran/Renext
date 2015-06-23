\name{summary.Rendata}
\alias{print.Rendata}
\alias{summary.Rendata}
\alias{print.summary.Rendata}
%%
\title{
  Summary and  print methods for "Rendata" objects
}
\description{
  Summary method for "Rendata" objects representing data
  to be used in renouvellement models.
}
\usage{

  \method{print}{Rendata}(x, \dots)
 
  \method{summary}{Rendata}(object, \dots)

  \method{print}{summary.Rendata}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{

    An object with class \code{"Rendata"}.

  }
  \item{x}{

    An object of class \code{"summary.Rendata"}, i.e.  a result of a
    call to \code{summary.Rendata}.

  }
  \item{\dots}{

    Further arguments passed to or from other methods.

  }
}
%%\details{}
%%\value{}
%%\references{}
%%\author{Yves Deville}
%%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{}
\examples{
## Brest example: no historical data
summary(Brest)

## Garonne example:  with historical data
summary(Garonne)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
