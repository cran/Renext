\name{Ren2gev}
\alias{Ren2gev}
\title{

  Translate a vector of coefficients from a Renewal-POT model with
  Pareto excesses into a vector of GEV parameters

}
\description{
  
  Translate a vector of coefficients from a Renewal-POT model with
  Pareto excesses into a vector of GEV parameters.

}
\usage{
Ren2gev(object,
        threshold = NULL,
        w = 1,
        distname.y = c("gpd", "GPD", "lomax", "maxlo"),
        jacobian = (length(w) == 1L),
        vcovRen = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{

    A named vector of parameters or an object of class \code{"Renouv"}.
    In the first case, the names of the vector element must conform to
    the distribution given in \code{distname.y}.

  }
  \item{threshold}{

    The threshold associated with the renewal-POT model.  This must be
    provided and be a non NA finite numeric value. It is the location
    parameter of the GPD.

  }
  \item{w}{
    
    The duration of the blocks.

  }
  \item{distname.y}{

    The distribution of the excesses in the renewal-POT model.  This
    is normally a \code{"gpd"} but can be a \code{"lomax"} or a
    \code{"maxlo"} distribution provided that the GEV parameters given
    in \code{object} specify a positive or a negative shape
    respectively.

  }
  \item{jacobian}{

    Logical. If \code{TRUE}, the jacobian matrix of the transformation
    is computed. This is only possible at the time when \code{w} has
    length 1.

  }
  \item{vcovRen}{

    A covariance matrix for the "Ren" vector of parameters.  If
    \code{object} has class \code{"Renouv"}, then the covariance matrix
    embedded in the object is used.

  }
  
}
\details{
  
  Given Renewal-POT parameters, it is possible to compute the
  distribution of block maxima. When the distribution is in the Pareto
  family, the marginal distribution of maxima is GEV. The location and
  the scale GEV parameters depend on the block duration \eqn{w}, while
  the GEV shape parameter is identical to that of the GPD input
  distribution.

}
\value{
  
  When \code{w} has length \code{1}, a named vector of GEV parameters as
  the one estimated by \code{\link[evd]{fgev}}. This vector has an
  elements named \code{"loc"}, \code{"scale"} and \code{"shape"}.
  
  When \code{w} has length \code{> 1}, a matrix with \code{length(w)}
  rows, each representing a vector of GEV parameters as before.
  
  The returned object has attributes named \code{"threshold"}.  and
  \code{"distname.y"} to recall how it was built.
  
}
%%\references{}
\author{
  
  Yves Deville

}
%%\note{}
\seealso{

  The \code{\link{gev2Ren}} function provides a reciprocal
  transformation.

}
\examples{
fit1 <- Renouv(Garonne, distname.y = "maxlo")
Ren2gev(fit1)
fit2 <- Renouv(Garonne, distname.y = "gpd")
Ren2gev(fit2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
