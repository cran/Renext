\name{gumbel2Ren}
\alias{gumbel2Ren}
\title{
  Translate a vector of Gumbel parameters into a vector of parameters
  for a renewal model
}
\description{
  
  Translate a (named) vector of Gumbel parameters into a vector of
  parameters for a renewal model.

}
\usage{
gumbel2Ren(parGumbel,
           threshold = NULL,
           lambda = NULL,
           w = 1,
           distname.y = c("exponential"),
           vcovGumbel = NULL,
           jacobian = TRUE,
           plot = FALSE)
}
\arguments{
  \item{parGumbel}{

    Named vector of Gumbel parameters, with name \code{"scale"} and
    \code{"shape"}.

  }
  \item{threshold}{

    The threshold for the target Renouv parameters.

  }
  \item{lambda}{
    
    The rate for the target Renouv parameters.

  }
  \item{w}{

    A block duration for the conversion.

  }
  \item{distname.y}{

    The distribution of the excesses.

  }
  \item{vcovGumbel}{

    A covariance matrix for the Gumbel parameters.

  }
  \item{jacobian}{

    Logical. If \code{TRUE}, the jacobian is used.

  }
  \item{plot}{
    
    Logical. If \code{TRUE}, a rough plot will be drawn.

  }
}
\details{
  
  Given a vector of Gumbel parameters and a block duration, there exits
  an infinity of \code{Renouv} models with exponential excesses leading to
  the prescribed Gumbel distributions for the maximum of the marks on a
  block with duration \code{w}.  One of these models may be chosen by
  specifying either a threshold or a rate \code{lambda}.

}
\value{
  
  A vector of Renouv parameters, which can be used with
  \code{\link{RenouvNoEst}}.

}
%%\references{}
\author{
  Yves Deville
}
\section{Caution}{
  
  All Renouv models lead to the same return level curve whatever be the
  choice of \code{threshold} or \code{lambda}. However, when a
  covariance matrix is given, the covariance matrix for the Renouv
  parameters and consequently the confidence bounds \bold{depend on the
  threshold or rate}. So the computed covariance matrix must in general
  be considered as putative.

}
%%\note{}
\seealso{
  
  \code{\link{gev2Ren}} for a similar translation from the GEV
  distribution.

}
%%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
