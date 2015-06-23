\name{spacings}
\encoding{UTF8}
\alias{spacings}
\alias{spacings.numeric}
\alias{spacings.data.frame}
\alias{spacings.Rendata}

\title{

  Methods computing spacings between Largest Order Statistics

}
\description{
  
  Methods computing the random \emph{spacings} for the Largest Order
  Statistics of the marks in a POT renewal.

}
\usage{
spacings(object, ...)

\S3method{spacings}{numeric}(object, wExp = TRUE, ...)

\S3method{spacings}{data.frame}(object, varName, wExp = TRUE, ...)

\S3method{spacings}{Rendata}(object, type = c("MAX", "OTS", "OT"), wExp = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{

    A object containing the marks \eqn{X_i}. This is can be a vector or
    a data frame for data aggregated by blocks. For an object of class
    \code{data.frame} one of the three data.frame slots: \code{OTdata},
    \code{MAXdata} or \code{OTSdata} can be used using the suitable
    value of \code{type}.

  }
  \item{varName}{

    Character vector of length 1 giving the name of the variable when
    \code{object} is a data.frame.

  }
  \item{wExp}{
    
    Logical. If \code{TRUE}, the spacings are weighted as explained in
    \bold{Details}.

  }
  \item{type}{

    Character specifying the data.frame to be used when object has class
    \code{"Rendata"}.

  }
  \item{\dots}{
    
    Not used yet.

  }
}
\details{

  The spacings are the differences between the Largest Order
  Statistics. They are useful for some estimation tasks or diagnostics.
  Given random variables \eqn{X_i}, the \eqn{i}-th spacing \eqn{Y_i} is
  the difference \eqn{X_{(i)}-X_{(i+1)}} between the \eqn{i}-th largest
  order statistic \eqn{X_{(i)}} and the next in the decreasing order
  i.e. \eqn{X_{(i+1)}}.

  When the r.vs \eqn{X_i} form a sample of an exponential or Gumbel
  distribution, the spacings associated with the largest order
  statistics are or tend to be independent and exponentially
  distributed. More precisely, the weighted spacings \eqn{i \times
  Y_i}{i * Y_i} have or tend to have the same exponential
  distribution. This can be used to estimate the shape parameter of the
  underlying distribution using only the largest order
  statistics. Moreover the \eqn{r-1} spacings \eqn{Y_i} built from the
  \eqn{r} largest order statistics \eqn{i \le r} are or tend to be
  independent from the \eqn{r}-th largest order statistic \eqn{X_{(r)}}.

  When \code{wExp} is \code{TRUE}, the returned values are the weighted
  spacings \eqn{i \times Y_i}{i * Y_i}.
  
}
\value{
  
  A list or vector containing the spacings. When the data is structured
  in blocks as is the \code{MAXdata} slot of an object of class
  \code{"Rendata"}, the spacings are computed form the order statistics
  \emph{within each block}, to maintain independence with the next order
  statistic in decreasing order.

}
\references{

  Embrechts P., Klüppelberg C. and Mikosch T. (1997) \emph{Modelling
  Extremal Events for Insurance and Finance}. Springer.

}
\author{
  Yves Deville
}
%%\note{}
\section{Caution}{
  
  By default, the spacings are scaled as explained above, thus assuming
  that the marks are exponentially distributed.

}
\examples{
sp <- spacings(rgumbel(200, loc = 0, scale = 1))
expplot(sp)
sp1 <- spacings(rgev(200, loc = 0, scale = 1, shape = 0.3))
expplot(sp1)
## spacings are computed by block
sp2 <- spacings(object = Garonne$MAXdata,
                varName = Garonne$info$varName)
expplot(unlist(sp2))
sp3 <- spacings(object = Garonne, type = "OT")
expplot(sp3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
