\name{PPplot}
\encoding{UTF8}
\alias{PPplot}
\alias{QQplot}
\alias{PPplot.Renouv}
\alias{QQplot.Renouv}

\title{
  Diagnostic plots for Renouv objects 
}
\description{
  Diagnostic plots for \code{Renouv} objects.
}
\usage{
PPplot(x, ...)

QQplot(x, ...)

\S3method{PPplot}{Renouv}(x,
       showHist = FALSE,
       legend = FALSE,
       par = NULL,
       posOptions = NULL,
       ...)

\S3method{QQplot}{Renouv}(x,
       showHist = FALSE,
       legend = FALSE,
       par = NULL,
       posOptions = NULL,
       ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{

    Object containing a fitted model.

  }
  \item{legend}{

    Should a legend be shown to identify historical blocks? NOT
    IMPLEMENTED YET.
    
  }
  \item{par}{

    A list of graphical parameters as returned by \code{\link{RLpar}}
    and used to control the appearance of points. NOT IMPLEMENTED YET.
    
  }
  \item{posOptions}{

    A pair list to be passed as list of formals to the \code{\link{SandT}}
    function computing the plotting positions.
    
  }
  \item{showHist}{
    
    If \code{TRUE}, historical information contained in the object
    \code{x} (if any) will be shown using special plotting positions
    as computed by \code{\link{SandT}}. 

  }
  \item{\dots}{
    Other arguments to be passed to \code{plot}.
  }
}
\value{
  No value returned.
}
\author{
  Yves Deville
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{SandT}} for the computation of the plotting
  positions used with historical data.}
%%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
