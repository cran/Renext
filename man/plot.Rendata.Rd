\name{plot.Rendata}
\alias{plot.Rendata}
\title{Plot a Rendata object}
\description{
  
   Plot 'Rendata' datasets with OT and historical data

 }
\usage{

    \method{plot}{Rendata}(x,
     textOver = quantile(x$OTdata[, x$info$varName], probs = 0.99),
     showHist = TRUE,
                 ...) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{

    Rendata object i.e. a list object as read with the \code{readXML}
    function.

  }
  \item{textOver}{

    Mark values of the variable in the \code{OTdata} part of
    \code{x}. Values above the \code{textOver} value (if any) will be
    marked with the character version of the block, typically a year

  }
  \item{showHist}{

    If \code{TRUE}, the historical periods (is any) are shown on the
    plot.

  }
  \item{\dots}{

    further args to be passed to \code{plot} function.

  }
}
\details{

  The plot shows the main data of the object \code{x} (the \code{OTdata}
  part) as well as historical data \code{MAXdata} or \code{OTSdata} if
  any. Different colours are used on the background. This function is
  not intended to produce nice plots to be printed.

}
\author{

  Yves Deville

}
\note{

  This function is mainly a companion function of \code{readXML}. Its
  goal is to check the content of the data read.

}
\seealso{

  \code{\link{readXML}}

}
\examples{
if (require(XML)) {
   ## use 'index.xml' file shipped with Renext
   dir1 <- system.file("Rendata", package = "Renext")
   BrestNew <- readXML(name = "Brest", dir = dir1)
   plot(BrestNew)
   GaronneNew <- readXML(name = "Garonne", dir = dir1)
   plot(GaronneNew)
   test1 <- readXML(name = "test1", dir = dir1)
   plot(test1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ }
%%\keyword{ }% __ONLY ONE__ keyword per line
