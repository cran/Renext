\name{readXML}
\alias{readXML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read data using an XML index file}
\description{
  Read one or several dataset(s) using an XML index file specifying the
  data sets to read and the structure of each
}
\usage{
   readXML(name,
           dir,
           index = "index.xml",
           trace = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{Name for the dataset that will be matched against the
    \code{name} attribute of datasets as they are given in the index file.}
  \item{dir}{Path to the directory where the index file and all data
    files should be found.}
  \item{index}{Name (short) of the index file. This file must be in the
    directory given by \code{dir}.}
  \item{trace}{Integer or logical to trace the successive steps by short
  indications.}
}
\details{
   The XML index file is parsed within R. Then according to the
   indications within the index file, other files are read (e.g. csv
   files). In this way, data returned as a list can contain heterogeneous data:
   Over Threshold (OT) data, missing periods, MAX data, etc.
   Various informations are also stored in list elements with name
   containing the \code{"info"} string.

   This function requires the CRAN package XML.
}
\value{
  A list with the data read.
  \item{info}{General information about the data: \code{varName},
    \code{varShorlLab} and \code{varUnit} give the variable name
    unit and short label.}
  \item{OTinfo}{Information for the Over the Threshold (OT).}
  \item{OTdata}{Over the Threshold (OT) data.}
  \item{OTmissing}{Missing periods within the OTdata period.}
  \item{MAXinfo}{Information for the MAX (\eqn{r}{r}-largest) supplement data.}
  \item{MAXdata}{MAX supplement data.}
  \item{OTSinfo}{Information for the Over the Threshold Supplement (OTS)
    data.}
  \item{OTSdata}{Over the Threshold (OT) supplement data.}
}
%\references{}
\author{Yves Deville}
\note{
  The flat files (usually \code{.csv} files) can also be read in a more
  conventional way e.g. through \code{read.table}. However, conform to
  the \code{index.xml} examples or to the \code{index.xsd} schema to see
  how to adjust the reading of parameters such as \code{sep}, etc. 
}
\seealso{
  See \code{\link{Brest}} for an example of such a list.
  
  See URL \url{https://gforge.irsn.fr/gf/project/renext} for data
  related to \code{Renext}.
  
}
\examples{
\dontrun{
library(XML)
## use 'index.xml' file shiped with Renext
dir1 <- system.file("Rendata", package = "Renext")
BrestNew1 <- readXML(name = "Brest", dir = dir1)
test1 <- readXML(name = "test1", dir = dir1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
