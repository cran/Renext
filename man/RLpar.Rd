\name{RLpar}
\alias{RLpar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Graphical parameters for Return Level plots
}
\description{
  Build a hierarchical list of graphical parameters that can
  be used in the methods plot or lines for the class \code{"Renouv"}. 
}
\usage{
RLpar(mono = TRUE,
      trace = 0L,
      ...)
}
\arguments{
  \item{mono}{

    Logical. The default \code{TRUE} is for plots possibly using colors
    but that can be printed in grayscale.  With the value \code{FALSE},
    curves or symbols will appear distinctly on a color device but not
    necessarily when printed in grayscale.

  }
  \item{trace}{

    Integer level of verbosity. The default value \code{0} prints
    nothing.

  }
  \item{\dots}{

    Arguments with names corresponding to the hierarchical structure
    and the graphical parameter to be changed.

  }
}
\details{
  
  The formals are in correspondence with the list hierarchy using a
  column \code{"."} as separator to define the tree. Thus a
  \code{quant.col} formal argument can be used to specify the color of
  the quantile (or return level) curve, while \code{conf.conf1.col} will
  be used for the first confidence limits (lower and upper).

}
\value{
  
  A list containing lists in a hierarchical fashion. At the root level,
  an element concerns a single curve (e.g. the return level curve), a
  single scatterplot (e.g. sample used in POT), a group of curves
  (e.g. the confidence limits) or a group of scatterplots (e.g. the
  collection of \code{MAX} historical blocks). For single elements
  (curve or scatterplot) the list contains graphical elements with
  values as they would be given in \code{plot} or \code{lines}
  calls. For group elements, each element is a list of such lists.

}
%%\references{}
\author{
 Yves Deville
}
\note{
  
  A list of default parameter values is built first using the model
  suitable for the \code{mono} value. Then the values provided by the
  user overwrite the existing. Thus a curve can be coloured even if
  \code{mono = TRUE}, if a colour specification is given for the
  corresponding element.

  When the same parameter name is used several times in \code{RLpar}, a
  warning is thrown.
  
}
\seealso{

  \code{\link{plot.Renouv}} and \code{\link{lines.Renouv}} with which
  \code{RLpar} is to be used.

}
\examples{
## change color for quantile curve and type for confidence
## limits #1 (with largest confidence level).
newRLpar <- RLpar(quant.col = "red", conf.conf1.lty = "dashed")
newRLpar$quant

## show the names of all possible editable parameters
names(unlist(RLpar()))



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
