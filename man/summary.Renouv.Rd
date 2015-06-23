\name{summary.Renouv}
%%
\alias{print.Renouv}
\alias{summary.Renouv}
\alias{print.summary.Renouv}
\alias{format.summary.Renouv}
%%
\title{Summary and print methods for "Renouv" objects}
%%
\description{
  Summary method for "Renouv" objects representing 'Renouvellement' (POT)
  fitted models.
}
\usage{
   \method{print}{Renouv}(x,
        digits = max(3L, getOption("digits") - 3L),
        \dots)

   \method{summary}{Renouv}(object,
        correlation = FALSE,
        symbolic.cor = FALSE,
        \dots)

   \method{print}{summary.Renouv}(x,
      coef = TRUE,
      pred = TRUE,
      probT = FALSE,
      digits = max(3, getOption("digits") - 3),
      symbolic.cor = x$symbolic.cor,
      signif.stars = getOption("show.signif.stars"),
      \dots)

   \method{format}{summary.Renouv}(x,
      \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{

    An object with class \code{"Renouv"}.

  }
  \item{x}{

    An object of class \code{"summary.Renouv"}, i.e.  a result of a call
    to \code{summary.Renouv}.

  }
  \item{correlation}{

    Logical; if \code{TRUE}, the correlation matrix of the estimated
    parameters is returned and printed.

  }
  \item{coef}{

    Logical. If \code{FALSE}, the table of coefficients and t-ratios'
    will not be printed.

  }
  \item{pred}{

    Logical. If \code{FALSE}, the table of return periods/levels will
    not be printed.

  }
  \item{probT}{

    If \code{FALSE}, the \eqn{p}{p}-values for the t-tests will not be
    printed nor displayed.

  }
  \item{digits}{

    the number of significant digits to use when printing.

  }
  \item{symbolic.cor}{

    logical. If \code{TRUE}, print the correlations in a symbolic form
    (see \code{\link{symnum}}) rather than as numbers.

  }
  \item{signif.stars}{

    logical. If \code{TRUE}, \sQuote{significance stars} are printed for
    each coefficient.

  }
  \item{\dots}{

    Further arguments passed to or from other methods.

  }
}
\details{
  
  \code{print.summary.Renouv} tries to be smart about formatting the
  coefficients, standard errors, return levels, etc.
  \code{format.summary.Renouv} returns as a limited content as a
  character string. It does not embed coefficients values nor
  predictions.

}
\value{
  
 The function \code{summary.RenOUV} computes and returns a list of
 summary statistics concerning the object of class \code{"Rendata"}
 given in \code{object}.  The returned list is an object with class
 \code{"summary.Renouv"}.
 
 The function \code{print.summary.Rendata} does not returns anything.
 
}
%%\references{}
%%\author{}
%%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  
  The model fitting function \code{\link{Renouv}} (to build
  \code{"Renouv"} model objects), \code{\link{summary}}.
  
}
\examples{
## use Brest data
fit <- Renouv(Brest)
summary(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
