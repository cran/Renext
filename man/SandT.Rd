\name{SandT}
\alias{SandT}
\title{

  Compute empirical survivals (S) and return periods (T)

}
\usage{
  SandT(object, points = c("p", "H"), a = 0, naive = FALSE)
}
\arguments{
  
  \item{object}{

    The object containing the data, with class \code{"Renouv"} or
    \code{"Rendata"}.

  }
  \item{points}{

    Option for the computation of the plotting positions.  When
    \code{points} is set to \code{"p"}, the \eqn{p}-points formula is
    used with the selected value of \code{a}. This formula is used to
    compute the survival from which the return period is computed. When
    instead \code{points} is set to \code{"H"}, \emph{Nelson's formula}
    is used to compute the return periods, the survival value still
    being given by the \eqn{p}-points formula. When the data is
    heterogeneous, i.e.  when \code{object} contains \code{MAX} and/or
    \code{OTS} data, Nelson's formula is used only to compute the return
    periods of the upper slice of levels.

  }
  \item{a}{

    Parameter used in the interpolation formula for the inverse return
    periods as in Hirsch and Stedinger (1987).
    
  }
  \item{naive}{

    Logical. When \code{TRUE}, naive plotting positions are used to
    display \code{MAX} or \code{OTS} data. These can be defined only
    when a main sample exists in \code{object} as a \code{x.OT} element.
    For each \code{MAX} or \code{OTS} block, the positions use the
    number of events predicted using the rate of events as estimated
    from the main sample. When the main sample has a small durations,
    such predictions are likely to be misleading.
  }
  
}
\value{

  A list with the following elements

  \item{x}{

    Numeric vector containing the ordered values from all the available
    sources in the object: main sample, historical periods either 'MAX' or
    'OTS'.
    
  }
  \item{group, groupNames}{

    Integer and character vectors giving the source of the values in
    \code{x}, in the same order of the values. For instance,
    \code{group[10]} gives the group form which \code{x[10]} was
    extracted, and the name of this group is \code{groupNames[group[10]]}.
    
  }
  \item{S, T}{
    
    Numeric vectors of the same length as \code{x} and containing the
    corresponding estimation of the survival value and of the return
    period.

  }
  \item{thresh, lambda.thresh, S.thresh, T.thresh}{

    Vector of thresholds and the corresponding estimation for the event
    rate, survival and return period. All the estimations are \emph{for
      the threshold values}. The value of \code{T.thresh[i]} for a threshold
    \code{thresh[i]} results from a simple computation: divide the sum of
    the durations for blocks with thresholds \code{>= thresh[i]} by the
    number of events for these blocks.
    
  }
}
\description{

  Compute the empirical survival values and the empirical return periods
  at the observations of an object. These are used as plotting positions
  in several plots.

}
\details{
  
  When the object contains historical information (\code{MAX} or
  \code{OTS}), the computation is an adaptation Hirsch and Stedinger
  (1987) for the Marked Process (MP) context. The original method is
  devoted to block maxima and interpolates the survival values at the
  thresholds which are computed first. For MP, the interpolation is
  done for the inverse return periods, and the survival values are
  deduced from those of the inverse return periods.

  Nelson's formula provides unbiased estimates for the values of the
  cumulative hazard \eqn{H(x)} at the order statistics, and thus
  can be used to estimate the log-return periods as required on the
  return level plot.
  
}
\section{seealso}{
  The \code{\link{ppoints}} and \code{\link{Hpoints}} functions.
}

\section{Warning}{

  When using \code{points = "H"} the estimated values of the survival
  returned in \code{S} and those for the return period \code{T} no
  longer verify \code{T=1/S/lambda}, where \code{lambda} is the
  estimated rate. In this case, the values in \code{T} should be used in
  the return level plot, while the values in \code{S} should be used in
  the PP-plot.
  
}
\examples{
## use an object with class "Rendata"
ST1 <- SandT(object = Garonne)
## basic return level plot
plot(ST1$T, ST1$x, col = ST1$group, log = "x")
## use an object with class "Renouv"
fit <- Renouv(x = Garonne, plot = FALSE)
ST2 <- SandT(object = fit)
plot(ST2$T, ST2$x, col = ST2$group, log = "x")
}
\author{
  Yves Deville
}
\references{
  The original method for block maxima is described in
\itemize{
  \item Hirsch R.M. and Stedinger J.R.(1887) Plotting Positions
  for Historical Floods and their precision. \emph{Water
    Ressources Research}, vol. 23, N. 4 pp. 715-727.
  
  \item Millard S. and Neerchal N. (2001). \emph{Environmental
    Statistics with S-Plus}. CRC Press
}
The adaptation for the Marked Process context is
described in the \emph{Renext Computing Details}
document.
}
