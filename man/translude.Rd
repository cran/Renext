\name{translude}
\alias{translude}
\title{Make translucient colors}
\usage{
translude(colors, alpha = 0.6)
}
\arguments{
  \item{colors}{

    A vector of colors in a format that can be understood by
    \code{\link{col2rgb}}.

  }
  \item{alpha}{

    Vector of level(s) of opacity between \code{0} and \code{1} (\code{0} means fully
    transparent and \code{1} means opaque). After recycling to reach the
    required length, this value or vector is used as \code{alpha} in
    \code{\link{rgb}}.  } } \value{ A vector of translucient (or
    semi-transparent) colors.

}
\note{

  Using the \bold{RColorBrewer} package might be a better option that
  transluding chosen colors!

}
\description{
Make translucient colors.
}