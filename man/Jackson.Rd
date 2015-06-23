\name{Jackson}
\alias{Jackson}
\title{Jackson's statistic}
\usage{
   Jackson(x, norm = FALSE)
}
\arguments{
  \item{x}{Numeric vector or matrix. In the second case,
    rows are considered as samples.}
  
  \item{norm}{Logical: if \code{TRUE}, the statistic is
    normalized by using the \emph{asymptotic} mean and
  standard deviation, respectively 2 and 1.}
}
\value{
  A numeric vector of length \code{1} when \code{x} is a
  vector, or with length \code{nrow(x)} when \code{x} is a
  matrix.
}
\description{
  Jackson's statistic for the exponentiality test.
}
\details{
  The value(s) of the statistic are the ratio of two weighted
  means of the order statistics.
}
\references{
  J. Beirlant and T. de Weit and Y. Goegebeur(2006) A
  Goodness-of-fit Statistic for Pareto-Type Behaviour,
  \emph{J. Comp. Appl. Math.}, 186(1), pp. 99-116
}
