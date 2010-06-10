\name{fweibull}
\alias{fweibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ML estimation of classical Weibull distribution}
\description{
  Fast Maximum Likelihood estimation of the classical two parameters
  Weibull distribution
}
\usage{
  fweibull(x,
           info.observed = FALSE,
           check.loglik = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Sample vector to be fitted}
  \item{info.observed}{Should the observed information matrix be used or
    the expected one be used? }
  \item{check.loglik}{If \code{TRUE}, the log-likelihood is recomputed
    using \code{dweibull} function with \code{log = TRUE}. The result is
    returned as a list element.}
  %% \item{\dots}{arguments XXX }
}
\details{
  The ML estimates are obtained thanks to a reparametrisation with
  \eqn{\eta = scale^{1/shape}}{eta = scale^{1/shape}}
  in place of \code{shape}. This allows the maximisation of a one-dimensional
  likelihood \eqn{L}{L} since the \eqn{\eta}{eta} parameter can be concentrated out
  of \eqn{L}{L}. This also
  allows the determination of the \emph{expected} information matrix for
  \eqn{[shape,\,\eta]}{[\code{shape}, \code{eta}]} rather than the usual \emph{observed} information. 
}
\value{
  A list
  \item{estimate}{parameter ML estimates}
  \item{sd}{(asymptotic) standard deviation for estimate}
  \item{cov}{(asymptotic) covariance matrix computed from theoretical or
    observed information matrix}
  \item{eta}{the estimated value for eta}
}
%%\references{ }
\author{Yves Deville}
\note{
 A well-accepted verdict is that expected information is better than
 observed information for problems with no mising data as it is the case
 here. We suspect  in the present estimation context that the expected information matrix
 has often a better condition number than its observed version.
 
}
\seealso{
  \code{\link{weibplot}} for Weibull plots
}
\examples{
library(MASS) 
n <- 1000
shape <- 2*runif(1)
x <- 100*rweibull(n, shape = 0.8, scale = 1)
res <- fweibull(x)
res2 <- fitdistr(x , "weibull")
est <- cbind(res$estimate, res2$estimate)
colnames(est) <- c("Renext", "MASS")
loglik <- c(res$loglik, res2$loglik)
est <- rbind(est, loglik)
est
weibplot(x,
         shape = c(res$estimate["shape"], res2$estimate["shape"]),
         scale = c(res$estimate["scale"], res2$estimate["scale"]),
         labels = c("Renext 'fweibull'", "MASS 'fitdistr'"),
         mono = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
