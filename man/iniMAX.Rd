\name{parIni.MAX}
\alias{parIni.MAX}
\alias{parIni.OTS}
\title{
  
  Initial estimation of GPD parameters for an aggregated renewal model

}
\description{

  Initial estimation for an aggregated renewal model with GPD marks.

}
\usage{
parIni.MAX(MAX, threshold, distname.y = "exp")
parIni.OTS(OTS, threshold, distname.y = "exp")
}
\arguments{
  \item{MAX}{

    A list describing partial observations of \code{MAX} type. These are
    block maxima or block \eqn{r}-largest statistics. The list must
    contain elements named \code{block}, \code{effDuration}, \code{r},
    and \code{data} (a by-block list of \eqn{r}-largest statistics).

  }
  \item{OTS}{
    
    A list describing partial observations of \code{OTS} type. These are
    observations above the block thresholds, each being not smaller than
    \code{threshold}. This list contains \code{block},
    \code{effDuration}, \code{threshold}, \code{r} and \code{data} (a
    by-block list of observations).

  }
  \item{threshold}{
    
    The threshold of the POT-renewal model. This is the location
    parameter of the marks from which Largest Order Statistics are
    observed.

  }
  \item{distname.y}{

    The name of the distribution. For now this can be \code{"exp"} or
    \code{"exponential"} for exponential excesses implying Gumbel
    block maxima, or \code{"gpd"} for GPD excesses implying GEV block
    maxima. The initialisation is the same in all cases, but the result
    is formatted according to the target distribution.

  }
  
}
\details{

  The functions estimate the Poisson rate \code{lambda} along with the
  shape parameter say \code{sigma} for exponential excesses.  If the
  target distribution is GPD, then the initial shape parameter is taken
  to be 0.
  
  In the "MAX" case, the estimates are obtained by regressing the maxima
  or \eqn{r}-Largest Order Statistics within the blocks on the
  log-duration of the blocks. The response for a block is the minimum of
  the \eqn{r} available Largest Order Statistics as found within the
  block, and \eqn{r} will generally vary across block.  When some blocks
  contain \eqn{r > 1} largest order statistics, the corresponding
  spacings are used to compute a spacings-based estimation of
  \eqn{\sigma}. This estimator is independent of the regression
  estimator for \eqn{\sigma} and the two estimators can be combined in a
  weighted mean.
  
  In the "OTS" case, the estimate of \code{lambda} is obtained by a
  Poisson regression using the log durations of the blocks as an
  offset. The estimate of \code{sigma} is simply the mean of all the
  available excesses, which by assumption share the same exponential
  distribution.
  
}
\value{
  
  A vector containing the estimate of the Poisson rate \eqn{\lambda},
  and the estimates of the parameters for the target distribution of the
  excesses. For exponential excesses, these are simply a
  \code{rate} parameter. For \code{GPD} excesses, these are the
  \code{scale} and \code{shape} parameters, the second taken as zero.

}
\references{

  See the document \emph{Renext Computing Details}.

}
\author{
  
  Yves Deville

}
\note{

  In the MAX case, the estimation is possible only when the number of
  blocks is greater than \eqn{1}, since otherwise no information about
  \eqn{\lambda} can be gained from the data; recall that the time at
  which the events occurred within a block is not known or used.

}
\seealso{
  
  The \code{\link{spacings}} methods for the spacings used in the
  estimation.

}
\examples{
set.seed(1234)
## initialisation for 'MAX' data
u <- 10
nBlocks <- 30
nSim <- 100
ParMAX <- matrix(NA, nrow = nSim, ncol = 2)
colnames(ParMAX) <- c("lambda", "sigma")

for (i in 1:nSim) {
  rd <- rRendata(threshold = u,
                 effDuration = 1,
                 lambda = 12,
                 MAX.effDuration = c(60, rexp(nBlocks)),
                 MAX.r = c(5, 2 + rpois(nBlocks, lambda = 1)),
                 distname.y = "exp", par.y = c(rate = 1 / 100))

  MAX <- Renext:::makeMAXdata(rd)
  pari <- parIni.MAX(MAX = MAX, threshold = u)
  ParMAX[i, ] <- pari   
}
## the same for OTS data
u <- 10
nBlocks <- 10
nSim <- 100
ParOTS <- matrix(NA, nrow = nSim, ncol = 2)
colnames(ParOTS) <- c("lambda", "sigma")
rds <- list()

for (i in 1:nSim) {
  rd <- rRendata(threshold = u,
                 effDuration = 1,
                 lambda = 12,
                 OTS.effDuration = rexp(nBlocks, rate = 1 / 10),
                 OTS.threshold = u + rexp(nBlocks, rate = 1 / 10),
                 distname.y = "exp", par.y = c(rate = 1 / 100))
  rds[[i]] <- rd
  OTS <-  Renext:::makeOTSdata(rd)
  pari <- parIni.OTS(OTS = OTS, threshold = u)
  ParOTS[i, ] <- pari   
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
