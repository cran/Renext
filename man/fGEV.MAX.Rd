\name{fGEV.MAX}
\alias{fGEV.MAX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  
  Fit a GEV distribution from block maxima or r largest order statistics
  using an aggregated Renewal POT process
  
}
\description{

  Fit a GEV distribution from block maxima or \eqn{r} largest order
  statistics using an aggregated Renewal POT process.

}
\usage{
fGEV.MAX(MAX.data, MAX.effDuration,
         MAX = NULL,
         control = list(maxit = 300, fnscale = -1),
         scaleData = TRUE,
         cov = TRUE,
         info.observed = TRUE,
         trace = 0)
}
\arguments{
  \item{MAX.data}{
    
    A list of block maxima or \eqn{r} largest statistics as in
    \code{\link[Renext]{Renouv}}.
    
  }
  \item{MAX.effDuration}{
    
    A vector of durations as in \code{\link[Renext]{Renouv}}. \emph{The
      durations must be identical} in order to have a common GEV
      distribution for the maxima.

  }
  \item{MAX}{
    
    A compact representation of the needed data as a list. This is
    typically created by using the (non exported)
    \code{Renext:::maxkeMAXdata} function.

  }
  \item{control}{

    List to control the optimisation in \code{\link[stats]{optim}}.

  }
  \item{scaleData}{
    
    Logical. If \code{TRUE}, the data in \code{MAX.data} are scaled
    before being used in the likelihood. The scaling operation is
    carried on the excesses (observations minus threshold). 

  }
  \item{cov}{

    Logical. If \code{TRUE} the standard deviation and the covariance
    matrix of the estimates are computed and returned as \code{sd} and
    \code{cov} elements of the list. However if the estimated shape
    parameter is \eqn{< - 0.5} the two elements are filled with
    \code{NA} because the regularity conditions can not be thought
    of as valid.
    
  }
  \item{info.observed}{

    Logical. If \code{TRUE} the covariance is computed from the
    \emph{observed} information matrix. If \code{FALSE}, the
    \emph{expected} information matrix is used instead. This is only
    possible for block maxima data, i.e.  when all the blocks contain
    only one observation. The computation relies on the formula given by
    Prescott and Walden. Note that the default value differs from that
    of the functions \code{\link{fGPD}}, \code{\link{flomax}} and
    \code{\link{fmaxlo}}, for historical reasons.
    
  }
  \item{trace}{

    Integer level of verbosity during execution. With the value \code{0},
    nothing is printed.
    
  }
  
}
\details{
  
  The data are assumed to provide maxima or \eqn{r} largest statistics
  arising from an aggregated renewal POT model with unknown event rate
  \eqn{\lambda} and unknown two-parameter Generalised Pareto
  Distribution for the excesses. A threshold \eqn{u} is fixed below the
  given data and the three unknown parameters \code{lambda},
  \code{scale} and \code{shape} of the POT model are found by maximising
  the likelihood.  Then a vector of the three parameters for the GEV
  distribution is computed by transformation. The covariance matrix and
  standard deviations are computed as well using the jacobian matrix of
  the transformation.

  The maximisation if for the log-likelihood with the rate \code{lambda}
  concentrated out, so it is a two-parameter optimisation.

}
\value{
  A list
  \item{estimate}{

    Named vector of the estimated parameters for the GEV distribution of
    maxima.

  }
  \item{opt}{

    Result of the optimisation.

  }
  \item{loglik}{

    Identical to \code{opt$value}. This is the maximised log-likelihood
    for the renewal POT model.
    
    %% It differs from the usual log-likelihood for GEV block maxima by
    %% a constant (w.r.t. parameters) depending on the number of blocks
    %% and the number of order statistics used.

  }
  \item{sd}{

    Standard deviation of the estimates (approximation based on the ML
    theory).

  }
  \item{cov}{

    Covariance matrix of the estimates (approximation based on the ML
    theory).

  }
}
\references{
  
  The \emph{Renext Computing Details} document.
  
  Prescott P. and Walden A.T. (1980) Maximum Likelihood Estimation of the 
  Parameters of the Generalized Extreme-Value Distribution.  
  \emph{Biometrika} \bold{67}(3), 723-724.

}
\author{
  
  Yves Deville

}
\section{Caution}{

  Whatever be the data, the log-likelihood is infinite (hence maximal)
  for any vector of GEV parameters with shape \eqn{\leq -1}{<= -1} and
  postive scale. Hence the log-likelihood should be maximised with a
  constraint on the shape, while unconstrained optimisation is used
  here. In practice, for numerical reasons, the estimate usually remains
  inside the \code{shape > -1} region. \emph{An estimation leading to shape
  \eqn{< -1} must be considered as meaningless}. An estimation with shape
  \eqn{< -0.5} should be considered with care.

}
\note{

  This function could get more arguments in the future.

}
\seealso{

  \code{\link[Renext]{Renouv}}.

}
\examples{
##====================================================================
## block maxima: simulated data and comparison with  the 'fgev'
## function from the 'evd' package
##====================================================================
set.seed(1234)
u <- 10
nBlocks <- 30
nSim <- 100   ## number of samples 
Par <- array(NA, dim = c(nSim, 3, 2),
             dimnames = list(NULL, c("loc", "scale", "shape"), c("MAX", "evd")))
LL <- array(NA, dim = c(nSim, 2),
            dimnames = list(NULL, c("MAX", "evd")))

for (i in 1:nSim) {
  rd <- rRendata(threshold = u,
                 effDuration = 1,
                 lambda = 12,
                 MAX.effDuration = rep(1, nBlocks),
                 MAX.r = rep(1, nBlocks),
                 distname.y = "exp", par.y = c(rate = 1 / 100))

  MAX <- Renext:::makeMAXdata(rd)
  fit.MAX <- fGEV.MAX(MAX = MAX)
  fit.evd <- fgev(x = unlist(MAX$data))
  Par[i, , "MAX"] <- fit.MAX$estimate
  Par[i, , "evd"] <- fit.evd$estimate
  LL[i, "MAX"] <- fit.MAX$loglik
  LL[i, "evd"] <- logLik(fit.evd)
}

##====================================================================
## r largest: use 'ismev::rlarg.fit' on the venice data set.
## NB 'venice' is taken from the 'evd' package here.
##====================================================================
\dontrun{ 
require(ismev);
fit1 <- ismev::rlarg.fit(venice)

## transform data: each row is a block
MAX.data <- as.list(as.data.frame(t(venice)))
## remove the NA imposed by the rectangular matrix format
MAX.data <- lapply(MAX.data, function(x) x[!is.na(x)])
MAX.effDuration <- rep(1, length(MAX.data))

fit2 <- fGEV.MAX(MAX.data = MAX.data,
                 MAX.effDuration = MAX.effDuration)

## estimates
est <- cbind(ismev = fit1$mle, RenextLab = fit2$estimate)
print(est)
# covariance
covs <- array(dim = c(2, 3, 3),
              dimnames = list(c("ismev", "RenextLab"),
                colnames(fit2$cov), colnames(fit2$cov)))
                
covs["ismev", , ] <- fit1$cov
covs["RenextLab", , ] <- fit2$cov
print(covs)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
