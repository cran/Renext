\name{fRenouv}
\alias{fRenouv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit a renewal ("Renouvellement") model}
\description{
  Fit a renewal ("Renouvellement") model to data
}
\usage{
   fRenouv(x.OT,
           sumw.BOT = 1.0,
           z.H = NULL,
           block.H = NULL,
           w.BH = NULL,
           x.U = NULL,
           w.U = NULL,
           distname.y = "exponential",
           fixed.par.y = NULL,
           start.par.y = NULL,
           force.start.H = FALSE,
           numDeriv = TRUE,
           threshold,
           trans.y = NULL,
           conf.pct = c(95, 70),
           prob = NULL,
           prob.max = 0.9995,
           pred.period = NULL,
           suspend.warnings = TRUE,
           control = list(maxit = 300, fnscale = -1),
           control.H = list(maxit = 300, fnscale = -1),
           trace = 0,
           plot = TRUE,
           main = "",
           ylim = NULL,
           ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  %%--------------------------------------------------------------------------
  \item{x.OT}{Vector giving the observations for the "Over the Threshold"
    part of the data. These are the true observations, and not the
    exceedances \code{y = x - threshold}.}
  %%--------------------------------------------------------------------------
  \item{sumw.BOT}{Effective duration.}
  %%--------------------------------------------------------------------------
  \item{z.H}{Vector of MAX data or historical data. These should be
    defined in accordance with the \code{block.H} arg: the values in
    \code{z.H} must represent the maxima observed during each
    block. Each \code{H}-block must contain at least one value.}
  %%--------------------------------------------------------------------------
  \item{block.H}{Factor (or integer vector) of length \code{length(z.H)}
    giving the blocks in correspondence with elements of \code{z.H}.
    Each element of \code{block.H} will usually describe wide
    periods (\code{"middle age"}, \code{"XVIth century"}, etc.).}
  %%--------------------------------------------------------------------------
  \item{w.BH}{Vector of time lengths (durations) for the blocks in
    \code{block.H}. Must be in the same unit as \code{sumw.OT}.}
  %%--------------------------------------------------------------------------
  \item{x.U}{A vector of block Upper or Unobserved levels. Each value
    corresponds to a time period with known duration given in
    \code{w.U}. This argument provides the information that the level in
    \code{x.U} was never exceeded during the period with corresponding
    duration in \code{w.U}.}
  %%--------------------------------------------------------------------------
  \item{w.U}{A vector of block durations for Upper or Unobserved
    levels. Must be of the same length as \code{x.U}.}
  %%--------------------------------------------------------------------------
  \item{distname.y}{Distribution name to fit to the threshold
    exceedances \code{y}.}
  %%--------------------------------------------------------------------------
  \item{fixed.par.y}{Named list of known (or fixed) parameter values for the
    \code{y}-distribution. This is still EXPERIMENTAL.}
  %%--------------------------------------------------------------------------
  \item{start.par.y}{Named list of parameter initial values for the
    \code{y}-distribution. Only used when the distribution does not
    belong to the list of special distributions.}
  %%--------------------------------------------------------------------------
  \item{force.start.H}{Logical. When \code{TRUE} the values in
    \code{start.par.y} (which must then be correct) will be used also as
    starting values in the maximisation of the global likelihood : OT data
    and historical data. This is useful e.g. when the historical data
    fall outside of the support for the distribution fitted without
    historical data. See below the \bold{Details} section.}
  %%--------------------------------------------------------------------------
  %%\item{distname.N}{NOT IMPLEMENTED YET. Character string specifying
  %%  the  distribution for the counts. It could be \code{"poisson"} or
  %%  \code{"negative binomial"} in the future but at present time ONLY \code{"poisson"}
  %%  is accepted.}
  %%--------------------------------------------------------------------------
  \item{numDeriv}{Logical: should the hessian be computed using the
    \code{numDeriv} package (value \code{TRUE}) or should it be taken
    from the results of \code{optim}?
  }
  %%--------------------------------------------------------------------------
  \item{threshold}{The threshold above which exceedances will be
    fitted. Exceedances are given (up to missing values problems) by
    \code{x.OT - threshold}.}
  %%--------------------------------------------------------------------------
  \item{trans.y}{Transformation of the levels \emph{before
      thresholding} (if not \code{NULL}). This is only possible with the \code{"exponential"}
    value \code{distname.y}. The two allowed choices are \code{"square"}
    and \code{"log"} meanin that the fitted (exponentially distributed)
    values are \code{x.OT^2
      -threshold^2} and \code{log(x.OT)-log(threshold)} respectively.
    The corresponding distributons for  \code{x.OT} may be called
    "square-exponential" and "log-exponential".
  }  
  %%--------------------------------------------------------------------------
  \item{conf.pct}{Character or numeric vector specifying the
    percentages for the confidence (bilateral) limits on quantiles.}
  %%--------------------------------------------------------------------------
  \item{prob}{Vector of probabilities.}
  %%--------------------------------------------------------------------------
  \item{prob.max}{Max value of probability for return level and
    confidence limits evaluations. This argument 'shortens' the default
    \code{prob} vector: values \code{> prob.max} in the default
    \code{prob} vector are omitted. Ignored when a \code{prob} argument
    is given.}
  %%--------------------------------------------------------------------------
  \item{pred.period}{A vector of "pretty" periods at which return level
    and probability will be evaluated and returned in the \code{pred}
    data.frame.}
  %%\item{gof}{a logical requiring goodness-of-fit tests and
  % diagnostices to be computed.}
  %%--------------------------------------------------------------------------
  \item{suspend.warnings}{If \code{TRUE}, the warnings will be suspended
    during optimisation steps. This is useful when the parameters
    are subject to constraints as is usually the case.}
  %%--------------------------------------------------------------------------
  \item{control}{A named list used in \code{\link{optim}} for the
    no-history stage (if any). Note that \code{fnscale = -1} says that 
    maximisation is required (not minimisation) and must not be changed!}
  %%--------------------------------------------------------------------------
  \item{control.H}{A named list used in \code{\link{optim}} for the
    historical stage (if any). }
  %%--------------------------------------------------------------------------
  \item{trace}{Level of verbosity. Value \code{0} prints nothing.}
  %%--------------------------------------------------------------------------
  \item{plot}{Draw a return level plot?}
  %%--------------------------------------------------------------------------
  \item{main}{Main title for the return level plot (defaults to empty title).}
  %%--------------------------------------------------------------------------
  \item{ylim}{Limits for the y axis (defaults to values computed from the data).}
  %%--------------------------------------------------------------------------
  \item{\dots}{Further args to be passed to \code{RLplot}. Should be
    removed in future versions.}
}
\details{
  Some distributions are considered in a special manner:
  \code{"exponential"}, \code{"weibull"} and
  \code{"gpd"} (possible only through the use of the \code{evd}
  package). For these distributions, it is not necessary to give
  starting values nor parameter names which are unambiguous. They will
  be called "special" distributions.
  
  Other distributions are possible, such as "gamma". Since the
  probability functions are then used in a "black-box" fashion, these
  distributions should respect the following
  \bold{formal requirements}:
  
  \enumerate{
    \item The name for the \emph{density}, \emph{distribution} and
    \emph{quantile} functions must obey to the \emph{classical "prefixing
    convention"}. Prefixes must be respectively: \code{"d"}, \code{"p"}, \code{"q"}.
    This rules applies for distribution of the \code{stats} package
    and those of many other packages such \code{evd}.
    
    \item \emph{First (main) argument must be vectorisable} in all three
    functions, i.e. a vector of \code{x}, \code{q} or \code{p} must be
    accepted by the density, the distribution and the quantile
    functions.
    
    \item \emph{The density must have a} \code{log} \emph{formal} argument. When
    \code{log} is \code{TRUE}, the log-density is returned instead of
    the density.
  }

  For such a distribution it is necessary to give arguments names
  in \code{start.par.y}. The arguments list must have exactly the required
  number of parameters for the family (e.g. \code{2} for \code{gamma}).
  In future version, some parameters will be allowed to be fixed (known)
  and the parameter set will be the reunion of those in \code{start.par.y}
  and those in \code{fixed.par.y}. Anyway, \emph{at least one parameter
  must be unknown} for the \code{y} part of model.
  
  \emph{Mathematical requirements} exist for a correct use of ML. They
  are refered to as "regularity conditions" in ML theory. Note that
  support of the distribution must be the set of non-negative real numbers.
  
  The estimation procedure differs according to the existence of
  historical (MAX) data.
  \enumerate{
  \item When no historical data is given, the whole set
  of parameters contains orthogonal subsets: a "point process" part
  concerning the process of events, and an "observation" part
  concerning the threshold exceedances part. The parameters can in this
  case be estimated separately. The rate of the Poisson process is estimated
  by the empirical rate, i.e. the number of events divided by the total
  duration given in \code{sumw.BOT}. The  "Over the Threshold"
  parameters are estimated from the exceedances computed as \code{x.OT} 
  minus the threshold.
  
  \item When historical data is given, the two parameter vectors must
  be coped with together in maximising the global likelihood. In this case, we
  begin the estimation ignoring the historical data and then use the
  estimates as starting values for the maximisation of the global
  likelihood. In some circumstances, the estimates obtained in the
  first stage can not be used with historical data because some of these
  fall outside the support of the distribution fitted. This can happen
  e.g. with a \code{distname.y = "gpd"} when historical data exceed
  \code{threshold} - \code{scale}/\code{shape} for the values of \code{shape}
  and \code{scale} computed in the first stage.
  }
}
\value{
  A list containing the results.
  \item{est.N}{Estimate(s) for the count \code{"N"} part. This estimate
    does not use the historical data, even if is available.}
  \item{est.y}{Estimate(s) for the exceedance \code{"y"} part. This
    estimate does not use the historical data, even if available.}
  \item{cov.N, cov.y}{are (co-)variances for the estimates above.}
  \item{estimate}{Estimate(s) for the whole set of parameters based on
    OT data \bold{and on historical data} if available.}
  \item{ks.test}{Kolmogorov-Smirnov goodness-of-fit test.}
  \item{ret.lev}{A data.frame containing return levels and confidence
    limits. The  corresponding probabilities are either provided by user
    or taken as default values.}
  \item{pred}{A data.frame similar to \code{ret.lev}, but with
    "pretty" return periods. These are taken as the provided values
    \code{pred.period} if any or are chosen as "round" multiples of the
    time unit (taken from \code{sumw.BOT}). The periods are chosen in order
    to cover periods ranging from 1/10 to 10 time units.}
  Other results are available. Use \code{names(result)} to see their
  list.
  
  Except in the the special case where \code{distname.y} is
  \code{"exponential"} and where no historical data are used, the
  inference on quantiles is obtained with the \emph{delta method} and
  using numerical derivatives. Confidence limits are unreliable for
  return levels much greater than the observation-historical period.
  
  Due to the presence of estimated parameters, the Kolmogorov-Smirnov
  test is unreliable when less than 30 observations are available.
  
}
\references{
Miquel J.(1984) \emph{Guide pratique d'estimation des
probabilites de crues}, Eyrolles (coll. EDF DER).
}
\author{Yves Deville}
\note{
  The model only concerns the "Over the Threshold" part of the
  distribution of the observations. When historical data is used,
  observations should all be larger than the threshold.

  The name of the list elements is indicative and is very likely to be
  changed in future versions. At the time, the effect of historical data
  on estimation (when such data exist) can be evaluated by comparing
  \code{c(res$est.N, res$est.y)} and \code{res$estimate} where
  \code{res} is the results list.

  Some warnings may indicate that missing values are met during the
  optimisation process. This is due to the evaluation of density at tail
  values. At the time the ML estimates are computed using an unconstrained
  optimisation, so invalid parameter values can be met during the
  maximisation or even be returned as (invalid) estimates. 

  Validity tests for the estimation in presence of historical have been
  limited at the time. Therefore this possibility should be regarded as
  experimental.
  
}
\section{Warning}{
  With some distributions or in presence of historical data, the
  estimation can fail due to some problem during the
  optimization or the determination of the (numerical) hessian can be impossible.
  This can happen if \emph{one or more parameter is too small}
  to compute a finite difference approximation of gradient. For instance
  the 'rate' parameter of the exponential distribution (= inverse mean)
  will be small when the mean of exceedances is large.

  A possible solution is then to \bold{rescale the data} e.g. dividing
  them  by 10 or 100. As a rule of thumb, an
  acceptable scaling leads to data (exceedances) of a few units to a few
  hundreds, but \bold{an order of magnitude of thousands or more should be
  avoided and reduced by scaling}. The rescaling is recommanded for the square
exponential distribution (obtained with \code{trans = "square"}) since
the observations are squared. 

  Another possible way to solve the problem is to change the
  \code{numDeriv} value. 

  This problem will be solved in a future version.
}
\seealso{\code{\link{rRenouv}} to simulate "renouvellement" data,
  \code{\link{RLplot}} for the \emph{return level plot}. See
  \code{\link{optim} for the tuning of optimisation.}
}
\examples{
## artificial data
test2 <- rRenouv(nb = 100,
                 threshold = 40,
                 par.N = list(lambda = 2),
                 densfun.y = "weibull",
                 par.y = mom2par("weibull", mean = 30, sd = 36))  

## exponential exceedances/no historical data/wrong threshold
fit.expon <- fRenouv(x.OT = test2$x, 
                     distname.y = "exponential",
                     threshold = 50,
                     conf.pct = c(70, 95, 99),
                     main = "\"exponential\" distr. (special)")

## The same but with regarding exponential distribution
## as an arbitrary (not in the listo of special distr.)
fit.exp <- fRenouv(x.OT = test2$x, 
                  distname.y = "exp",
                  start.par.y = list(rate = 1),
                  threshold = 50,
                  main = "\"exp\" distr. (non-special)")

## The same but with another non-target distribution 
## (thus requiering inital values in start.par.y) 
fit.gamma <- fRenouv(x.OT = test2$x,
                     distname.y = "gamma",
                     threshold = 50,
                     main = "\"gamma\" distr. (non-special)")


## Weibull exceedances/no historical data/true threshold
fit.weibull <- fRenouv(x.OT = test2$x,
                       distname.y = "weibull",
                       threshold = 50,
                       main = "\"weibull\" distr. (special). NO history")

## Add historical data 400 "years"
fit.weibull2  <-
   fRenouv(x.OT = test2$x,
           z.H = 300, block.H = 1, w.BH = 400,  ## HISTORY
           distname.y = "weibull",
           threshold = 40,
           main = "\"weibull\" distr. (special). History: 1 block of 400 years")
## comparison...
comp <- cbind(c(fit.weibull$est.N, fit.weibull$est.y), fit.weibull2$estimate)
colnames(comp) <- c("No hist.", "With hist.")
print(comp)

## gpd exeedances, like POT, but with history
\dontrun{
require(evd)
fit.gpd <-
    fRenouv(x.OT = test2$x, 
            z.H = 300, block.H = 1, w.BH = 400,     ## HISTORY
            distname.y = "gpd",
            threshold = 40,
            main = "\"gpd\" distr. (special). History: 1 block of 400 years")
}

\dontrun{

## Coles'book  rain example. Blocks are guessed and might not be exact
## due to missing values

library(ismev)
data(rain)
  
block <- 1:length(rain)
block <- 1914 + floor(block/365)
  
fit.rain <- fRenouv(x.OT = rain, 
                  distname.y = "gpd",
                  threshold = 30,
                  main = "\"gpd\" distr. (target) on ismev 'rain' dataset")
  
fit.rain.evd <- fpot(x = rain, threshold = 30, model ="gpd")
    
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
