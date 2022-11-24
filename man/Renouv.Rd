\name{Renouv}
\alias{Renouv}
\encoding{UTF-8}
\title{
  Fit a 'Renouvellement' model
}
\description{
  
  Fit a 'renouvellement' POT model using Over the Threshold data and
  possibly historical data of two kinds.

}
\usage{
Renouv(x,
       threshold = NULL,
       effDuration = NULL,
       distname.y = "exponential",
       MAX.data = NULL,
       MAX.effDuration = NULL,
       OTS.data = NULL,
       OTS.effDuration = NULL,
       OTS.threshold = NULL,
       fixed.par.y = NULL,
       start.par.y = NULL,
       force.start.H = FALSE,
       numDeriv = TRUE,
       trans.y = NULL,
       jitter.KS = TRUE,
       pct.conf = c(95, 70),
       rl.prob = NULL,
       prob.max = 1.0-1e-04 ,
       pred.period = NULL,
       suspend.warnings = TRUE,
       control = list(maxit = 300, fnscale = -1),
       control.H = list(maxit = 300, fnscale = -1),
       trace = 0,
       plot = TRUE,
       label = "",
       ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{x}{
    
    Can be a numeric vector, an object of the class \code{"Rendata"} or
    \code{NULL}.  In the first case, \code{x} contains all the levels
    above the threshold for a variable of interest.  In the second case,
    most formal arguments take values in accordance with the object
    content, and can be by-passed by giving the formal explicitly.  When
    \code{x} is \code{NULL}, the model is fitted using the data provided
    using the \code{OTS} and \code{MAX} formals.

  }
  
  \item{threshold}{

    Value of the threshold for the OT data.

  }

  \item{effDuration}{

    Effective duration, i.e. duration of the OT period.

  }
  \item{distname.y}{

    Name of the distribution for the excesses over the threshold. See
    \bold{Details} below.

  }
  \item{MAX.data}{

    Either a numeric vector or a list of numeric vectors representing
    historical data \eqn{r}{r}-max by blocks. When \emph{a vector} is
    given, there is only one block, and the data are the corresponding
    \eqn{r}{r}-max observed levels where \eqn{r}{r} is the vector
    length; the block duration is given in \code{MAX.effDuration}. When
    \emph{a list} is given, each list element contains the data for one
    block, and the effective duration are in \code{MAX.effDuration}

  }
  \item{MAX.effDuration}{

    Vector of (effective) durations, one by block MAX data.

  }
  \item{OTS.data}{
    
    A numeric vector or a list of numeric vectors representing
    supplementary Over Threshold data in blocks. When a \emph{vector} is
    given, there is only one block, and the data contain all the
    'historical' levels over the corresponding threshold given in
    \code{OTS.threshold}. The block duration is given in
    \code{OTS.effDuration}. When a \emph{list} is given, each list
    element contains the data for one block, and the threshold and
    effective duration are in \code{OTS.threshold} and
    \code{OTS.effDuration}.

  }
  \item{OTS.effDuration}{
    
    A numeric vector giving the (effective) durations for the OTS
    blocks.

  }
  \item{OTS.threshold}{

    A vector giving the thresholds for the different OTS blocks. The
    given values must be greater than or equal to the value of
    \code{threshold}.

  }
  \item{fixed.par.y}{

    Named list of known (or fixed) parameter values for the
    \code{y}-distribution.

  }
  \item{start.par.y}{

    Named list of parameter initial values for the
    \code{y}-distribution. Only used when the distribution does not
    belong to the list of special distributions.

  }
  \item{force.start.H}{
    
    Logical. When \code{TRUE}, the values in \code{start.par.y} (which
    must then be correct) will be used also as starting values in the
    maximisation of the global likelihood : OT data and historical
    data. This is useful e.g. when the historical data fall outside of
    the support for the distribution fitted without historical data. See
    below the \bold{Details} section.

  }
  \item{numDeriv}{
    
    Logical: should the hessian be computed using the \code{numDeriv}
    package (value \code{TRUE}) or should it be taken from the results
    of \code{optim}?

  }
  \item{trans.y}{
    
    Transformation of the levels \emph{before thresholding} (if not
    \code{NULL}).  This is only possible with the \code{"exponential"}
    value \code{distname.y}. The two allowed choices are \code{"square"}
    and \code{"log"} meaning that the fitted (exponentially distributed)
    values are \code{x.OT^2} \code{-threshold^2} and \code{log(x.OT)}
    \code{-log(threshold)} respectively.  The corresponding
    distributions for \code{x.OT} may be called "square-exponential" and
    "log-exponential".

  }  
  \item{jitter.KS}{

    Logical. When set to \code{TRUE}, a small amount of noise is added
    to the "OT" data used in the Kolmogorov-Smirnov test in order to
    remove ties. This is done using the \code{\link{OTjitter}} function.

  }
  \item{pct.conf}{
    
    Character or numeric vector specifying the percentages for the
    confidence (bilateral) limits on quantiles.

  }
  \item{rl.prob}{
    
    Vector of probabilities for the computation of return levels.  These
    are used in plots (hence must be dense enough) and appear on output
    in the data.frame \code{ret.lev}.
    
  }
  \item{prob.max}{

    Max value of probability for return level and confidence limits
    evaluations. This argument 'shortens' the default \code{prob}
    vector: values \code{> prob.max} in the default \code{prob} vector
    are omitted. Ignored when a \code{prob} argument is given.

  }
  \item{pred.period}{
    
    A vector of "pretty" periods at which return level and probability
    will be evaluated and returned in the \code{pred} data.frame.

  }
  \item{suspend.warnings}{
    
    Logical. If \code{TRUE}, the warnings will be suspended during
    optimisation steps. This is useful when the parameters are subject
    to constraints as is usually the case.

  }
  \item{control}{
    
    A named list used in \code{\link{optim}} for the no-history stage
    (if any). Note that \code{fnscale = -1} says that maximisation is
    required (not minimisation) and must not be changed!

  }
  \item{control.H}{
    
    A named list used in \code{\link{optim}} for the historical stage
    (if any).

  }
  \item{trace}{

    Level of verbosity. Value \code{0} prints nothing.

  }
  \item{plot}{

    Draw a return level plot?

  }
  \item{label}{

    Label to be used in the legend when \code{plot} is \code{TRUE}.

  }
  %%\item{main}{Main title for the return level plot (defaults to empty title).}
  %%\item{ylim}{Limits for the y axis (defaults to values computed from the data).}
  \item{\dots}{

    Arguments passed to \code{\link{plot.Renouv}}, e.g. \code{main},
    \code{ylim}.

  }
}
\details{
  
  The model is fitted using Maximum Likelihood (ML).
  
  Some distributions listed below and here called "special" are
  considered in a special manner.  For these distributions, it is not
  necessary to give starting values nor parameter names which are
  unambiguous.

  \tabular{ll}{
    distribution \tab parameters \cr %
    \code{exponential} \tab \code{rate} \cr%
    \code{weibull} \tab \code{shape}, \code{scale}\cr%
    \code{GPD} \tab \code{scale}, \code{shape}\cr%
    \code{gpd} \tab \code{scale}, \code{shape}\cr%
    \code{lomax} \tab \code{scale}, \code{shape}\cr%
    \code{maxlo} \tab \code{scale}, \code{shape}\cr%
    \code{log-normal} \tab \code{meanlog},
    \code{sdlog} \cr%
    \code{gamma} \tab \code{shape}, \code{scale}\cr%
    \code{mixexp2} \tab \code{prob1}, \code{rate1}, \code{delta}
  }
  
  Other distributions can be used. Because the probability functions are
  then used in a "black-box" fashion, these distributions should respect
  the following \emph{formal requirements}:
  
  \enumerate{

    \item The name for the \emph{density}, \emph{distribution} and
    \emph{quantile} functions must obey to the \emph{classical
    "prefixing convention"}. Prefixes must be respectively: \code{"d"},
    \code{"p"}, \code{"q"}.  This rules applies for distribution of the
    \code{stats} package and those of many other packages such
    \code{evd}.
    
    \item \emph{The first (main) argument must be vectorisable} in all
    three functions, i.e. a vector of \code{x}, \code{q} or \code{p}
    must be accepted by the density, the distribution and the quantile
    functions.
    
    \item \emph{The density must have a} \code{log} \emph{formal}
    argument. When \code{log} is \code{TRUE}, the log-density is
    returned instead of the density.

  }
  
  For such a distribution, it is necessary to give arguments names in
  \code{start.par.y}. The arguments list must have exactly the required
  number of parameters for the family (e.g. \code{2} for \code{gamma}).
  Some parameters can be fixed (known); then the parameter set will be
  the reunion of those appearing in \code{start.par.y} and those in
  \code{fixed.par.y}. Anyway, in the present version, \emph{at least one
  parameter must be unknown} for the \code{y} part of the model.
  
  \emph{Mathematical requirements} exist for a correct use of ML. They
  are referred to as "regularity conditions" in ML theory. Note that the
  support of the distribution must be the set of non-negative real
  numbers.
  
  The estimation procedure differs according to the existence of
  the different types of data: main sample, MAX and OTS.
  
  \enumerate{
    
    \item When no historical data is given, the whole set of parameters
    contains orthogonal subsets: a "point process" part
    concerning the process of events, and an "observation" part
    concerning the excesses over the threshold. The parameters can in this
    case be estimated separately. The rate of the Poisson process is estimated
    by the empirical rate, i.e. the number of events divided by the total
    duration given in \code{effDuration}. The  "Over the Threshold"
    parameters are estimated from the excesses computed as \code{x.OT} 
    minus the threshold.
    
    \item When historical data is given, the two parameter vectors must be
    coped with together in maximising the global likelihood. In this case,
    we begin the estimation ignoring the historical data and then use the
    estimates as starting values for the maximisation of the global
    likelihood. In some circumstances, the estimates obtained in the first
    stage can not be used with historical data because some of these fall
    outside the support of the distribution fitted. This can happen
    e.g. with a \code{distname.y = "gpd"} when historical data exceed
    \code{threshold} - \code{scale}/\code{shape} for the values of
    \code{shape} and \code{scale} computed in the first stage.
    
    \item From version 2.1-1 on, it is possible to use \code{OTS} and/or
    \code{MAX} data with no \code{OT} data by specifying \code{x =
      NULL}. Yet at the time this is only possible \code{distname.y} takes
    one of the two values: \code{"exp"}, or \code{"gpd"}.  The initial
    values for the parameter are then obtained by using the
    \code{\link{parIni.OTS}}, \code{\link{parIni.MAX}} functions and
    possibly by combining the two resulting initial parameter
    vectors. This possibility can be used to fit a model from \emph{block
      maxima} or \eqn{r}-\emph{largest} classical data but with more
    flexibility since the duration of the blocks may here not be constant.
    
  }
  
  The returned \code{Renouv} object contains a \code{MAX} element
  concerning the distribution of block maxima in the two following
  cases.
  
  \enumerate{ 

    \item When \code{distname.y} is \code{"exponential"} or \code{"exp"},
    the distribution of the maximum is Gumbel. The estimated parameters
    can be used with the \code{gumbel} function of the \bold{evd} package.
    
    \item When  \code{distname.y} is \code{"gpd"}, \code{"lomax"},
    \code{"maxlo"} or \code{"GPD"}  the distribution of the maximum
    is a Generalised Extreme Values distribution. The estimated parameters
    can be used with the \code{gev} function of the \bold{evd} package.
    
  }





 
}
\value{
  
  An object with class \code{"Renouv"}. This is mainly a list with the
  various results.

  \item{est.N}{

    Estimate(s) for the count \code{"N"} part. This estimate
    does not use the historical data, even if is available.

  }
  \item{est.y}{

    Estimate(s) for the excess \code{"y"} part. This estimate does
    not use the historical data, even if available.

  }
  \item{cov.N, cov.y}{

    The (co-)variances for the estimates above.

  }
  \item{estimate}{

    Estimate(s) for the whole set of parameters based on OT data
    \bold{and on historical data} if available.

  }
  \item{ks.test}{

    Kolmogorov-Smirnov goodness-of-fit test.

  }
  \item{ret.lev}{

    A data frame containing return levels and confidence limits. The
    corresponding probabilities are either provided by user or taken as
    default values.

  }
  \item{pred}{

    A data frame similar to \code{ret.lev}, but with "pretty" return
    periods. These are taken as the provided values \code{pred.period}
    if any or are chosen as "round" multiples of the time unit (taken
    from \code{effDuration}). The periods are chosen in order to cover
    periods ranging from 1/10 to 10 time units.

  }
  \item{MAX}{

    A list providing the estimated distribution of the maximum of the
    marks over a block of unit duration. This list element only exists when
    this distribution can be deduced from the fit, which
    is the case when \code{distname.y} is a GPD in a broad sense,
    see \bold{Details}.
      
  }

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
  \itemize{
    
    \item Miquel J. (1984) \emph{Guide pratique d'estimation des
      \enc{probabilités}{probabilites} de crues}, Eyrolles (coll. EDF DER).
    
    \item Coles S. (2001) \emph{Introduction to Statistical Modelling of Extremes
      Values}, Springer.

    \item Embrechts P., Klüppelberg C. and Mikosch T. (1997) \emph{Modelling
      Extremal Events for Insurance and Finance}. Springer.
  }
}
\author{

  Yves Deville

}
\note{

  The model only concerns the "Over the Threshold" part of the
  distribution of the observations. When historical data is used,
  observations should all be larger than the threshold.

  The name of the elements in the returned list is indicative, and is
  likely to be changed in future versions. At the time, the effect of
  historical data on estimation (when such data exist) can be evaluated
  by comparing \code{c(res$est.N, res$est.y)} and \code{res$estimate}
  where \code{res} is the results list.

  Some warnings may indicate that missing values are met during the
  optimisation process. This is due to the evaluation of the density at tail
  values. At the time the ML estimates are computed using an unconstrained
  optimisation, so invalid parameter values can be met during the
  maximisation or even be returned as (invalid) estimates. 

}
\section{Warning}{
  
  With some distributions or in presence of historical data, the
  estimation can fail due to some problem during the optimisation. Even
  when the optimisation converges, the determination of the (numerical)
  hessian can be impossible: This can happen if \emph{one or more
  parameter is too small} to compute a finite difference approximation
  of gradient. For instance the 'rate' parameter of the exponential
  distribution (= inverse mean) will be small when the mean of
  the excesses is large.

  A possible solution is then to \bold{rescale the data} e.g. dividing
  them by 10 or 100. As a rule of thumb, an acceptable scaling leads to
  data (excesses) of a few units to a few hundreds, but \bold{an
  order of magnitude of thousands or more should be avoided and reduced
  by scaling}. The rescaling is recommended for the square exponential
  distribution (obtained with \code{trans =} \code{"square"}) since the
  observations are squared.

  Another possible way to solve the problem is to change the
  \code{numDeriv} value.

}
\seealso{
  
  \code{\link{RLplot}} for the \emph{return level plot}. See
  \code{\link[stats]{optim}} for the tuning of the optimisation. The
  \code{\link{RenouvNoEst}} can be used to create an object with S3
  class \code{"Renouv"} from known parameters.

}
\examples{
## Garonne data. Use a "Rendata" object as 'x'. Historical data are used!
fit <- Renouv(x = Garonne, distname = "weibull", trace = 1,
              main = "'Garonne' data")
summary(fit)

## generates a warning because of the ties
fit2 <- Renouv(x = Garonne, distname = "GPD",
               jitter.KS = FALSE,
               threshold = 2800, trace = 1,
               main = "'Garonne' data with threshold = 2800 and GPD fit")

## use a numeric vector as 'x'
fit3 <-
    Renouv(x = Garonne$OTdata$Flow,
           threshold = 2500,
           effDuration = 100,
           distname = "GPD",
           OTS.data = list(numeric(0), c(6800, 7200)),
           OTS.effDuration = c(100, 150),
           OTS.threshold = c(7000, 6000), 
           trace = 1,
           main = "'Garonne' data with artificial \"OTS\" data")
## Add historical (fictive) data
fit4 <- Renouv(x = Garonne$OTdata$Flow,
               threshold = 2500,
               effDuration = 100,
               distname = "weibull",
               fixed.par.y = list(shape = 1.1),
               OTS.data = list(numeric(0), c(6800, 7200)),
               OTS.effDuration = c(100, 150),
               OTS.threshold = c(7000, 6000),
               trace = 0,
               main = "'Garonne' data with artificial \"OTS\" data")

##============================================================================
## use the 'venice' dataset in a r-largest fit from the 'evd' package
##============================================================================
## transform data: each row is a block
MAX.data <- as.list(as.data.frame(t(venice)))
## remove the NA imposed by the rectangular matrix format
MAX.data <- lapply(MAX.data, function(x) x[!is.na(x)])
MAX.effDuration <- rep(1, length(MAX.data))

## fit a Renouv model with no OT data. The threshold
## must be in the support of the gev distribution
u <- 66
fit.gpd <- Renouv(x = NULL,
                  MAX.data = MAX.data,
                  MAX.effDuration = MAX.effDuration,
                  distname.y = "GPD",
                  threshold = u,
                  numDeriv = FALSE,
                  trace = 0,
                  plot = FALSE)
\dontrun{
  require(ismev)
  ## compare with results from the ismev package 
  fit.gev <- rlarg.fit(venice)
  est.gev <- fit.gev$mle
  names(est.gev) <- c("loc", "scale", "shape")
  
  ## transform the 'gev' fit into a Ren parameter set.
  cov.gev <- fit.gev$cov
  rownames(cov.gev) <- colnames(cov.gev) <-  c("loc", "scale", "shape")
  trans <- gev2Ren(est.gev,
                   threshold = u,
                   vcovGev = cov.gev)
  est <- cbind(ismev = trans, RenextLab = coef(fit.gpd))
  colnames(est) <- c("ismev", "RenextLab")
  print(est)
  
  ## fill a 3d array with the two gpd covariance matrices
  cov2 <- attr(trans, "vcov")[c(1, 3, 4), c(1, 3, 4)]
  
  ## covariance
  covs <-
    array(dim = c(2, 3, 3),
          dimnames = list(c("ismev", "RenextLab"),
            colnames(fit.gpd$cov), colnames(fit.gpd$cov)))
  
  covs["ismev", , ] <- cov2
  covs["RenextLab", , ] <- fit.gpd$cov
  print(covs)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
