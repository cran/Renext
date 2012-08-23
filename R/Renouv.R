##===================================================================
## Retrieves colnames for confidence or prediction bands
## from a dataframe given in 'x'. The levels are given in
## the colnames
##
## NEW IN VERSION  1.3-5
##===================================================================

predNames <- function(x, prefix = c("L", "U")) {

  cn <- names(x)
  mn <- character(0)
  pct <- character(0)
  type <- character(0)
  
  ## find the narrowest confint
  if ("L" %in% prefix) {
    ##pct.L <- grep("L.[0-9]*$", cn, value = TRUE)
    pct.L <-  grep("^L.[0-9]*\\.?[0-9]*$", cn, value = TRUE)
    if (ln <- length(pct.L)) {
      mn <- c(mn, pct.L)
      type <- c(type, rep("L", ln))
      pct <- c(pct, gsub("L.", "", pct.L))
    }
  }
  if ("U" %in% prefix) {
    pct.U <- grep("^U.[[:digit:]]*\\.?[[:digit:]]*$", cn, value = TRUE)
    if (ln <- length(pct.U)) {
      mn <- c(mn, pct.U)
      type <- c(type, rep("U", ln))
      pct <- c(pct, gsub("U.", "", pct.U))
    }
  }
    
  list(names = mn,
       type = type,
       pct = as.numeric(pct))
  
}

##===================================================================
## Compute limits for levels using confidence or prediction bands
## as well as historical data.
##
## NEW IN VERSION 1.4-0
##
##===================================================================

rangeLev.Renouv <- function(x,
                            show.MAX = TRUE,
                            show.OTS = TRUE,
                            Tlim = NULL) {

  pn <- predNames(x$ret.lev, prefix = c("L", "U"))
  pnn <- c("quant", pn$names)

  ## cat("pn = ", pnn, "\n")
  Tlim <- range(x$pred$period, Tlim)
  ## cat("x$pred = \n"); print(x$pred)
  ## cat("Tlim = \n"); print(Tlim)
  
  if (!is.null(Tlim)) {
    period <- x$ret.lev$period
    ind <- (period >= Tlim[1]) & (period <= Tlim[2])
    r <- range(x$ret.lev[ind, pnn], x$x.OT)
  } else {
    r <- range(x$ret.lev[ , pnn], x$x.OT)
  }
    
  if (x$history.MAX$flag && show.MAX) {
    x1 <- unlist(x$history.MAX$data)
    if (length(x1)) r <- range(r, range(x1))
  } 
  
  if (x$history.OTS$flag && show.OTS) {
    x1 <- unlist(x$history.OTS$threshold)
    if (length(x1)) r <- range(r, range(x1))
    x1 <- unlist(x$history.OTS$data)
    if (length(x1)) r <- range(r, range(x1))
  }
  
  r 
   
}

##=====================================================================
## round quantiles and confidence limits using a 
## suitable number of digits
##=====================================================================

roundPred <- function(pred, dig.quant = NA) {

  cn <- colnames(pred)
  ## find the narrowest confint
  pct.L <- grep("L.[0-9]", cn, value = TRUE)
  pct.U <- grep("U.[0-9]", cn, value = TRUE)
  
  if (is.na(dig.quant) || length(dig.quant) == 0) {
    
    disp <- FALSE
    if (length(pct.L)) {
      i <- which.min(substring(pct.L, first = 3))
      vn1 <- pct.L[i]
      disp <- TRUE
    } else vn1 <- "quant"
    if (length(pct.U)) {
      i <- which.min(substring(pct.U, first = 3))
      vn2 <- pct.U[i]
      disp <- TRUE
    } else vn2 <- "quant"
    
    prec <-  mean(pred[ ,"quant"])/100
    if (disp)  prec <- pmin(prec, min(pred[ ,vn2] - pred[, vn1]))
    dig.quant <- -floor(log(prec, base = 10))
    if (dig.quant < 0) dig.quant <-  0 
  } 

  ## roudn all selected variables
  pred.mod <- pred
  ind <- c("quant", pct.L, pct.U)
  pred.mod[ , ind] <- round(pred.mod[ , ind], digits = dig.quant) 
  pred.mod
  
}

##===================================================================
## add a small amount of noise to x
##
##===================================================================

OTjitter <- function(x, threshold = NULL) {
  d <- diff(sort(x))
  signoise <- pmin(min(d[d>0]) / 5, mean(abs(x))/500)
  mynoise  <- rnorm(length(x), mean = 0, sd = signoise)
  x.noised <- x + mynoise
  if (length(threshold) > 0) {
    if (any(x < threshold)) stop("'threshold' must be a numeric <= min(x)")
    ind <- x.noised < threshold
    x.noised[ind] <- x[ind] + abs(mynoise[ind])
  }
  x.noised  
}

##======================================================================
## summary method for 'Renouv'
##======================================================================

summary.Renouv <- function(object,
                           correlation = FALSE,
                           symbolic.cor = FALSE,
                           ...) {

  ## distribution
  ans <- object

  ## improve info on degree of freedom
  ans$df <- as.integer(c(object$p, object$df))
  names(ans$df) <- c("par", "obs")
  ## coefficients: take care for fixed ones!
  
  est <- object$estimate
  se <- sqrt(diag(object$cov))
  tval <- est/se
    
  ## correction odf df due to historical data
  rdf <- object$df - object$p
  
  ans$coefficients <-
    cbind(est, se, tval, 2*pt(abs(tval), rdf, lower.tail = FALSE))
  
  dimnames(ans$coefficients) <-
    list(names(est), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  
  ans$pred <- roundPred(object$pred)

  if (correlation) {
    ans$correlation <- object$corr
    ans$symbolic.cor <- symbolic.cor  
  }
 
  class(ans) <- "summary.Renouv"
  ans

}

print.summary.Renouv <-
  function(x,
           coef = TRUE,
           pred = TRUE,
           probT = FALSE,
           digits = max(3, getOption("digits") - 3),
           symbolic.cor = x$symbolic.cor,
           signif.stars = getOption("show.signif.stars"),
           ...) {
  
  ## cat(sprintf("o Number of OT observations : %d\n", length(x$y.OT)))
  cat(sprintf(paste("o Main sample 'Over Threshold'\n",
                    "   . Threshold        %8.2f\n",
                    "   . Effect. duration %8.2f years\n",
                    "   . Nb. of exceed.   %5d\n\n",
                    collapse = ""),
              x$threshold, x$effDuration, length(x$y.OT)))
 
  cat(sprintf("o Estimated rate 'lambda' for Poisson process (events): %5.2f evt/year.\n\n",
              x$estimate["lambda"]))
 
  cat(sprintf("o Distribution for exceedances y: \"%s\", with %d par. ",
              x$distname.y, x$p.y))
  cat(paste(sprintf("\"%s\"", x$parnames.y), collapse = ", "), "\n\n")

  if (x$transFlag){
    cat(sprintf("o Transformation applied: \"%s\"\n\n",
                x$trans.y))
  } else {
    cat("o No transformation applied\n\n")
  }
    
  if (coef) {
    cat("o Coefficients\n\n")

    if (!probT) {
      print(x$coefficients[ , 1:3])

    } else {
      printCoefmat(x$coefficients,
                   digits = digits,
                   signif.stars = signif.stars,
                   na.print = "NA", ...)
    }
    cat("\n")
    cat(sprintf("Degrees of freedom: %d (param.) and %d (obs)\n", x$df["par"], x$df["obs"]))
    cat("\n")
  }
  
  if (any(x$fixed)) {
    cat("The following coef. were fixed\n")
    print(names(x$fixed)[x$fixed])
    cat("\n")
  }
  
  cat(sprintf("o Inference method used for return levels\n\"%s\"\n\n",
              x$infer.method))
  
  if (pred) {
    cat("o Return levels\n\n")
    print(x$pred)
    cat("\n\n")
  }
  
  if (x$history.MAX$flag) {
    cat(sprintf("o 'MAX' historical info: %d  blocks, %d obs., total duration = %5.2f years\n\n",
                nlevels(x$history.MAX$block),
                length(unlist(x$history.MAX$data)),
                sum(x$history.MAX$effDuration)))
    
    for ( i in 1L:nlevels(x$history.MAX$block) ) {

      Zi <- x$history.MAX$data[[i]]
      ri <- length(Zi)
      
      if ( ri <= 12L ) {
        obs.str <- paste(format(Zi),  collapse = ", ")
      } else {
        obs.str <- paste(paste(format(Zi[1L:3L]), collapse = ", "),
                         "...",
                         paste(format(Zi[(ri-2L):ri]), collapse = ", "),
                         sep = ", ")
      }
      cat(sprintf("  * block %d, %5.2f years, %d obs.\n\n     %s\n\n",
                  i, x$history.MAX$effDuration[[i]], ri, obs.str))
      
    }
    
  } else {
    cat("o no 'MAX' historical data\n\n")
  }
  
  if (x$history.OTS$flag) {
    cat(sprintf("o 'OTS' historical info: %d  blocks, %d obs., total duration = %5.2f years\n\n",
                nlevels(x$history.OTS$block),
                length(unlist(x$history.OTS$data)),
                sum(x$history.OTS$effDuration)))
    
    for ( i in 1L:nlevels(x$history.OTS$block) ) {
      
      Zi <- x$history.OTS$data[[i]]
      ri <- length(Zi)
      
      if ( (ri > 1L) && (ri <= 12L) ) {
        obs.str <- paste(format(Zi),  collapse = ", ")
      } else {
        obs.str <- paste(paste(format(Zi[1L:3L]), collapse = ", "),
                         "...",
                         paste(format(Zi[(ri-2L):ri]), collapse = ", "),
                         sep = ", ")
      }
      cat(sprintf("  * block %d, %5.2f years, thresh. %8.2f, %d obs.\n\n    %s\n\n",
                  i, x$history.OTS$effDuration[[i]],
                  x$history.OTS$threshold[[i]],
                  ri, obs.str))
      
    }

  } else {
    cat("o no 'OTS' historical data\n\n")
  }

  ## 2012-12-22 PREVIOUS CODE
  if (FALSE) {
    
    if (x$history.MAX$flag) {
      cat(sprintf("o 'MAX' historical info: %d blocks, %d obs., total duration = %5.2f years\n\n",
                  nlevels(x$history.MAX$block),
                  length(unlist(x$history.MAX$data)),
                  sum(x$history.MAX$effDuration)))
      
    } else  cat("o no 'MAX' historical data\n\n")
    
  
    if (x$history.OTS$flag) {
      cat(sprintf("o 'OTS' historical info: %d blocks, %d obs., total duration = %5.2f years\n\n",
                  nlevels(x$history.OTS$block),
                  length(unlist(x$history.OTS$data)),
                sum(x$history.OTS$effDuration)))
    } else  cat("o no 'OTS' historical data\n\n")

  } 
  
  
  cat("o Kolmogorov-Smirnov test\n")
  print(x$KS.test)
  cat("\n")

  ## copied from 'print.summary.lm' in the 'stats' package
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1L) {
      cat("o Correlation of Coefficients:\n")
      if(is.logical(symbolic.cor) && symbolic.cor) {# NULL < 1.7.0 objects
        print(symnum(correl, abbr.colnames = NULL))
        cat("\n")
      } else {
        correl <- format(round(correl, 2), nsmall = 2, digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
        cat("\n")
      }
    }
  }
  
  
}

##===========================================================================
##
##===========================================================================

format.summary.Renouv <-
  function(x,
           ...) {
  
  ## cat(sprintf("o Number of OT observations : %d\n", length(x$y.OT)))
  text <- sprintf(paste("o Main sample 'Over Threshold'\n",
                    "   . Threshold        %8.2f\n",
                    "   . Effect. duration %8.2f years\n",
                    "   . Nb. of exceed.   %5d\n\n",
                    collapse = ""),
              x$threshold, x$effDuration, length(x$y.OT))
 
  text <- paste(text,
                sprintf("o Estimated rate 'lambda' for Poisson process (events): %5.2f evt/year.\n\n",
                        x$estimate["lambda"]))
 
  text <- paste(text, sprintf("o Distribution for exceedances y: \"%s\", with %d par. ",
                              x$distname.y, x$p.y))

  text <- paste(text, paste(sprintf("\"%s\"", x$parnames.y), collapse = ", "), "\n\n")

  if (any(x$fixed)) {
    cat("The following coef. were fixed\n", names(x$fixed)[x$fixed], "\n")
  }

  if (x$history.MAX$flag) {
    text <- paste(text,
                  sprintf("o 'MAX' historical info: %d blocks, %d obs., total duration = %5.2f years\n\n",
                          nlevels(x$history.MAX$block), length(unlist(x$history.MAX$data)),
                          sum(x$history.MAX$effDuration)))
    
  } else  text <- paste(text, "o no 'MAX' historical data\n\n")
  
  
  if (x$history.OTS$flag) {
    text <- paste(text, sprintf("o 'OTS' historical info: %d blocks, %d obs., total duration = %5.2f years\n\n",
                                nlevels(x$history.OTS$block), length(unlist(x$history.OTS$data)),
                                sum(x$history.OTS$effDuration)))
  } else  text <- paste(text, "o no 'OTS' historical data\n\n")
  
  text <- paste(text, sprintf("o Kolmogorov-Smirnov test\n   D = %6.4f, p-value = %6.4f\n",
                              x$KS.test$statistic, x$KS.test$p.value))
  
  text
  
}

##====================================================================
## coefficients method
##
##====================================================================

coef.Renouv <- function(object, ...) {
  object$estimate
}

##====================================================================
## Author: Y. Deville
##
## The black-box maximization of the log-likelhood is inspired
## 'fitdistr' of Brian Ripley (package MASS)
##
##====================================================================

Renouv <- function(x,
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
                   prob.max = 1.0-1e-4,
                   pred.period = NULL,
                   suspend.warnings = TRUE,
                   control = list(maxit = 300, fnscale = -1),
                   control.H = list(maxit = 300, fnscale = -1),
                   trace = 0,
                   plot = TRUE,
                   ## main = "",
                   ## ylim = NULL,
                   ...) {    
  
  ## for numerical differentiation  
  eps <- sqrt(.Machine$double.eps) ## seems too small (for 2nd order diff)
  eps <- 1e-6

  mc <- match.call()
  
  if (is(x, "Rendata")) {
    if(trace) cat("processing the 'Rendata' object 'x'\n\n")
    vn <- x$info$varName
    x.OT <- x$OTdata[ , vn]
    if (is.null(effDuration)) effDuration <- x$OTinfo$effDuration
    if (is.null(threshold)) threshold <- x$OTinfo$threshold
  } else{
    x.OT <- x
    if (length(effDuration) == 0 || is.na(effDuration) || (effDuration < 0))
      stop("a valid 'effDuration' must be given")
    if (length(threshold) == 0 || is.na(threshold) )
      stop("a valid 'threshold' must be given")
  }
  
  ## check and make MAXdata info
  if (!missing(MAX.data)) {
    MAX <- makeMAXdata(x = x,
                       data = MAX.data,
                       effDuration = MAX.effDuration)
  } else {
    MAX <- makeMAXdata(x = x, effDuration = MAX.effDuration)
  }
  
  ## check and make OTSdata info
  if (!missing(OTS.data)) {
    OTS <- makeOTSdata(x = x,
                       data = OTS.data,
                       threshold = OTS.threshold,
                       effDuration = OTS.effDuration)
  } else {
    OTS <- makeOTSdata(x = x,
                       threshold = OTS.threshold,
                       effDuration = OTS.effDuration)
  }
  
  if ( OTS$flag && (any(OTS$threshold < threshold)) )
    stop("OTS thresholds must be >= threshold")
  
  ##=================================================================
  ## prepare transforms if wanted/possible
  ##=================================================================
  
  if (!is.null(trans.y)) {
    if( !is.character(trans.y) || !(trans.y %in% c("square", "log")) ) 
      stop("trans.y must be NULL or be character in c( \"square\', \"log\")")
    else if (distname.y != "exponential") {
      stop("non-null value for 'trans.y' is only allowed when distname.y == \"exponential\"") 
    } else {
      transFlag <- TRUE
      if (trans.y == "square") {
        transfun <- function(x) x^2 
        invtransfun <- get("sqrt", mode = "function")   
      }
      if (trans.y == "log") {
        transfun <- get("log", mode = "function")   
        invtransfun <- get("exp", mode = "function")   
      }
    }
  } else {
    transFlag <- FALSE
    transfun <- NULL
    invtransfun <- NULL
  }

  ##=================================================================
  ## default rl.prob for quantile and confidence lims 
  ## should densify near 0 and 1
  ##=================================================================
  
  if (is.null(rl.prob)) {

    rl.prob <- c(0.0001,
                 seq(from = 0.01, to = 0.09, by = 0.01),
                 seq(from = 0.10, to = 0.80, by = 0.10),
                 seq(from = 0.85, to = 0.99, by = 0.01),
                 0.995, 0.996, 0.997, 0.998, 0.999,
                 0.9995, 0.9996, 0.9997, 0.9998, 0.9999,
                 0.99999, 0.999999)
    
    rl.prob <- rl.prob[rl.prob <= prob.max]
    
  } else {
    if (any(is.na(rl.prob))) stop("'rl.prob' values can not be NA") 
    if ( any(rl.prob <= 0.0) || any(rl.prob >= 1.0) ) stop("'rl.prob' values must be >0 and <1") 
    rl.prob <- sort(rl.prob)
  }
  
  if (is.null(pred.period)) {
    rr <- round(log(effDuration)/log(10))
    pred.period <- (10^rr)*c(0.1, 0.2, 0.5, 1:10)
  } else {
    if (any(is.na(pred.period))) stop("'pred.period' values can not be NA") 
    pred.period <- sort(pred.period)
  }
 
  ##===============================================================
  ## CODING RULES
  ## .OT    an object related to OT data
  ## .BOT   an object related to  BLOCK OT, typically a vector with
  ##        length = number of blocks
  ##===============================================================
  
  ## First some essential elements
  nb.x <- length(x.OT)
  ind <- (x.OT > threshold) & !is.na(x.OT)
  
  if (!transFlag) {
    y.OT <- as.double(x.OT[ind]) - threshold
  } else {
    threshold.trans <- transfun(threshold)
    dth <- threshold.trans - threshold
    y.OT <- transfun(as.double(x.OT[ind])) - threshold.trans
  }
  
  nb.OT <- length(y.OT)
  if(!nb.OT) stop("no data above threshold")

  if(trace) cat("Number of obs > threshold", nb.OT, "\n")
  
  hist.MAX <- MAX$flag
  
  if (hist.MAX) {
    z.MAX <- unlist(MAX$data)
    if ( any(z.MAX <= threshold) ) stop("all historical 'MAX' data must exceed the threshold")
    w.MAX <- MAX$effDuration
    block.MAX <- MAX$block
    r.MAX <- MAX$r
    nblock.MAX <- length(r.MAX)
  } else {
    z.MAX <- numeric(0)
    threshold.max <- numeric(0)
  }
    
  hist.OTS <- OTS$flag
  
  if (hist.OTS) {
    z.OTS <- unlist(OTS$data)
    if ( any(z.OTS <= threshold) ) stop("all historical 'OTS' data must exceed the threshold")
    w.OTS <- OTS$effDuration
    threshold.OTS <- OTS$threshold
    if ( any(threshold.OTS <= threshold) ) stop("all 'OTS' thresholds must exceed the threshold")
    block.OTS <- OTS$block
    r.OTS <- OTS$r
    nblock.OTS <- length(r.OTS)
  } else {
    z.OTS <- numeric(0)
    threshold.OTS <- numeric(0)
  }

  
  ##=================================================================
  ## scale.OT is a rounded quantity used to scale the parameters
  ## during optimisation
  ##
  ## y.OT should be divided by scale.OT. Instead, scale parameters
  ## of distributions should be divided by scale.y 
  ## 
  ##==================================================================
  
  mexpon.OT <- floor(log(quantile(y.OT, prob = 0.75), base = 10))

  if (mexpon.OT >= 2) {
    scale.need <- TRUE
    scale.OT <- 10^mexpon.OT
  } else {
    scale.need = FALSE
    scale.OT <- 1.0
  }
  
  ##=================================================================
  ## Analysis of the chosen distribution for the exceedances
  ## 
  ## distname.y   a code, not necessarily related to the fun names
  ## funname.y    the suffix of called prob. functions
  ## parnames.y   vector of names for ALL the parameters
  ## parnb.y      total number of parameters (fixed parameters
  ##              included if any).  
  ##
  ## parnames.all is c("lambda", parnames.y) since the evt rate
  ##              is called here "lambda".
  ## 
  ## The special  distribution have perdefined parnames. 
  ##
  ##=================================================================

  special <- TRUE
  
  if (distname.y == "exponential") {
    funname.y <- "exp"
    parnames.y <- "rate"
    scale.y <- 1/scale.OT
  } else if (distname.y == "weibull") {
    funname.y <- "weibull"
    parnames.y <- c("shape", "scale")
    scale.y <- c(1, scale.OT)
  } else if (distname.y == "gpd") {
    require(evd)      ## At the time, evd is necessary
    funname.y <- "gpd"
    parnames.y <- c("scale", "shape")
    scale.y <- c(scale.OT, 1)
  } else if (distname.y %in% c("log-normal", "lognormal")) {
    distname.y <- "log-normal"
    funname.y <- "lnorm"
    parnames.y <- c("meanlog", "sdlog")
    scale.y <- c(1, 1)
  } else {
    
    special <- FALSE
    
    ## Will LATER join the "special" group
    if (distname.y == "gamma"){
      funname.y <- "gamma"
      parnames.y <- c("shape", "scale")
      scale.y <- c(1, scale.OT)
      ## Compute initial values
      if ( is.null(start.par.y) ) {
        start.par.y <-
          unlist(mom2par("gamma", mean = mean(y.OT), sd = sd(y.OT)))
      }
    } else if (distname.y %in% c("MixExp2", "mixexp2")) {

      warning("With the 'mixexp2' ML estimation may fail to converge") 
      
      distname.y <- "mixexp2"
      funname.y <- "mixexp2"
      parnames.y <- c("prob1", "rate1", "delta")
      scale.y <- c(1, 1/scale.OT, 1/scale.OT)
      ## compute initial values
      if ( is.null(start.par.y) ) {
         ini <- ini.mixexp2(y.OT, plot = FALSE)$estimate
         start.par.y <- c(ini["prob1"], ini["rate1"], ini["rate2"]- ini["rate1"])
         names(start.par.y) <- parnames.y
      }

    } else{

      ##-------------------------------------------------------------------
      ## The parameters names are obtained by concatenating the
      ## fixed.par.y and start.par.y list in this order (Fixed First)
      ## keeping the order in both source
      ## This order must be respected when storing estimation results
      ##-------------------------------------------------------------------
      
      warning("warning: distribution not in target list. Still EXPERIMENTAL")
      funname.y <- distname.y
      
      test <- intersect(names(fixed.par.y), names(start.par.y))
      if (length(test)>0) stop("parameter name(s) can not be in fixed.par.y AND start.par.y") 
      if ( (length(fixed.par.y)>0) && any(is.na(fixed.par.y)) )
        stop("'fixed.par.y' can not contain NA values") 
      
      parnames.y <- c(names(fixed.par.y), names(start.par.y))
      if (trace) cat("parnames.y=", parnames.y, "\n")
      scale.y <- rep(1, length.out = length(parnames.y))
    }
  }

  names(scale.y) <- parnames.y
  if (trace) {
    cat("scaling for parameters\n")
    print(scale.y)
  }
  
  if (trace) {
    if (special) cat("o The distribution for exceedances is recognized as 'special'\n")
    else cat("o The distribution for exceedances is not 'special'\n")
    cat("  Its parameter names are\n")
    cat("  ", parnames.y, "\n")
  }

  parnames.all <- c("lambda", parnames.y)
  parnb.y <- length(parnames.y)
  
  ##=================================================================
  ## Analysis of provided values for parameters (FIXED parameters)
  ##
  ## 'fixed.y'    is a flag indicating which parameter is fixed.
  ##              It must have the same length and the same order as
  ##              'parnames.y'
  ## 'fixed.all'  is c(FALSE, fixed.y) since the evt rate can not
  ##              be fixed.
  ## 'p.y'        is the number of parameters for y in estimation
  ##              1 <= p.y <= parnb.y.
  ## 'pf.y'       is the number of fixed parameters.
  ##
  ##=================================================================
  
  fixed.par.y <- unlist(fixed.par.y)
  
  m <- match(names(fixed.par.y), table = parnames.y)
  
  if ( any(is.na(m)) ) stop("name(s) not understood in 'fixed.par.y'")
  
  if ( length(fixed.par.y) && any(is.na(fixed.par.y)) )
    stop("fixed.par.y must contain only non-missing values")
  
  fixed.y <- rep(FALSE, length(parnames.y))
  names(fixed.y) <- parnames.y
  fixed.y[m] <- TRUE
  
  ## cat("XXX parnames.y", parnames.y, "\n"); print(fixed.y)
  
  p.y <- length(parnames.y[!fixed.y])
  if(p.y == 0) stop("No parameter to estimate for y!")
  
  if ( any(fixed.y)) {
    pf.y <- sum(fixed.y)
    if (trace) {
      cat("o Fixed parameters for y\n")
      print(fixed.y)
     }
  } else pf.y <- 0
  
  if (trace) {
    cat("o Number of param for exceedances estimated from data\n")
    cat("  p.y = ", p.y, "\n")
  }

  fixed.all <- c(FALSE, fixed.y)
  names(fixed.all) <- parnames.all

  ##=================================================================
  ## compute the degree of freedom
  ## For historical blocks, each observation is considered as a
  ## df, and each empty OTS block must be so.
  ##=================================================================
  
  df <- length(y.OT)
  if (MAX$flag) 
    df <- df + length(z.MAX) 
  if (OTS$flag)  
    df <- df + length(z.OTS) + sum(r.OTS == 0)  
  
  ##=================================================================
  ## Event rate estimation
  ##=================================================================
  
  lambda.hat <- nb.OT / effDuration 
  est.N <- lambda.hat

  ## There wa a bug here previous to 0-5.0!!!
  cov.N <- lambda.hat / effDuration 
  
  names(est.N) <- "lambda"
  names(cov.N) <- "lambda"
  
  if (trace) cat("o Estimate of the process rate (evt/bloc length)", est.N, "\n")
  
  ##====================================================================
  ## Wrapper functions
  ## 
  ## In this part, formals of functions are changed following the
  ## ideas in "fitdistr" of the MASS package. See therein.
  ##
  ## Note: the code could be rearranged using a loop on function names
  ## 
  ##====================================================================
  
  ## reorder arguments to densfun and co
  dfun.y <- get(paste("d", funname.y, sep = ""), mode = "function")
  pfun.y <- get(paste("p", funname.y, sep = ""), mode = "function")
  qfun.y <- get(paste("q", funname.y, sep = ""), mode = "function")
  
  fms <- formals(dfun.y)
  args <- names(fms)
  m <- match(parnames.y, args)  
  if(any(is.na(m)))
    stop("'parnames.y' specifies names which are not arguments to 'densfun'")

  ## 'x' is maintened in pole position '1', then come the params in parnames,
  ## then the other if any (e.g. surrogate parameters)
  formals(dfun.y) <- c(fms[c(1, m)], fms[-c(1, m)])
  
  ## Caution: in function call, remember to use  log = TRUE
  
  logf.y <- function(parm, x) dfun.y(x, parm, log = TRUE)
  
  ## Same thing for 'pfun'. Note that although the main arg of
  ## distributions functions is usually names "q", we force it 
  ## to be "x" here because f.y and F.y are usaed in the same
  ## manner in the log-likelihood!
  
  fms <- formals(pfun.y)
  args <- names(fms)
  m <- match(parnames.y, args)
  if(any(is.na(m)))
    stop("parnames.y specifies names which are not arguments to 'pfun'")
  
  formals(pfun.y) <- c(fms[c(1, m)], fms[-c(1, m)])
  
  F.y <- function(parm, x) pfun.y(x, parm)
  
  ## reorder arguments to densfun
  fms <- formals(qfun.y)
  args <- names(fms)
  m <- match(parnames.y, args)
  if(any(is.na(m)))
    stop("parnames.y specifies names which are not arguments to 'qfun'")
  
  formals(qfun.y) <- c(fms[c(1, m)], fms[-c(1, m)])
  
  q.y <- function(parm, p) qfun.y(p, parm)

  
  ##=================================================================
  ## Hack formals and body for the wrapper functions
  ## The case p.y == 0 (no estimation) is not possible at the time
  ##=================================================================
  
  if (p.y >= 1) {

    ## to remove later!!!
    str <- paste(paste("parm[", 1:p.y, "]", collapse = ", ", sep = ""), ")")
    
    nfn <- parnames.y[!fixed.y]
    
    str <- paste(paste(paste(nfn, paste("parm[\"", nfn, "\"]", sep = ""), sep = " = "),
                       collapse = ", "), sep = "")

  } else {
    if (pf.y == 0)
      stop("no parameter to estimate and no parameter fixed. Check distribution")
  }
    
  ## Add fixed values if necessary
  ## This is done by modifying the body of the functions by calling
  ## the relevant function in it with suitable NAMED args. 
  
  if (pf.y >= 1) {
    ## modif du 2010-01-25 for fixed par
    if (p.y >= 1) {
      strf <- paste(paste(paste(names(fixed.par.y), fixed.par.y, sep = " = "),
                          collapse = ", "), sep = "")
      str <- paste(str, strf, sep = ", ")
    } else {
      str <- paste(paste(paste(names(fixed.par.y), fixed.par.y, sep = " = "),
                         collapse = ", "), sep = "")
    }
    
  } 
  
  body(logf.y) <- parse(text = paste("dfun.y(x,", str,", log = TRUE)") )     
  body(F.y) <- parse(text = paste("pfun.y(x,", str, ")"))
  body(q.y) <- parse(text = paste("qfun.y(p,", str, ")"))
  
  ##==================================================================
  ## a list of functions  to be exported.
  ## Note that the definition of 'dfun.y', 'qfun.y' and 'pfun.y'
  ## is taken from the environment where the funs are defined, i.e.
  ## here (lexical scopting).
  ##
  ## CAUTION
  ##
  ## Do not re-use the functions logf.y, q.y ot F.y in any environment
  ## where the definition of dfun.y, pfun.y or qfun.y could be
  ## different
  ##
  ##==================================================================
  
  funs <- list(transfun = transfun,
               invtransfun = invtransfun,
               dfun.y = dfun.y,
               pfun.y = pfun.y,
               qfun.y = qfun.y,
               logf.y = logf.y,
               q.y = q.y, 
               F.y = F.y)

  ##  modif 2010-01-25 for the fixed parms case
  ## (did not exist before)
  ind.est.y <- (1:parnb.y)[!fixed.y]   
  
  ##=================================================================
  ## ML estimation ignoring historical data (if any). The main goal
  ## is here to compute
  ##
  ## 'est.y'    estimates for th y part including fixed parms if any
  ## 'cov0.y'   covariance matrix for the y parameters (filled with
  ##            zeros for fixed params).
  ##
  ## For some distributions, a special routine may be used.
  ## 
  ##=================================================================

  opt0 <- NULL
  
  if(distname.y == "exponential") {
    
    ## Explicit maximum likelihood
    if (fixed.y) {
      rate.hat <- fixed.par.y
      cov0.y <- matrix(0,
                       nrow = 1, ncol = 1)
    } else {
      rate.hat <- 1 / mean(y.OT)
      est.y <- rate.hat
      
      cov0.y <- matrix(rate.hat^2 / nb.OT,
                       nrow = 1, ncol = 1)
    }
    names(est.y) <- parnames.y
    colnames(cov0.y) <- rownames(cov0.y) <- "rate"
    
  } else if(distname.y == "weibull") {
    
    ## Explicit maximum likelihood
    if (fixed.y["shape"]) {
      
      shape.hat <- fixed.par.y["shape"]
      if (shape.hat < 0) stop("'fixed.par.y' must specify a 'shape' > 0")
 
      if (fixed.y["scale"]) {
        scale.hat <- fixed.par.y["scale"]
        cov0.y <- matrix(0, nrow = 2, ncol = 2)
      } else {
        scale.hat <- mean(y.OT^shape.hat)^{1/shape.hat}
        cov0.y <- matrix(0, nrow = 2, ncol = 2)
        cov0.y[2, 2] <- ( (scale.hat/shape.hat)^2 ) / nb.OT
      }
      
    } else {

      if (fixed.y["scale"]) {
         if (shape.hat < 0) stop("fixed scale in Weibull not allowed in this version")
      } else {
        ## concentrated maximum likelihood
        fit.OT <- fweibull(x = y.OT, info.observed = FALSE)
        est.OT <- fit.OT$estimate
        shape.hat <- est.OT["shape"]
        scale.hat <- est.OT["scale"]
        cov0.y <- fit.OT$cov
      }
      
    }

    ## CAUTION: here respect the order specified above for Weibull!
    est.y <- c(shape.hat, scale.hat)
    names(est.y) <- parnames.y
    colnames(cov0.y) <- rownames(cov0.y) <- parnames.y

  } else if (distname.y =="log-normal") {
    
    ly.OT <- log(y.OT)
    meanlog.hat <- mean(ly.OT)
    sdlog.hat <- sd(ly.OT)
    
    est.y <- c(meanlog.hat, sdlog.hat)
    names(est.y) <- parnames.y
    
    sig2 <- sdlog.hat^2
    ## Do not forget that 
    cov0.y <- matrix(c(sig2/nb.OT, 0, 0, sig2/2/nb.OT), nrow = 2, ncol = 2)
    colnames(cov0.y) <- rownames(cov0.y) <- names(est.y)

  } else if(distname.y == "gpd") {
    
    if (pf.y) stop("Fixed parameters not allowed yet for gpd special distribution")
   
    ## XXX NB: fixed parameters case to be done. Write a fgpd???
    ## Then CAUTION with the parameter names order.
    
    fit.OT <- fpot(x = y.OT, threshold = 0, model = "gpd")
    est.y <- fit.OT$estimate
    cov0.y <-  fit.OT$var.cov
    colnames(cov0.y) <- rownames(cov0.y) <- names(est.y)
    
    scale.hat <- est.y["scale"] 
    shape.hat <- est.y["shape"]
    
  } else {

    optim0 <- TRUE
    
    ## Arbitrary distribution: perform a general ML estimation 
    
    loglik0 <- function(parms) {
      logL <- sum(logf.y(parm = parms, x = y.OT))     
    }
    
    ## let's go...
    if (trace) {
      cat("o Optimization\n")
      cat("  initial values\n")
      print(start.par.y)
    }

    ## Scale the parameters.
    if (!is.null(control$parscale)) {
      pn <- parnames.y[!fixed.y]
      cp <-  scale.y[!fixed.y]
      cpn <- names(control$parscale)
      if (!all(cpn %in% pn))
        stop("'control$parscale' must have names in", pn)
      cp[cpn] <- control$parscale[cpn]
      control$parscale <- cp
    } else {
      control$parscale <- scale.y[!fixed.y]
    }
    
    control$ndeps <- rep(eps, length(control$parscale))
     
    if (suspend.warnings) opt.old <- options(warn = -1)
    
    if (trace) {
      cat("  parscale used in 'control'\n"); print(control$parscale)
      cat("  fnscale used in 'control'\n"); print(control$fnscale)
    }   
    
    opt0 <- optim(par = start.par.y,
                  fn = loglik0,
                  ## method = "BFGS",
                  hessian = !numDeriv,
                  control = control)
    
    if (suspend.warnings) opt.old <- options(opt.old)

    if (opt0$convergence!=0) {
      print(opt0)
      stop("convergence not reached in optimisation")
    }

    if (trace) {
      cat("  optim ended normally.\n")
      print(start.par.y)
    }

    ## Now build est.y and cov0.y
    est.y <- rep(NA, length(parnames.y))
    names(est.y) <- parnames.y

    ## Fill the two parts separately
    est.y[fixed.y] <- fixed.par.y
    est.y[!fixed.y] <- opt0$par
    
    ## modified on 2010-01-25
    cov0.y <- matrix(0, nrow = parnb.y, ncol = parnb.y)
    colnames(cov0.y) <- rownames(cov0.y) <- parnames.y

    ## This was added in versions > 0.5-2 due to
    ## repeted problems with hessian
    if (numDeriv) {
      require(numDeriv)
      opt0Hessian <- numDeriv::hessian(fun = loglik0, x = opt0$par)
    } else {
      opt0Hessian <- opt0$hessian
    }
    
    if (trace) {
      cat("Hessian for ordinary logL\n")
      print(opt0Hessian)
    }
    
    cov.prov0 <- -solve(opt0Hessian)
    eig  <- eigen(cov.prov0, symmetric = TRUE)$values 
    if (any(eig <= 0)) {
      warning("hessian not negative definite (1-st optim). Confidence limits may be misleading")
    }
    cov0.y[ind.est.y, ind.est.y] <- cov.prov0
    
  }

  ## modified on 2010-01-25
  if (trace) {
    cat("o Estimated values / covariance for the exceedances part\n")
    print(est.y)
    print(cov0.y)
  }
  
  param.N <- "lambda"
  
  ##=================================================================
  ## When historical data is present, the global log-likelihood must
  ## be maximised with "optim"
  ##=================================================================
  
  if ( hist.MAX || hist.OTS ) {

    if (hist.MAX) {   ## homogenise with other obs.
      if (!transFlag) zMod.MAX   <- z.MAX - threshold
      else zMod.MAX   <- transfun(z.MAX) - threshold.trans
      
      zrMod.MAX <- tapply(zMod.MAX, block.MAX, min)

      if (trace) 
        cat("\nTake into account MAX historical data of time-length", sum(w.MAX), "units\n")
    }
    
    if (hist.OTS) {   ## homogenise with other obs. 
      if (!transFlag) {
        zMod.OTS   <- z.OTS - threshold
        thresholdMod.OTS <- threshold.OTS - threshold
      } else {
        zMod.OTS   <- transfun(z.OTS) - threshold.trans
        thresholdMod.OTS <- transfun(threshold.OTS) - threshold.trans
      }
      
      if (trace) 
        cat("\nTake into account OTS historical data on time-length", sum(w.OTS), "units \n")
    }

    ##-----------------------------------------------------------------
    ## General Log-Likelihood: its formal args are obtained by
    ## concatenating "lambda" and the vector of parameters for the "y"
    ## part (exceedances).
    ##-----------------------------------------------------------------
    
    loglik <- function(parms) {

      lambda <- parms[1]  
      lw <- lambda*effDuration
      
      ## Ordinary OT part caution all blocks OT are grouped as one 
      logL  <- dpois(nb.OT, lambda = lw, log = TRUE)
      logL  <- logL + sum(logf.y(parm = parms[-1], x = y.OT))     

      if (hist.MAX) {
        lw.MAX <- lambda*w.MAX   
        logL <- logL + sum(r.MAX*log(lw.MAX)) - sum(lw.MAX * ( 1 - F.y(x = zrMod.MAX, parm = parms[-1]) ))
        logL <- logL + sum(logf.y(parm = parms[-1], x = zMod.MAX))
      }

      if (hist.OTS) {
        lw.OTS <- lambda * w.OTS * (1.0 - F.y(x = thresholdMod.OTS, parm = parms[-1]))
        logL <- logL - sum(lw.OTS)
        if (sum(r.OTS) > 0) logL <- logL + sum(r.OTS *log(lw.OTS)) + sum(logf.y(parm = parms[-1], x = zMod.OTS))
        ## cat("LogL =", logL, "u.OTS = ", thresholdMod.OTS, "mod = ", -sum(lw.OTS), "\n")
      }
        
      logL

    }
    
    ## let's go...
    if (trace) cat("o Optimisation\n")
    if (suspend.warnings) opt.old <- options(warn = -1)

    if (force.start.H) {

      if (is.null(start.par.y))
        stop("You must provide initial values in 'start.par.y' when 'force.start.H' is TRUE")

      par.ini <- c(est.N, as.numeric(start.par.y))

    } else {

      ## use the estimation without historical data. Note that the elts need to be named 
      par.ini.all <- c(est.N, est.y)
      names(par.ini.all) <- parnames.all
      par.ini <- par.ini.all[!fixed.all]
      
      if ( is.na(loglik(par.ini)) || !is.finite(loglik(par.ini)) ) {
        
        ## GPD "Weibullian" case for gpd. "doctorize" the identified problem
        ## The initial values may fail to be acceptable
        
        if ( (distname.y == "gpd") && (par.ini["shape"] < 0) ) {

          ## if one of the historical levels  is over the quantile with prob. 0.995
          ## of the distribution, try to give a larger value to the scale parameter. 
          
          warning(paste("initial parameters for the GPD distribution",
                        "lead to some values outside of support. These",
                        "values are modified"))
          
          upmax <- q.y(parm = est.y, p = 0.999)
          
          ## known doctorizable situations
          if ( (hist.MAX || hist.OTS) && any(c(z.MAX, z.OTS, threshold.OTS) > upmax + threshold) ) {
            if (trace) {
              mmax <- max(c(z.MAX, z.OTS, threshold.OTS))
              cat("mmax = ", mmax, "\n")
              
              uu <- threshold - par.ini["scale"] / par.ini["shape"]
              cat("   Upper limit of the support (estimation phase 1)", uu,"\n")
            }

            if (trace)
              cat("   Old par.ini[\"scale\"]", par.ini["scale"],"\n")
            
            par.ini["scale"] <- par.ini["scale"] * ( max(c(z.MAX, z.OTS, threshold.OTS)) - threshold ) / upmax
            
            if (trace)
              cat("   New par.ini[\"scale\"]", par.ini["scale"],"\n")
            
          }
          ## else if (hist.OTS && any(z.OTS > upmax + threshold) ) {
          ##  par.ini["scale"] <- par.ini["scale"] * (max(z.OTS) - threshold) / upmax 
          ##}

          if (trace) {
            cat("   Initial values 2nd stage (modified)\n"); print(par.ini)
          }
          
        } else {
          stop(paste("loglik is NA or infinite at the initial values.",
                     "Give correct initial values in 'start.par.y' and set 'force.start.H' to TRUE"))    
        }
        
      }
      
    }

    if (trace) {
      cat("o Initial values for 2nd stage optimization\n")
      print(par.ini)
    }
    
    ## Scale the parameters.
    if (!is.null(control.H$parscale)) {
      pn <- c("est.N", parnames.y[!fixed.y])
      cp <- c(est.N, scale.y[!fixed.y])
      cpn <- names(control.H$parscale)
      if (!all(cpn %in% pn))
        stop("'control$parscale' must have names in", pn)
      cp[cpn] <- control.H$parscale[cpn]
      control.H$parscale <- cp
    } else {
      control.H$parscale <- c(est.N, scale.y[!fixed.y])
    }
    
    control.H$ndeps <- rep(eps, length(control.H$parscale))
     
    if (trace) {
      cat("  parscale used in 'control.H'\n")
      print(control.H$parscale)
      cat("  ndeps used in 'control.H'\n")
      print(control.H$ndeps)
      cat("  fnscale used in 'control.H'\n")
      print(control.H$fnscale)
    }   
    
    opt <- optim(par = par.ini,
                 fn = loglik,
                 hessian = !numDeriv,
                 control = control.H)           

    if (suspend.warnings) opt.old <- options(opt.old)

    if (opt$convergence != 0) {
      print(opt) 
      stop("convergence not reached in optimisation")
    }

    ## XXX 
    estimate <- rep(NA, parnb.y + 1)
    names(estimate) <- parnames.all
    estimate[fixed.all] <- fixed.par.y
    estimate[!fixed.all] <- opt$par

    if (numDeriv) {
      require(numDeriv)
      optHessian <- numDeriv::hessian(fun = loglik, x = opt$par) 
    } else {
      optHessian <- opt$hessian
    }
      
    if (trace) {
      cat("   Fixed parameters in historical optimization and hessian\n")
      print(fixed.all)
      cat("   Hessian\n")
      print(optHessian)
    }
    
    ## XXX Caution: only parameters which are not fixed.
    cov.all <- matrix(0, nrow = parnb.y + 1, ncol = parnb.y + 1)
    colnames(cov.all) <- rownames(cov.all) <- parnames.all
    cov.prov <- -solve(optHessian)
    
    eig  <- eigen(cov.prov, symmetric = TRUE)$values 
    if (any(eig <= 0)) {
      warning("hessian not negative definite (2-nd optim). Confidence limits may be misleading")
    }
    cov.all[!fixed.all, !fixed.all] <- cov.prov
    
    res <- list(call = mc,
                x.OT = x.OT,
                y.OT = y.OT,
                effDuration = effDuration,
                threshold = threshold,
                distname.y = distname.y,
                p.y = p.y,
                parnames.y = parnames.y,
                fixed.y = fixed.y,
                trans.y = trans.y,
                est.N = est.N,
                cov.N = cov.N,
                est.y = est.y,
                cov.y = cov0.y,
                corr.y = cov2corr(cov0.y),
                estimate = estimate,
                fixed = fixed.all,
                df = df,
                p = p.y + 1,
                opt0 = opt0,
                opt = opt,
                sigma = sqrt(diag(cov.all)),
                cov = cov.all,
                corr = cov2corr(cov.all),
                history.MAX = MAX,
                history.OTS = OTS,
                transFlag = transFlag,
                funs = funs)
    
    ## For later use...
    est.y <- estimate[-1]
    names(est.y) <- parnames.y
    
    ind <- 1 + ((1L):(parnb.y))
    cov.y <- cov.all[ind, ind, drop = FALSE]
    
  }  else {
    
    ##=============================================================
    ## No historical data
    ##=============================================================
    
    estimate <- c(est.N, est.y)
    
    ## modif 2010-01-25 for the fixed parms case
    ## cov.all <- matrix(0, nrow = p.y+1, ncol = p.y+1)
    cov.all <- matrix(0, nrow = parnb.y + 1, ncol = parnb.y + 1)
    colnames(cov.all) <- rownames(cov.all) <- c("lambda", parnames.y)
    
    cov.all[1, 1] <- cov.N
   
    ind <- (1L):parnb.y
    
    cov.all[1+ind, 1+ind] <- cov0.y[ind, ind, drop = FALSE]    
    cov.y <- cov0.y
    
    res <- list(call = mc,
                x.OT = x.OT,
                y.OT = y.OT,
                effDuration = effDuration,
                threshold = threshold,
                distname.y = distname.y,
                p.y = p.y,
                parnames.y = parnames.y,
                fixed.y = fixed.y,
                trans.y = trans.y,
                est.N = est.N,
                cov.N = cov.N,
                est.y = est.y,
                cov.y = cov0.y,
                corr.y = cov2corr(cov0.y),
                estimate = estimate,
                fixed = fixed.all,
                df = df,
                p = p.y + 1,
                sigma = sqrt(diag(cov.all)),
                cov = cov.all,
                corr = cov2corr(cov.all),
                history.MAX = MAX,
                history.OTS = OTS,
                funs = funs,
                transFlag = transFlag)
  }

  ##=============================================================
  ## Inference on return levels using the delta method
  ##
  ## Prepare a matrix object ret.lev
  ##
  ##=============================================================

  rl.period <- 1 / estimate[1] / (1 - rl.prob)
  
  ## restrict 'pred.period' if necessary
  pred.prob <- 1.0 - 1.0 / estimate[1] / pred.period
  ind <- (pred.prob > 0.0) & (pred.prob < 1.0)
  if (any(!ind)) warning("some return periods in 'pred.period' out of range")

  pred.period <- pred.period[ind]
  pred.prob <- pred.prob[ind]
  
  rl.sort <- sort(c(rl.period, pred.period), index.return = TRUE)
  rl.period <- rl.sort$x
  
  rl.prob <- c(rl.prob, pred.prob)
  rl.prob <- rl.prob[rl.sort$ix]

  ## remove dupplicates
  ind <- !duplicated(rl.period)
  rl.period <- rl.period[ind]
  rl.prob <- rl.prob[ind]
  
  ind.pred <- rl.period %in% pred.period
  
  nc <- 3 + 2*length(pct.conf)
  cnames <-
    c("period", "prob", "quant",
      paste(rep(c("L", "U"), length(pct.conf)), rep(pct.conf, each = 2), sep = "."))
  
  ## Return levels (e.g. to be used in RLplot)
  ret.lev <- matrix(NA, nrow = length(rl.period), ncol = nc)

  ## 2012-12-22 comment out ret.lev
  ## rownames(ret.lev) <- format(rl.period)
  colnames(ret.lev) <- cnames

  ret.lev[ , "period"] <- 1 / estimate[1] /(1 - rl.prob)
  ret.lev[ , "prob"] <- rl.prob
  ret.lev[ , "quant"] <- threshold + q.y(parm = estimate[-1], p = rl.prob)
  
  if (transFlag) {
    ret.lev[ , "quant"] <- invtransfun(dth + ret.lev[ , "quant"])
  }
  
  if ( (distname.y == "exponential") && !hist.MAX && !hist.OTS ) {

    if (trace) cat("Special inference for the exponential case without history\n")
    
    ##--------------------------------------------------------------
    ## Use the sampling distribution to derive confidence
    ## limits on quantiles. If lambda was known, these limits
    ## would be exact.
    ##--------------------------------------------------------------
    
    for (ipct in 1:length(pct.conf)) {
      
      alpha.conf <- (100 - pct.conf[ipct])/100
      
      theta.L <- 2*nb.OT / est.y / qchisq(1 - alpha.conf/2, df = 2 * nb.OT)
      theta.U <- 2*nb.OT / est.y / qchisq(alpha.conf/2, df = 2 * nb.OT)
      
      ret.lev[ , 2*ipct + 2] <- threshold + qexp(p = rl.prob, rate = 1/theta.L)
      ret.lev[ , 2*ipct + 3] <- threshold + qexp(p = rl.prob, rate = 1/theta.U)
      
    }
    
    if (transFlag) {

      dth <- threshold.trans - threshold
      
      for (ipct in 1:length(pct.conf)) {
        
        ret.lev[ , 2*ipct + 2] <- invtransfun(dth + ret.lev[ , 2*ipct + 2])
        ret.lev[ , 2*ipct + 3] <- invtransfun(dth + ret.lev[ , 2*ipct + 3])

      }
      
    }
    
    ##=============================================================
    ## perform Bartlett's gof test
    ##==============================================================

    res$infer.method <- "chi-square for exponential distribution (no historical data)"
    
  } else {
 
    ##==============================================================
    ## DELTA METHOD
    ## matrices for numerical derivation
    ## The column ip of Parmat contains the parameter value with
    ## a tiny modification of its ip component.
    ##==============================================================

    parmMat <- matrix(est.y[!fixed.y], nrow = p.y, ncol = p.y)
    
    rownames(parmMat) <- parnames.y[!fixed.y]
    colnames(parmMat) <- parnames.y[!fixed.y]
    
    dparms <- abs(est.y[!fixed.y])*eps
    dparms[dparms < eps] <- eps
    
    for (ip in 1:p.y) parmMat[ip, ip] <- parmMat[ip, ip] + dparms[ip]
    
    delta <- rep(NA, p.y)
    sig <- rep(NA, length(rl.prob))
        
    ##=============================================================
    ## Compute the delta's  and sig's
    ## delta conains derivative w.r.t. unknown params
    ##
    ## CAUTION
    ##
    ## The quantile function q.y takes a vector as first arg which
    ## should morally be of length p.y but is here of length
    ## parnb.y. This works because q.y uses the elements in named
    ## form irrespective of their position.
    ##
    ##=============================================================
    
    for (i in 1:length(rl.prob)) {
      
      for (ip in 1:p.y) {

        est.prov <- est.y
        est.prov[!fixed.y] <-  parmMat[ , ip]
        
        dd <- ( q.y(est.prov, p = rl.prob[i]) - q.y(est.y, p = rl.prob[i]) ) / dparms[ip]
        delta[ip] <- dd
        
      }
      
      sig[i] <- sqrt(t(delta)%*%cov.y[!fixed.y, !fixed.y]%*%delta)

    }
    
    ##=============================================================
    ## Apply "delta method" for the quantiles
    ##
    ## Note that q.y MUST accept a vectorized prob 
    ## Use a loop on i to remove this constraint ???
    ##==============================================================
      
    for (ipct in 1:length(pct.conf)) {
      alpha.conf <- (100 - pct.conf[ipct])/100
      z.conf <- qnorm(1 - alpha.conf/2)
      ret.lev[ , 2*ipct + 2] <- threshold + q.y(est.y, p = rl.prob) - z.conf * sig
      ret.lev[ , 2*ipct + 3] <- threshold + q.y(est.y, p = rl.prob) + z.conf * sig
    }
    
    if (transFlag) {
      for (ipct in 1:length(pct.conf)) {
        ret.lev[ , 2*ipct + 2] <- invtransfun(dth + ret.lev[ , 2*ipct + 2])
        ret.lev[ , 2*ipct + 3] <- invtransfun(dth + ret.lev[ , 2*ipct + 3])
      }
    }
      
    res$infer.method <- "delta-method with numerical derivative"
    
  }

  ##======================================================================
  ## Add ret.lev and pred to the list of returned items
  ##======================================================================
  
  res$pct.conf <- pct.conf
  
  ret.lev <- as.data.frame(ret.lev)
  res[["ret.lev"]] <- ret.lev

  pred <- as.data.frame(ret.lev[ind.pred, , drop = FALSE])

  ## change rownames to have integers, which is generally not the case for ret.lev!
  rownames(pred) <- format(pred$period)
    
  res[["pred"]] <- pred
  
  ##======================================================================
  ## perform a Kolmogorov Smirnov test Note that a mix positional
  ## matching and name matching in the call!!! This is because
  ## F.y has parm as first arg
  ##======================================================================

  if (jitter.KS) {
    KS <- ks.test(OTjitter(y.OT, threshold = 0.0),
                  F.y, parm = estimate[-1])
  } else KS <- ks.test(y.OT, F.y, parm = estimate[-1])
  
  res$KS.test <- KS

  if (distname.y == "exponential")  res$expon.test <- gofExp.test(x = y.OT)

  ##======================================================================
  ## Return level plot
  ##======================================================================

  class(res) <- "Renouv"

  if (plot) plot(res, ...)

  return(res)
  
}
