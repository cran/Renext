
##====================================================================
## Prefer corr and sd to covariance?
##====================================================================

"cov2corr" <- function(cov) {
  s <- sqrt(diag(cov))
  res <- cov / outer(X = s, Y = s, FUN = "*")
  res
}

##====================================================================
## Author: Y. Deville
##
## The black-box maximization of the log-likelhood is inspired
## fitdistr of B. Ripley or 
##
##====================================================================

"fRenouv" <- function(x.OT,
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
                      ...) {

  ## for numerical differentiation
  eps <- sqrt(.Machine$double.eps)
  eps <- 1e-6
  
  ##=================================================================
  ## prepare transforms if wanted/possible
  ##=================================================================
  
  if (!is.null(trans.y)) {
    if( !is.character(trans.y) || !(trans.y %in% c("square", "log")) ) 
      stop("trans.y must be NULL or be character in c( \"square\', \"log\")")
    else if (distname.y != "exponential") {
      stop("non null value for 'trans.y' is only allowed when distname.y == \"exponential\"") 
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
    ## transfun <- function(x) x 
    ## invtransfun <- function(x) x 
  }


  ##=================================================================
  ## default prob for quantile and confidence lims 
  ## should densify near 0 and 1
  ##=================================================================
  
  if (is.null(prob)) {
    prob <- c(0.0001,
              seq(from = 0.01, to = 0.09, by = 0.01),
              seq(from = 0.10, to = 0.80, by = 0.10),
              seq(from = 0.85, to = 0.99, by = 0.01),
              0.995, 0.996, 0.997, 0.998, 0.999,
              0.9995, 0.9996, 0.9997, 0.9998, 0.9999,
              0.99999, 0.999999)
    prob <- prob[prob <= prob.max]
  } else {
    if (any(is.na(prob))) stop("'prob' values can not be NA") 
    if ( any(prob <= 0) || any(prob >= 1) ) stop("'prob' values must be >0 and <1") 
    prob <- sort(prob)
  }

  if (is.null(pred.period)) {
    rr <- round(log(sumw.BOT)/log(10))
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

  if(trace) cat("Number of obs > threshod", nb.OT, "\n")

  history <- length(z.H)>0
  
  if ( history && any(z.H <= threshold) ) stop("all historical data must exceed the threshold")

  if (history) {
    block.H <- as.factor(block.H)
    if (nlevels(block.H) != length(w.BH)) stop("w.BH must be of length nlevels(block.H)")
    r.BH <- table(block.H)
    if (trace) {
      cat("Historical data\n")
      cat("   number of blocks", nlevels(block.H), "\n")
      tapply(z.H, block.H, print)
    }   
  }
  
  Udata <- length(x.U)>0

  if (Udata) {
    if ( any(is.na(x.U)) || any(!is.finite(x.U)) ) {
      stop("x.U must contain only non NA finite values")
    }
    if ( (length(x.U) != length(w.U)) || any(w.U <= 0) ) {
      stop("w.U must have the same length as x.U  and contain only positive values")
    }
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
  
  ## cat("XXX parnames.y", parnames.y, "\n")
  ## print(fixed.y)
  
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
  ## Event rate estimation
  ##=================================================================
  
  lambda.hat <- nb.OT / sumw.BOT 
  est.N <- lambda.hat

  ## There wa a bug here previous to 0-5.0!!!
  cov.N <- lambda.hat / sumw.BOT 
  
  names(est.N) <- "lambda"
  names(cov.N) <- "lambda"
  
  if (trace) cat("o Estimate of the process rate (evt/bloc length)", est.N, "\n")
  
  ##====================================================================
  ## Wrapper functions
  ## 
  ## In this part, formals of functions are changed following the
  ## ideas in "fitdistr" of the MASS package. See therein.
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
  
  ## Atention: in function call, remember to use  log = TRUE
  
  logf.y <- function(parm, x) dfun.y(x, parm, log = TRUE)
  
  ## Same thing for pfun. Note that although the main arg of
  ## ditributions functions is usually names "q", we force it 
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
  ##
  ## The cas p.y == 0 (no estimation) is not possible at the time
  ##
  ##=================================================================
  
  if(p.y >= 1) {

    ## to remove later!!!
    str <- paste(paste("parm[", 1:p.y, "]", collapse = ", ", sep = ""), ")")
    
    nfn <- parnames.y[!fixed.y]

    ## Previously (before 0.5 ...)
    ## str <- paste(paste(paste(nfn, paste("parm[", 1:p.y, "]", sep = ""), sep = " = "),
    ##                   collapse = ", "), sep = "")
    str <- paste(paste(paste(nfn, paste("parm[\"", nfn, "\"]", sep = ""), sep = " = "),
                       collapse = ", "), sep = "")
      
    ## Add fixed values if necessary
    ## This is done by modifying the body of the functions by calling
    ## the relevant function in it with suitable NAMED args. 
    
    if (pf.y) {
      ## modif du 2010-01-25 for fixed par
      strf <- paste(paste(paste(names(fixed.par.y), fixed.par.y, sep = " = "),
                          collapse = ", "), sep = "")
      str <- paste(str, strf, sep = ", ")
    }
    
    body(logf.y) <- parse(text = paste("dfun.y(x,", str,", log = TRUE)") )     
    body(F.y) <- parse(text = paste("pfun.y(x,", str, ")"))
    body(q.y) <- parse(text = paste("qfun.y(p,", str, ")"))
    
  }
  
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

    ## Arbitrary distribution: 
    ## Perform a general maximum-likelihood estimation 
    
    loglik0 <- function(parms) {
      ## cat(parms, "\n")
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
      cat("  parscale used in 'control'\n")
      print(control$parscale)
      cat("  fnscale used in 'control'\n")
      print(control$fnscale)
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
    if (any(eig < 0)) warning("hessian matrix not negative definite first optimization")
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
  
  if ( history || Udata ) {

    if (history) {
      ## homogenise with other obs.
      if (!transFlag) z.HOT   <- z.H - threshold
      else z.HOT   <- transfun(z.H) - threshold.trans
      
      zr.BHOT <- tapply(z.HOT, block.H, min)
      if (trace) {
        cat("\n")
        cat("Take into account historical data of time-length", sum(w.BH), "units\n")
        cat("====================================================================\n")
      }
    }
    if (Udata) {
      ## homogenise with other obs.
      if (!transFlag) x.UOT   <- x.U - threshold
      else x.UOT   <- transfun(x.U) - threshold.trans
      
      if (trace) {
        cat("\n")
        cat("Take into account unoberved high levels on time-length", sum(w.U), "units \n")
        cat("===========================================================================\n")
      }
    }

    ##-----------------------------------------------------------------
    ## General Log-Likelihood: its formal args are obtained by
    ## concatenating "lambda" and the vector of parameters for the "y"
    ## part (exceedances).
    ##-----------------------------------------------------------------
    
    loglik <- function(parms) {

      lambda <- parms[1]  
      lw.BOT <- lambda*sumw.BOT
      
      ## Ordinary OT part caution all blocks OT are
      ## grouped as one in versions >= 0.3.0
      logL  <- dpois(nb.OT, lambda = lw.BOT, log = TRUE)
      logL  <- logL + sum(logf.y(parm = parms[-1], x = y.OT))     

      if (history) {
        ## Historical part
        lw.BH <- lambda*w.BH
        logL <- logL + sum(r.BH*log(lw.BH)) - sum(lw.BH * ( 1 - F.y(x = zr.BHOT, parm = parms[-1]) ))
        logL <- logL + sum(logf.y(parm = parms[-1], x = z.HOT))
        ## cat("H = ", zr.BHOT, F.y(x = zr.BHOT, parm = parms[-1]), "\n")
      }

      if (Udata) {
        ## Unobserved part
        lw.U <- lambda*w.U
        logL <- logL - sum(lw.U * ( 1 - F.y(x = x.UOT, parm = parms[-1]) ))
        ## cat("U = ", x.UOT, F.y(x = x.UOT, parm = parms[-1]), "\n")
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

      ## use the estimation without historical data. Note that the
      ## elemnts need to be named 
      par.ini.all <- c(est.N, est.y)
      names(par.ini.all) <- parnames.all
      par.ini <- par.ini.all[!fixed.all]
      
      if ( is.na(loglik(par.ini)) || !is.finite(loglik(par.ini)) ) {
        
        ## GPD "Weibullian" case for gpd. "doctorize" the identified problem
        ## The initial values may fail to be acceptable
        
        if ( (distname.y == "gpd") && (par.ini["shape"] < 0) ) {

          ## if one of the z.H is over the quantile with prob. 0.995 of the distribution,
          ## try to give a larger value to the scale parameter. 
          ## if (any(z.H) > q
          
          warning(paste("initial parameters for the GPD distribution",
                        "lead to some values outside of support. These",
                        "values are modified"))
          
          upmax <- q.y(parm = est.y, p = 0.999)
          
          ## known doctorizable situations
          if ( history && any(z.H > upmax + threshold) ) {
            if (trace) {
              uu <- threshold - par.ini["scale"] / par.ini["shape"]
              cat("   Upper limit of the support (estimation phase 1)", uu,"\n")
            }
            
            par.ini["scale"] <- par.ini["scale"] * ( max(z.H) - threshold ) / upmax
            
            if (trace)
              cat("   New par.ini[\"scale\"]", par.ini["scale"],"\n")
            
          } else if ( Udata && any(x.U > upmax + threshold) ) {
            par.ini["scale"] <- par.ini["scale"] * (max(x.U) - threshold) / upmax 
          }

          if (trace) {
            cat("   Initial values 2nd stage (modified)\n")
            print(par.ini)
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
      cat("Fixed parameters in historical optimization and hessian\n")
      print(fixed.all)
      print(optHessian)
    }
    
    ## XXX Caution: only parameters which are not fixed.
    cov.all <- matrix(0, nrow = parnb.y + 1, ncol = parnb.y + 1)
    colnames(cov.all) <- rownames(cov.all) <- parnames.all
    cov.prov <- -solve(optHessian)
    
    eig  <- eigen(cov.prov, symmetric = TRUE)$values 
    if (any(eig < 0)) warning("hessian matrix not negative definite 2nd optimization")
      
    cov.all[!fixed.all, !fixed.all] <- cov.prov
    
    res <- list(y.OT = y.OT,
                threshold = threshold,
                distname.y = distname.y,
                trans.y = trans.y,
                est.N = est.N,
                cov.N = cov.N,
                est.y = est.y,
                cov.y = cov0.y,
                corr.y = cov2corr(cov0.y),
                estimate = estimate,
                opt0 = opt0,
                opt = opt,
                sigma = sqrt(diag(cov.all)),
                cov = cov.all,
                corr = cov2corr(cov.all))
    
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
    
    res <- list(y.OT = y.OT,
                threshold = threshold,
                distname.y = distname.y,
                trans.y = trans.y,
                est.N = est.N,
                cov.N = cov.N,
                est.y = est.y,
                cov.y = cov0.y,
                corr.y = cov2corr(cov0.y),
                estimate = estimate,
                sigma = sqrt(diag(cov.all)),
                cov = cov.all,
                corr = cov2corr(cov.all))
  }

  ##=============================================================
  ## Inference on return levels using the delta method
  ##
  ## Prepare a matrix object ret.lev
  ##
  ##=============================================================
  
  nc <- 3 + 2*length(conf.pct)
  cnames <-
    c("prob", "period", "quant",
      paste(rep(c("L", "U"), length(conf.pct)), rep(conf.pct, each = 2), sep = "."))


  ## Return levels (e.g. to be used in RLplot)
  ret.lev <- matrix(NA, nrow = length(prob), ncol = nc)
  rownames(ret.lev) <- prob
  colnames(ret.lev) <- cnames
  
  ret.lev[ , "prob"] <- prob
  ret.lev[ , "period"] <- 1 / estimate[1] /(1 - prob)
  ret.lev[ , "quant"] <- threshold + q.y(parm = estimate[-1], p = prob)
    
  ## Prepare a 'pred' matrix 
  ## Return levels (e.g. to be used in RLplot)
  
  cnames <-
    c("period", "prob", "quant",
      paste(rep(c("L", "U"), length(conf.pct)), rep(conf.pct, each = 2), sep = "."))

  pred.prob <- 1 - 1/estimate[1]/pred.period
  ind <- (pred.prob>0) & (pred.prob <1)
  pred.period <- pred.period[ind]
  pred.prob <- pred.prob[ind]
  
  pred <- matrix(NA, nrow = length(pred.period), ncol = nc)
  rownames(pred) <- pred.period
  colnames(pred) <- cnames
  
  pred[ , "period"] <- pred.period
  pred[ , "prob"]  <- pred.prob
  pred[ , "quant"] <- threshold + q.y(parm = estimate[-1], p = pred.prob)

  if (transFlag) {
    ret.lev[ , "quant"] <- invtransfun(dth + ret.lev[ , "quant"])
    pred[ , "quant"] <- invtransfun(dth +  pred[ , "quant"])
  } 
  
  if ( (distname.y == "exponential") && (!history) && (!Udata) ) {

    if (trace) cat("Special inference for the exponential case without history\n")
    
    ##--------------------------------------------------------------
    ## Use the sampling distribution to derive confidence
    ## limits on quantiles. If lambda was known, these limits
    ## would be exact.
    ##--------------------------------------------------------------
    
    for (ipct in 1:length(conf.pct)) {
      
      alpha.conf <- (100 - conf.pct[ipct])/100
      
      theta.L <- 2*nb.OT / est.y / qchisq(1 - alpha.conf/2, df = 2 * nb.OT)
      theta.U <- 2*nb.OT / est.y / qchisq(alpha.conf/2, df = 2 * nb.OT)
      
      ret.lev[ , 2*ipct + 2] <- threshold + qexp(p = prob, rate = 1/theta.L)
      ret.lev[ , 2*ipct + 3] <- threshold + qexp(p = prob, rate = 1/theta.U)
      
      pred[ , 2*ipct + 2] <- threshold + qexp(p = pred.prob, rate = 1/theta.L)
      pred[ , 2*ipct + 3] <- threshold + qexp(p = pred.prob, rate = 1/theta.U) 
      
    }
    
    if (transFlag) {

      dth <- threshold.trans - threshold
      
      for (ipct in 1:length(conf.pct)) {

        ret.lev[ , 2*ipct + 2] <- invtransfun(dth + ret.lev[ , 2*ipct + 2])
        ret.lev[ , 2*ipct + 3] <- invtransfun(dth + ret.lev[ , 2*ipct + 3])
        
        pred[ , 2*ipct + 2] <- invtransfun(dth + pred[ , 2*ipct + 2])
        pred[ , 2*ipct + 3] <- invtransfun(dth + pred[ , 2*ipct + 3]) 

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

    ## modif 2010-01-25 for the fixed parms case
    parmMat <- matrix(est.y[!fixed.y], nrow = p.y, ncol = p.y)
    
    ## rownames(parmMat) <- parnames.y
    ## colnames(parmMat) <- parnames.y
    
    rownames(parmMat) <- parnames.y[!fixed.y]
    colnames(parmMat) <- parnames.y[!fixed.y]

    ##  modif 2010-01-25 for the fixed parms case
    ## (ligne inexistante avant)
    ## eps <- sqrt(.Machine$double.eps)
    
    dparms <- abs(est.y[!fixed.y])*eps
    dparms[dparms < eps] <- eps
    
    for (ip in 1:p.y) parmMat[ip, ip] <- parmMat[ip, ip] + dparms[ip]
    
    delta <- rep(NA, p.y)
    sig <- rep(NA, length(prob))
    
    
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
    
    for (i in 1:length(prob)) {
      
      for (ip in 1:p.y) {
        
        ##  modif 2010-01-25 for the fixed parms case 
        ## dd <- ( q.y(parmMat[ , ip], p = prob[i]) -
        ##       q.y(est.y, p = prob[i]) ) / dparms[ip]

        est.prov <- est.y
        est.prov[!fixed.y] <-  parmMat[ , ip]
        
        dd <- ( q.y(est.prov, p = prob[i]) -
               q.y(est.y, p = prob[i]) ) / dparms[ip]
        
        delta[ip] <- dd
        
      }
      
      ## modif 2010-01-25 for the fixed parms case
      ## sig[i] <- sqrt(t(delta)%*%cov.y%*%delta)
       sig[i] <- sqrt(t(delta)%*%cov.y[!fixed.y, !fixed.y]%*%delta)
    }
    
    ##=============================================================
    ## Apply "delta method" for the quantiles
    ##
    ## Note that q.y MUST accept a vectorized prob 
    ## Use a loop on i to remove this constraint ???
    ##==============================================================
      
    for (ipct in 1:length(conf.pct)) {
      alpha.conf <- (100 - conf.pct[ipct])/100
      z.conf <- qnorm(1 - alpha.conf/2)
      ret.lev[ , 2*ipct + 2] <- threshold + q.y(est.y, p = prob) - z.conf * sig
      ret.lev[ , 2*ipct + 3] <- threshold + q.y(est.y, p = prob) + z.conf * sig
    }
    
    if (transFlag) {
      for (ipct in 1:length(conf.pct)) {
        ret.lev[ , 2*ipct + 2] <- invtransfun(dth + ret.lev[ , 2*ipct + 2])
        ret.lev[ , 2*ipct + 3] <- invtransfun(dth + ret.lev[ , 2*ipct + 3])
      }
    }
      
    ##=============================================================
    ## Compute the delta's  and sig's
    ## delta contains derivative w.r.t. unknown params
    ##=============================================================

    pred.sig <- rep(NA, length(pred.prob))
    
    for (i in 1:length(pred.prob)) {
      for (ip in 1:p.y) {
        est.prov <- est.y
        est.prov[!fixed.y] <-  parmMat[ , ip]
        dd <- ( q.y(est.prov, p = pred.prob[i]) -
               q.y(est.y, p = pred.prob[i]) ) / dparms[ip]
        delta[ip] <- dd
      }
      ## modif 2010-01-25 for the fixed parms case
      ##pred.sig[i] <- sqrt(t(delta)%*%cov.y%*%delta)
      pred.sig[i] <- sqrt(t(delta)%*%cov.y[!fixed.y, !fixed.y]%*%delta)
    }
    
    ##=============================================================
    ## Apply "delta method" for the quantiles
    ##
    ## Note that q.y MUST accept a vectorized prob 
    ## Use a loop on i to remove this constraint ???
    ##==============================================================
    
    
    for (ipct in 1:length(conf.pct)) {
      alpha.conf <- (100 - conf.pct[ipct])/100
      z.conf <- qnorm(1 - alpha.conf/2)
      pred[ , 2*ipct + 2] <- threshold + q.y(est.y, p = pred.prob) - z.conf * pred.sig
      pred[ , 2*ipct + 3] <- threshold + q.y(est.y, p = pred.prob) + z.conf * pred.sig
    }

    if (transFlag) {
      for (ipct in 1:length(conf.pct)) {
        pred[ , 2*ipct + 2] <- invtransfun(dth  + pred[ , 2*ipct + 2])
        pred[ , 2*ipct + 3] <- invtransfun(dth  + pred[ , 2*ipct + 3])
      }
    }

    res$infer.method <- "delta-method with numerical derivative"
    
  }

  ret.lev <- as.data.frame(ret.lev)
  res[["ret.lev"]] <- ret.lev

  pred <- as.data.frame(pred)
  res[["pred"]] <- pred
  
  ##======================================================================
  ## perform a Kolmogorov Smirnov test Note that a mix positional
  ## matching and name matching in the call!!! This is because
  ## F.y has parm as first arg
  ##======================================================================
  
  KS <- ks.test(y.OT, F.y, parm = estimate[-1])
  res$KS.test <- KS

  if (distname.y == "exponential")  res$expon.test <- gofExp.test(x = y.OT)

  ##======================================================================
  ## Return level plot
  ##======================================================================
  
  if (plot) {
    
    RLplot(data = ret.lev,
           x = sort(x.OT),
           duration = sumw.BOT,
           lambda = estimate[1],
           conf.pct = conf.pct,
           mono = TRUE,
           main = main,
           ylim = ylim,
           ...)

    if (history) {

      ## compute the predicted number on each block
      ## Should be at least r.BH
      N.pred <- est.N[1]*w.BH
      r.BH <- as.numeric(r.BH)
      
      ## Il doit y a ooir au moins r.BH données prédites!
      N.pred[N.pred < r.BH] <- r.BH
      
      ## "tapply" would be more efficient less readable here
      for (ib in 1:nlevels(block.H)) {

        ## z.H[as.integer(block.H) == ib] conains the wantedr
        ## odinates, but maybe  not in the right order.
        z.H.ib <- z.H[as.integer(block.H) == ib]
        ind.ib <- rev(order(z.H.ib))
        a <- ( N.pred[ib] + 1 ) / N.pred[ib]
        ww <- a * w.BH[ib] / (1:r.BH[ib]) 
          
        points(x = log(ww),
               y = z.H.ib[ind.ib],
               pch = 24,
               col = "red3",
               bg = NA,
               lwd = 2,
               cex = 1.1)
        
      }
      
    }
    
    ## "U data" x.U are given a return period w.U 
    if (Udata) {
      
      ## points(x = log(w.U), y = x.U,
      ##        pch = 21, col = "purple",
      ##        bg = NA, cex = 1.2)
      
      abline(h = x.U, col = "purple")
      
    }
     
  }
  
  return(res)
  
}


