##=============================================================
## 'Prediction' and Inference on return levels using the
## so called 'delta method'
##
## 'newdata' contain return periods.
##
##=============================================================

predict.Renouv <- function(object,
                           newdata = c(10, 20, 50, 100, 200, 500, 1000),
                           cov.rate = FALSE,
                           level = c(0.95, 0.70),
                           trace = 1,
                           eps = 1e-6,
                           ...) {
  
  if ( names(object$estimate)[1L] != "lambda") {
    stop("element 1 of 'estimate' must have name \"lambda\"")
  }
  
  pct.conf <- level*100
  pred.period <- newdata
  
  nc <- 3L + 2L*length(pct.conf)
  
  ## Prepare a 'pred' matrix with dims and names
  
  cnames <- c("period", "prob", "quant",
      paste(rep(c("L", "U"), length(pct.conf)),
            rep(pct.conf, each = 2L), sep = "."))
  
  pred.prob <- 1.0 - 1/object$estimate[1L]/pred.period
  ind <- (pred.prob > 0.0) & (pred.prob < 1.0)
  pred.period <- pred.period[ind]
  pred.prob <- pred.prob[ind]
  
  pred <- matrix(NA, nrow = length(pred.period), ncol = nc)
  rownames(pred) <- pred.period
  colnames(pred) <- cnames
  
  pred[ , "period"] <- pred.period
  pred[ , "prob"]  <- pred.prob
  pred[ , "quant"] <- object$threshold +
    object$funs$q.y(parm = object$estimate[-1], p = pred.prob)

  
  if (object$transFlag) {
    threshold.trans <- object$funs$transfun(object$threshold)
    dth <- threshold.trans - object$threshold
    pred[ , "quant"] <- object$funs$invtransfun(dth +  pred[ , "quant"])
  }

  p.y <- object$p.y
  est.y <- object$estimate[-1]
  fixed.y <- object$fixed.y
  cov.y <- object$cov[-1, -1, drop = FALSE]

  nb.OT <- object$nb.OT
  if ( is.null(nb.OT) ) nb.OT <- length(object$y.OT)
  
  ## cat("nb.OT = ", nb.OT, "\n")
  ## nb.OT <- object$nb.OT
  
  if ( (object$distname.y == "exponential") && !object$history.MAX$flag
      && !object$history.OTS$flag ) {
    
    if (trace) cat("Special inference for the exponential case without history\n")
    
    ##--------------------------------------------------------------
    ## Use the sampling distribution to derive confidence
    ## limits on quantiles. If 'lambda' was known, these limits
    ## would be exact.
    ##--------------------------------------------------------------
    
    for (ipct in 1:length(pct.conf)) {
      
      alpha.conf <- (100 - pct.conf[ipct])/100
      
      theta.L <- 2*nb.OT / est.y / qchisq(1 - alpha.conf/2, df = 2 * nb.OT)
      theta.U <- 2*nb.OT / est.y / qchisq(alpha.conf/2, df = 2 * nb.OT)
      
      pred[ , 2L * ipct + 2L] <-
        object$threshold + qexp(p = pred.prob, rate = 1.0/theta.L)
      pred[ , 2L * ipct + 3L] <-
        object$threshold + qexp(p = pred.prob, rate = 1.0/theta.U) 
      
    }
    
    if (object$transFlag) {

      ## dth was defined above
      
      for (ipct in 1L:length(pct.conf)) { 
        pred[ , 2L * ipct + 2L] <-
          object$funs$invtransfun(dth + pred[ , 2L * ipct + 2L])
        pred[ , 2L * ipct + 3L] <-
          object$funs$invtransfun(dth + pred[ , 2L * ipct + 3L]) 
      }
      
    }
    
    
    ##==============================================================
    ## perform Bartlett's gof test
    ##==============================================================

    infer.method <-
      "chi-square for exponential distribution (no historical data)"
    
  } else {
     
    ##==============================================================
    ## DELTA METHOD
    ## matrices for numerical derivation
    ## The column ip of Parmat contains the parameter value with
    ## a tiny modification of its ip component.
    ##==============================================================

    sig <- rep(0, length(pred.prob))

    if (p.y > 0L) { ## ADDED 2013-09-07
      
      parmMat <- matrix(est.y[!fixed.y], nrow = p.y, ncol = p.y)
      
      rownames(parmMat) <- object$parnames.y[!fixed.y]
      colnames(parmMat) <- object$parnames.y[!fixed.y]
      
      dparms <- abs(est.y[!fixed.y])*eps
      dparms[dparms < eps] <- eps
      
      for (ip in 1:p.y) parmMat[ip, ip] <- parmMat[ip, ip] + dparms[ip]
    
      delta <- rep(NA, p.y)
      ## sig <- rep(NA, length(pred.prob)) ## MOVED 
    
      ##=============================================================
      ## Compute the delta's  and sig's
      ## delta contains derivative w.r.t. unknown params
      ##
      ## CAUTION
      ##
      ## The quantile function 'q.y' takes a vector as first arg which
      ## should morally be of length 'p.y', but is here of length
      ## parnb.y. This works because 'q.y' uses the elements in named
      ## form, irrespective of their position.
      ##
      ##=============================================================
      
      for (i in 1:length(pred.prob)) {
        
        for (ip in 1:p.y) {
          
          est.prov <- est.y
          est.prov[!fixed.y] <-  parmMat[ , ip]
          
          dd <- ( object$funs$q.y(est.prov, p = pred.prob[i]) -
                 object$funs$q.y(est.y, p = pred.prob[i]) ) / dparms[ip]
          delta[ip] <- dd
          
        }
        
        sig[i] <- sqrt(t(delta)%*%cov.y[!fixed.y, !fixed.y]%*%delta)
        
      }
    }  ## ADDED 2013-09-07
      
    
    ##=============================================================
    ## Compute the delta's  and sig's
    ## delta contains derivative w.r.t. unknown params
    ##=============================================================

    pred.sig <- rep(0, length(pred.prob))

    if (p.y > 0L){
      for (i in 1:length(pred.prob)) {
        for (ip in 1:p.y) {
          est.prov <- est.y
          est.prov[!fixed.y] <-  parmMat[ , ip]
          dd <- ( object$funs$q.y(est.prov, p = pred.prob[i]) -
                 object$funs$q.y(est.y, p = pred.prob[i]) ) / dparms[ip]
          delta[ip] <- dd
        }
        ## modif 2010-01-25 for the fixed parms case
        ##pred.sig[i] <- sqrt(t(delta)%*%cov.y%*%delta)
        pred.sig[i] <- sqrt(t(delta) %*% cov.y[!fixed.y, !fixed.y]%*%delta)
      }
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
      pred[ , 2*ipct + 2] <-
        object$threshold + object$funs$q.y(est.y, p = pred.prob) - z.conf * pred.sig
      pred[ , 2*ipct + 3] <-
        object$threshold + object$funs$q.y(est.y, p = pred.prob) + z.conf * pred.sig
    }

    if (object$transFlag) {
      for (ipct in 1L:length(pct.conf)) {
        pred[ , 2L*ipct + 2L] <- object$funs$invtransfun(dth  + pred[ , 2L*ipct + 2L])
        pred[ , 2L*ipct + 3L] <- object$funs$invtransfun(dth  + pred[ , 2L*ipct + 3L])
      }
    }
    
    infer.method <- "Delta method"
  }

  pred <- as.data.frame(pred)
  attr(pred, "infer.method") <- infer.method
  pred

}
