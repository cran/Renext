##*****************************************************************************
## AUTHOR: Yves Deville
##
## CONTAINS
##
## o The 'plot' and 'lines' methods for the "Renouv" class.
##
## o Misc utility functions to build legends on Return levels plot
##   in a semi-automatic fashion.
##
##*****************************************************************************

`RLpar0` <- function(mono = TRUE) {

  ##=======================================================================
  ## Set default graphical parameters for Return Level plots
  ##========================================================================

  if (mono) {
    l.cols <- c("black", "black")
    p.cols <- "black"
    l.typs <- c("solid", "dashed", "dotted")
  } else {
    l.cols <- c("SteelBlue4", "orangered", "SpringGreen3", "purple", "firebrick")
    p.cols <- "black"
    l.typs <- c("solid", "solid", "solid")
  }
  block.cols <- c("yellow", "SpringGreen1", "mistyrose", "cyan", "SteelBlue1",
                  "gold")
 
  l.cols <- rep(l.cols, length.out = 10L)
  p.cols <- rep(p.cols, length.out = 10L)
  l.typs <- rep(l.typs, length.out = 10L)
  block.cols <- rep(block.cols, length.out = 10L)
  
  .RLpar <- list()
  
  ## quantile (return level) curve
  .RLpar[["quant"]] <- list(type = "l", col = "black", lwd = 2, lty = "solid")
  
  ## sample data
  .RLpar[["OT"]] <- list(col = "black", pch = 16, cex = 0.8, bg = "black")
  
  ## confidence levels
  prov <- list()
  nconf <- 6L
  for (i in 1L:nconf) {
    prov[[paste("conf", i, sep = "")]] <-
      list(lty = i+1L, col = l.cols[i], lwd = 2)
  }

  .RLpar[["conf"]] <- prov

  ## historical MAX data
  prov <- list()
  nMAX <- 10L
  for (i in 1L:nMAX) {
    prov[[paste("block", i, sep = "")]] <-
      list(col = "black", pch = 24, cex = 1.1, lwd = 2,  bg = block.cols[i])
  }
  .RLpar[["MAX"]] <- prov

  ## OTS data 
  prov <- list()
  nOTS <- 10L
  for (i in 1L:nOTS) {
    prov[[paste("block", i, sep = "")]] <-
      list(col = "black", pch = 22, cex = 1.1, lwd = 2, bg = block.cols[i])
  }
  .RLpar[["OTS"]] <- prov

  ## Modify the default values ???
  
  .RLpar
  
}

##*****************************************************************************

`RLpar` <- function(mono = TRUE, ...) {
  
  ##=======================================================================
  ## Set the graphical parameters for Return Level plots OR CHANGE THEM
  ##========================================================================
  
  mc <- match.call(expand.dots = TRUE)
  
  newPar <- as.list(mc)
  newPar[[1]] <- NULL
  newPar[["mono"]] <- NULL
  
  oldPar <- RLpar0(mono = mono)
  
  dp <- duplicated(names(newPar))

  if (any(dp)) {
    warning("dupplicated par names: ",
            paste("'", names(newPar)[dp], "'", sep = "", collapse = ", "))
  }
  
  uOldPar <- unlist(oldPar)
  nm0 <- match(names(newPar), names(uOldPar))
  if (any(is.na(nm0))) {
    warning("'from' contains unused names:\n",
            sprintf("\"%s\"", names(newPar)[is.na(nm0)]))
  }
  
  ## Replacment of the values
  Names0 <- names(newPar)[!is.na(nm0)]
  
  ## This can not be vectorized since nmVec is a vector with elements for
  ## the hierarchical levels 1, 2, ...
  for (Name0 in Names0){
    nmVec <- unlist(strsplit(Name0, split = "\\."))
    oldPar[[nmVec]] <- newPar[[Name0]]
  }
  
  oldPar
 
}

##*****************************************************************************

`rReplace` <- function(from, to, trace = 0) {

  ##===========================================================================
  ## Replace values in the list 'to' using the values
  ## from the list 'from' in elements having the same
  ## (hierarchical) name 
  ##===========================================================================
  
  ulFrom <- unlist(from)
  ulTo <- unlist(to)
  nm0 <- match(names(ulFrom), names(ulTo))
    
  if (any(is.na(nm0))) {
    warning("'from' contains unused names:\n",
            sprintf("\"%s\"", names(ulFrom)[is.na(nm0)]))
  }
  
  Names0 <- names(ulFrom)[!is.na(nm0)]

  ## This can not be vectorized since nmVec is a vector with elements for
  ## the hierarchical levels 1, 2, ...
  
  for (Name0 in Names0){
    nmVec <- unlist(strsplit(Name0, split = "\\."))
    to[[nmVec]] <- from[[nmVec]]
  }
  
  to 
  
}

##*****************************************************************************
`RLlegend.ini` <- function(x = "topleft", bty = "n", envir, ...) {
  .RLlegend <- list(x = x, bty = bty, ...)
  assign(".RLlegend", .RLlegend, envir = envir)
  .RLlegend
}

##*****************************************************************************
`RLlegend.show` <- function(envir) {
  .RLlegend <- get(".RLlegend", envir = envir, mode = "list")
  do.call("legend", .RLlegend)
}

##*****************************************************************************

`plot.Renouv` <-
  function(x,
           pct.conf = x$pct.conf,
           show = list(OT = TRUE, quant = TRUE, conf = TRUE, MAX = TRUE,
             OTS = TRUE),
           mono = TRUE,
           predict = FALSE, 
           par = NULL,
           legend = TRUE,
           legendEnvir = NULL,
           label = NULL,
           problim = NULL,
           Tlim = NULL,
           main = NULL,
           xlab = "periods",
           ylab = "level",
           ...) {
    
    ##==========================================================================
    ## The plot method
    ##==========================================================================
    
    mc <- match.call( , expand.dots = TRUE)
    RLpar <- par 
    
    ## This will NOT be given as 'ylim' but used as such
    ## because the user is allowed to pass ylim =
    lambda <- x$estimate["lambda"]
    
    yLim <- rangeLev.Renouv(x, Tlim = Tlim)
    
    labs <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
    
    if (is.null(main)) main <- ""

    ## compute suitable xlim 
    if (!is.null(problim)) {
      
      if ( !is.numeric(problim) || length(problim) != 2L ||
          any(is.na(problim)) || any(problim <= 0) || any(problim >= 1) ) {
        stop("invalid limits in 'problim'.")
      }
      
      if (!is.null(Tlim)) {
        stop("only one of 'problim' and 'Tlim' arguments can be provided")
      }

      xLim <-  -log(lambda*c(1-problim))

    } else {
      
      if (!is.null(Tlim)) {
        
        if ( !is.numeric(Tlim) || length(Tlim) != 2 || any(is.na(Tlim)) ||
            any(Tlim < 0) ) {
          stop("invalid limits in 'Tlim'.")
        }
        
        if (Tlim[1] < 0.01) Tlim[1] <- 0.01
        xLim <-  log(Tlim)

      } else {
        xMin <- max(c(0,  log(min(x$ret.lev$period))), na.rm = TRUE)
        xLim <- c(xMin, log(max(x$pred$period)))  
      }
      
    }
    
    ## prepare (empty) plot
    plot(x = xLim,
         y = yLim,
         type = "n",
         main = main,
         xlab = xlab,
         ylab = ylab,
         xaxt = "n",
         ...)

    ## add x-axis and grid
    axis(side = 1, at = log(labs), labels = labs)
    abline(v = log(labs), col = "gray", lty = "dotted")
    abline(h = pretty(par()$usr[3:4]), col = "gray", lty = "dotted")

    ## prepare legend
    if (legend) {
      if (is.null(legendEnvir)) legendEnvir <- new.env()
      RLlegend.ini(envir = legendEnvir)
      if (is.null(label)) {
        label <- deparse(substitute(x))
      }
    }
    
    ## Call 'lines' method
    lines.Renouv(x = x,
                 pct.conf = x$pct.conf,
                 show = show,
                 mono = mono,
                 predict = predict,
                 RLpar = RLpar,
                 legend = FALSE,
                 legendEnvir = legendEnvir,
                 label = label,
                 ...)

    ## draw legend
    if (legend) {
      RLlegend.show(envir = legendEnvir)
    }
    
  }

##*****************************************************************************

`lines.Renouv` <-
  function(x,
           pct.conf = x$pct.conf,
           show = NULL,
           mono = TRUE,
           predict = FALSE,
           par = NULL,
           legend = FALSE,
           legendEnvir = NULL,
           label = NULL,
           ...) {

    DEBUG <- FALSE
    RLpar <- par
   
    ##*************************************************************************
    ## lines method for objects with class "Renouv" 
    ##*************************************************************************
    nx.OT <- length(x$x.OT)   
    ## opar <- par(no.readonly = TRUE)
    
    lambda <- x$estimate["lambda"]
    
    ##=========================================================================
    ## shown objects
    ##========================================================================
    
    .show <- list(quant = TRUE, conf = FALSE,
                  OT = FALSE, MAX = FALSE, OTS = FALSE)

    if (!is.null(show)) {
      if (!is.list(show)) stop("'show' must be a list")
      .show <- rReplace(from = show, to = .show)
    }
    
    .par <- RLpar(mono = mono)
    
    if (!is.null(RLpar)) {
      if (!is.list(RLpar)) stop("'par' must be a list")
      .par <- rReplace(from = RLpar, to = .par)
    }

    ##=========================================================================
    ## analyse 'label' and prepare legend
    ##=========================================================================
    
    if (is.null(label)) {
      label <- deparse(substitute(x))
    }
    ## 
    if (is(label, "character")) {
      Label <- list(quant = paste(label, "quant"))
      if ( .show[["conf"]] && length(pct.conf) ) {  ## pct
        Label$conf <- list()
        for (il in 1L:length(pct.conf)) {
          ## cat("il = ", il, "\n")
          Label$conf[[paste("conf", il, sep = "")]] <-
            paste(label, " ", pct.conf[il], "%", sep = "")
          ## cat(paste(label, pct.conf[il]), "\n")
        }
      }
      
      if (.show[["OT"]] && (nx.OT > 0L) ) {
        Label$OT <- paste(label, "sample")
      }

      if (.show[["MAX"]] && x$history.MAX$flag ) {
        
        Label$MAX <- list()
        
        if (nlevels(x$history.MAX$block) > 1L) { 

          if (nlevels(x$history.MAX$block) > 10L) {
            stop("too many (> 10) 'MAX' blocks for plotting.",
                 " Use show(MAX = FALSE).")
          }
            
          for (ib in 1L:10L) {
            Label$MAX[[paste("block", ib, sep = "")]] <-
              paste(label, " MAX block", ib,  sep = "")
          }
          
        } else {
          Label$MAX[["block1"]] <- "hist. MAX"
        }

      }
      
      if (.show[["OTS"]] && x$history.OTS$flag ) {
     
        Label$OTS <- list()
        
        if (nlevels(x$history.OTS$block) > 1L) {

          if (nlevels(x$history.MAX$block) > 10L) {
            stop("Too many (> 10) 'OTS' blocks for plotting. ",
                 "Use 'show(OTS = FALSE)'.")
          }
          
          for (ib in 1L:10L) {
            Label$OTS[[paste("block", ib, sep = "")]] <-
              paste(label, " OTS block", ib,  sep = "")
          }
        } else {
          Label$OTS[["block1"]] <- "hist. OTS"
        }
        
      }
      
    } else if (is(label, "list")) {
      Label <- label
    }

    ##=========================================================================
    ## predict if necessary
    ##
    ##=========================================================================
    
    ## limits in years with min >= 1 year
   
    xLim <- par()$usr[c(1L, 2L)]
    if (xLim[1L] < 0)  xLim[1L] <- 0
    
    ## if (length(pct.conf) > 0) pct.conf <- rev(sort(pct.conf))
    
    if (predict) {
      x.g <- seq(from = xLim[1L], to = xLim[2L], length.out = 100L)
      if (length(pct.conf) > 0) lev.conf0 <- pct.conf / 100
      else lev.conf0 <- 0.95
      Data <- predict.Renouv(object = x, newdata = exp(x.g), level = lev.conf0)
    } else {
      Data <- x$ret.lev
    }
    
    x.g <- log(Data$period) ## x.g is replaced

    if (!is.null(legendEnvir)) {
      .RLlegend <- get(".RLlegend", legendEnvir)
    } else {
      .RLlegend <- list() 
    }
    
    ##=========================================================================
    ## plot the quantile curve
    ##=========================================================================
    
    if ( .show[["quant"]]) {

      parList <- c(list(x = x.g, y = Data$quant), .par[["quant"]])
      do.call(lines, parList)
      
      ## XXX
      ## lines(x = x.g, y = Data$quant, col = "green", type = "o", lty = 2)

      ##-----------------------------------------------------------------------
      ## update legend
      ##-----------------------------------------------------------------------
      
      if (is.null(Label[["quant"]])) {
        stop("'label' must be a character or a list with a 'quant' element")
      }
      
      .RLlegend$legend <- c(.RLlegend$legend, Label[["quant"]])
      
      lty.prov <- .par[["quant"]][["lty"]]
      par(lty = lty.prov)
      lty.prov <- par()$lty
      .RLlegend[["lty"]] <- c(.RLlegend[["lty"]], lty.prov)
      .RLlegend[["pch"]] <- c(.RLlegend[["pch"]], NA)
      .RLlegend[["pt.bg"]] <- c(.RLlegend[["pt.bg"]], NA)
      
      for (nm in c("col", "lwd")) {
        .RLlegend[[nm]] <- c(.RLlegend[[nm]], .par[["quant"]][[nm]])
      }
    }

    ##=========================================================================
    ## plot confidence limits
    ##=========================================================================
    
    if ( .show[["conf"]] && length(pct.conf) ) {
      
      if (DEBUG) cat("showing 'conf' lines ... \n")

      cnames <- colnames(Data)
      candLnames <- match(paste("L", pct.conf, sep = "."), cnames)
      candUnames <- match(paste("U", pct.conf, sep = "."), cnames)
      ind <- !is.na(candLnames) & !is.na(candUnames)
      
      for (i in 1L:length(pct.conf))  {
  
        if (ind[i]) {
          iL <- candLnames[i]
          iU <- candUnames[i]
          parList <- c(list(x = x.g, y = Data[, iL]), .par[["conf"]][[i]])
          do.call(lines, parList)
          parList$y <- Data[ , iU]
          do.call(lines, parList)

          ## update legend
          nm <- paste("conf", i, sep = "")
      
          if (is.null(Label$conf[[nm]])) {
            stop("'label' must be a character or a list with a '",
                 nm, "' element")
          }
          
          .RLlegend$legend <- c(.RLlegend$legend, Label$conf[[nm]])

          ## translate 'lty' if needed using 'par' function
          lty.prov <- .par[["conf"]][[nm]][["lty"]]
          par(lty = lty.prov)
          lty.prov <- par( )$lty
          .RLlegend[["lty"]] <- c(.RLlegend[["lty"]], lty.prov)
          .RLlegend[["pch"]] <- c(.RLlegend[["pch"]], NA)
          .RLlegend[["pt.bg"]] <- c(.RLlegend[["pt.bg"]], NA)
          
          for (nm2 in c("col", "lwd")) {
            .RLlegend[[nm2]] <- c(.RLlegend[[nm2]], .par[["conf"]][[nm]][[nm2]])
          }

          if (DEBUG) {
            cat("YYY\n")
            print(.RLlegend)
          }
          
        } else {
          warning("confidence limits for level ",
                  pct.conf[i], "% not found in data")   
        }
      }
   
    }
    
    ##=========================================================================
    ## show sample points if wanted, and if there are some
    ##=========================================================================
    
    if (.show[["OT"]] && (nx.OT > 0L) ) {
      
      if (DEBUG) cat("showing 'OT' data ... \n")

      x.OT <- sort(x$x.OT)
      f.emp <- (1 - (1:nx.OT)/(nx.OT + 1)) * nx.OT / x$effDuration
      
      parList <- c(list(x = -log(f.emp), y = x.OT), .par[["OT"]])
      do.call(points, parList)
      
      .RLlegend$legend <- c(.RLlegend$legend, Label$OT)
      .RLlegend[["lty"]] <- c(.RLlegend[["lty"]], NA)
      .RLlegend[["col"]] <- c(.RLlegend[["col"]], .par[["OT"]][["col"]])
      .RLlegend[["lwd"]] <- c(.RLlegend[["lwd"]], NA)
      .RLlegend[["pch"]] <- c(.RLlegend[["pch"]], .par[["OT"]][["pch"]])
      .RLlegend[["pt.bg"]] <- c(.RLlegend[["pt.bg"]], .par[["OT"]][["bg"]])
      
      if (DEBUG) {
        print(.RLlegend)
        cat("... 'OT' data shown\n")
      }
    }

    
    ##=======================================================================
    ## show 'MAX' historical data if wanted and if there are some
    ##=======================================================================
    
    if (.show[["MAX"]] && x$history.MAX$flag ) {
      
      if (DEBUG) {
        cat("showing 'MAX' data ... \n",
            "==================     \n")
      }
       
      N.pred <- lambda * x$history.MAX$effDuration 
      w.MAX <- unlist(x$history.MAX$effDuration)
      z.MAX <- unlist(x$history.MAX$data)
      r.MAX <- as.numeric(x$history.MAX$r)
      block.MAX <- x$history.MAX$block
      
      indc <- (N.pred < r.MAX)
      if (any(indc)) N.pred[indc] <- r.MAX[indc]
      
      ## "tapply" would be more efficient less readable here
      for (ib in 1:nlevels(block.MAX)) {
        
        z.MAX.ib <- z.MAX[as.integer(block.MAX) == ib]
        ind.ib <- rev(order(z.MAX.ib))
        a <- ( N.pred[ib] + 1 ) / N.pred[ib]
        ww <- a * w.MAX[ib] / (1:r.MAX[ib])
        parList <- c(list(x = log(ww), y = z.MAX.ib[ind.ib]), .par[["MAX"]][[ib]])
        do.call(points, parList)
        
        ## update legend
        nm <- paste("block", ib, sep = "")
        
        if (is.null(Label$MAX[[nm]])) {
          stop("'label' must be a character or a list with a \"MAX\"",
               " element containing a '", nm, "' element")
        }

        .RLlegend$legend <- c(.RLlegend$legend, Label$MAX[[nm]])
        
        ## points properties.
        .RLlegend[["lty"]] <- c(.RLlegend[["lty"]], NA)
        .RLlegend[["lwd"]] <- c(.RLlegend[["lwd"]], NA)
        
        ## points properties 'lwd' is for empty symbols pch = 21 to 26 
        .RLlegend[["col"]] <- c(.RLlegend[["col"]], .par[["MAX"]][[nm]][["col"]])
        .RLlegend[["pch"]] <- c(.RLlegend[["pch"]], .par[["MAX"]][[nm]][["pch"]])
        .RLlegend[["pt.bg"]] <- c(.RLlegend[["pt.bg"]], .par[["MAX"]][[nm]][["bg"]]) 
        .RLlegend[["pt.lwd"]] <- c(.RLlegend[["pt.lwd"]], .par[["MAX"]][[nm]][["lwd"]])
    
      }
      
      if (DEBUG) {
        print(.RLlegend)
        cat("... 'MAX' data shown\n")
      }
     
    }
    
    ##=======================================================================
    ## show 'OTS' historical data if wanted and if there are some
    ##=======================================================================
    
    if (.show[["OTS"]] && x$history.OTS$flag) {

      if (DEBUG) {
        cat("showing 'OTS' data ... \n",
            "==================     \n")
      }
      
      est.y <- x$estimate[-1]
      
      N.pred <- lambda * x$history.OTS$effDuration
      w.OTS <- unlist(x$history.OTS$effDuration)  ## is 'unlist' really useful here?
      z.OTS <- unlist(x$history.OTS$data) 
      r.OTS <- as.numeric(x$history.OTS$r)
      block.OTS <- x$history.OTS$block
      threshold.OTS <- x$history.OTS$threshold

      ## 'Npred' can be less than the observed number!
      indc <- (N.pred < r.OTS)
      if (any(indc)) N.pred[indc] <- r.OTS[indc]
      
      for (ib in 1:nlevels(block.OTS)) {
        
        ## cat("XXX block OTS", ib, "r.OTS[ib] = ", r.OTS[ib], "\n")
        if (r.OTS[ib]>0) {
          
          z.OTS.ib <- z.OTS[as.integer(block.OTS) == ib]
          ind.ib <- rev(order(z.OTS.ib))
          a <- ( N.pred[ib] + 1L ) / N.pred[ib]
          ww <- a * w.OTS[ib] / (1L:r.OTS[ib]) 

          parList <- c(list(x = log(ww), y = z.OTS.ib[ind.ib]), .par[["OTS"]][[ib]])
          do.call(points, parList)
          
          ## update legend
          nm <- paste("block", ib, sep = "")
        
          if (is.null(Label$OTS[[nm]])) {
            stop("'label' must be a character or a list with a \"OTS\"",
                 " element containing a '", nm, "' element")
          }

          .RLlegend$legend <- c(.RLlegend$legend, Label$OTS[[nm]])
          
          ## points properties.
          .RLlegend[["lty"]] <- c(.RLlegend[["lty"]], NA)
          .RLlegend[["lwd"]] <- c(.RLlegend[["lwd"]], NA)
          
          ## points properties 'lwd' is for empty symbols pch = 21 to 26 
          .RLlegend[["col"]] <- c(.RLlegend[["col"]], .par[["OTS"]][[nm]][["col"]])
          .RLlegend[["pch"]] <- c(.RLlegend[["pch"]], .par[["OTS"]][[nm]][["pch"]])
          .RLlegend[["pt.bg"]] <- c(.RLlegend[["pt.bg"]], .par[["OTS"]][[nm]][["bg"]]) 
          .RLlegend[["pt.lwd"]] <- c(.RLlegend[["pt.lwd"]], .par[["OTS"]][[nm]][["lwd"]])
          

        } else {
          ## draw an horizontal segment.
          segments(x0 = par()$usr[1], x1 = log(w.OTS[ib]), 
                   y0 = threshold.OTS[ib], y1 = threshold.OTS[ib],
                   col = "SpringGreen1", lwd = 2)
          points(x = log(w.OTS[ib]), y = threshold.OTS[ib],
                 col = "SpringGreen1", pch = 21)
        }
      }

      if (DEBUG) {
        print(.RLlegend)
        cat("... 'OTS' data shown\n")
      }
      
    }

    if (!is.null(legendEnvir)) {
      assign(".RLlegend", .RLlegend, envir = legendEnvir)
    }
      
    ## par(opar)
    
  }


