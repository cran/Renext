##===================================================================
## Author : Y. Deville
##
## For an object of class "Renouv", a Return Level plot is built.
##
##===================================================================

plot.Renouv <- function(x,
                        pct.conf = NULL,
                        show.MAX = TRUE,
                        show.OTS = TRUE,
                        mono = TRUE,
                        rl.mark = NULL,
                        labels.mark = rl.mark,
                        col.mark = NULL,
                        main = NULL,
                        problim = NULL,
                        Tlim = NULL,
                        xlab = "periods",
                        ylab = "level",
                        ...) {

  mc <- match.call( , expand.dots = TRUE)
  ##print(names(mc))
  
  if (mono) {
    l.cols <- c("black", "black")
    p.cols <- "black"
    l.typs <- c("solid", "dashed", "dotted")
  } else {
    l.cols <- c("SteelBlue4", "orangered", "SpringGreen3", "purple", "firebrick")
    p.cols <- "black"
    l.typs <- c("solid", "solid", "solid")
  }

  ## This will NOT be given as ylim but used as such
  ## because, the user is allowed to pass ylim =
  
  yLim <- rangeLev.Renouv(x, Tlim = Tlim)
  ## cat("yLim = ", yLim, "\n")
  
  duration <- x$effDuration
  data <- x$ret.lev
  cnames <- colnames(data)
  
  if (is.null(pct.conf))  pct.conf <- x$pct.conf
      
  l.cols <- rep(l.cols, length.out = length(pct.conf) + 1)
  l.typs <- rep(l.typs, length.out = length(pct.conf) + 1)

  candLnames <- match(paste("L", pct.conf, sep = "."), cnames)
  candUnames <- match(paste("U", pct.conf, sep = "."), cnames)

  threshold <- x$threshold
    
  x.OT <- sort(x$x.OT)
  nx.OT <- length(x$x.OT)

  lambda <- x$estimate["lambda"]
  freq.g <- lambda*(1 - data$prob)
  x.g <- -log(freq.g)

  labs <- c(1, 2,  5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)

  if (FALSE) {
    ry <- c(min(data$L.95, na.rm = TRUE),
            max(data$U.95, na.rm = FALSE))
    
    if (!is.null(x.OT)) {
      ry2 <- range(x.OT, na.rm = TRUE)
    if (ry2[1] < ry[1]) ry[1] <- ry2[1]
      if (ry2[2] > ry[2]) ry[2] <- ry2[2]
    }
  }

  ## if (is.null(ylim)) ylim <- ry

  if (is.null(main)) main <- ""

  if (!is.null(problim)) {
    if ( !is.numeric(problim) || length(problim) != 2 ||
        any(is.na(problim)) ||
        any(problim<= 0) || any(problim>= 1) ) stop("invalid limits in 'problim'.")
    if (!is.null(Tlim)) stop("only one of 'problim' and 'Tlim' can be provided")
    xLim <-  -log(lambda*c(1-problim))
  } else {
    if (!is.null(Tlim)) {
      if ( !is.numeric(Tlim) || length(Tlim) != 2 ||
          any(is.na(Tlim)) ||
          any(Tlim < 0) ) stop("invalid limits in 'Tlim'.")
      if (Tlim[1] < 0.01) Tlim[1] <- 0.01
      xLim <-  log(Tlim)
    } else {
      xLim <- c(min(x.g), log(max(x$pred$period)))  
    }
  }

  ## prepare plot
  plot(x = xLim,
       y = yLim,
       type = "n",
       main = main,
       ## ylim = ylim,
       xlab = xlab,
       ylab = ylab,
       xaxt = "n",
       ##pch = 21, col = l.cols[1],
       ...)

  lines(x = x.g,
        y = data$quant,
        type = "l", lwd = 2,
        col = l.cols[1])
  
  ind <- !is.na(candLnames) & !is.na(candUnames)

  for (i in 1:length(pct.conf))  {
    
    if (ind[i]) {
      iL <- candLnames[i]
      iU <- candUnames[i]
      lines(x = x.g, y = data[, iL],
            type = "l", lwd = 2, col = l.cols[1+i], lty = l.typs[1+i])      
      lines(x = x.g, y = data[ ,iU],
            type = "l", lwd = 2, col = l.cols[1+i], lty = l.typs[1+i])
    } else warning("confidence limits for level ", pct.conf[i], "% not found in data")
    
  }
  
  if (!is.null(x.OT)) {
    f.emp <- (1 - (1:nx.OT)/(nx.OT + 1)) * nx.OT / x$effDuration
    points(x = -log(f.emp),
           y = x.OT,
           pch = 16, col = p.cols[1])

  }

  ## add x-axis
  axis(side = 1, at = log(labs), labels = labs)

  ## add upper x axis ??
  ## axis(side = 3, at = x.g, labels = data$prob)

  abline(v = log(labs), col = "gray", lty = "dotted")
  ## abline(h = pretty(ry), col = "gray", lty = "dotted")
  abline(h = pretty(par()$usr[3:4]), col = "gray", lty = "dotted")
  
  Cols <- c("SteelBlue3", "orangered")

  ## Legend position changed in versions > 1.2
  if (sum(ind)) {
    legend("topleft",
           ## x = x.g[1], y = ry[2] -(ry[2] - ry[2])/10,
           legend = c("theo.", paste(pct.conf[ind], "%", sep = "")),
           col = l.cols[1:(1+sum(ind))],
           lwd = 2,
           lty = l.typs[1:(1+sum(ind))])
  } else {
    legend("topleft",
           ## x = x.g[1], y = ry[2] -(ry[2] - ry[2])/10,
           legend = "theo.",
           col = l.cols[1],
           lwd = 2,
           lty = l.typs[1])
  }

  est.y <- x$estimate[-1]
  
  if (show.MAX && x$history.MAX$flag) {
    
    N.pred <- lambda * x$history.MAX$effDuration 
    w.MAX <- unlist(x$history.MAX$effDuration)
    z.MAX <- unlist(x$history.MAX$data)
    r.MAX <- as.numeric(x$history.MAX$r)
    block.MAX <- x$history.MAX$block

    ## Change by Yves 2012-01-21 old code next line
    ## N.pred[N.pred < r.MAX] <- r.MAX
    indc <- (N.pred < r.MAX)
    if (any(indc)) N.pred[indc] <- r.MAX[indc]
    

    ## identify blocks by background
    Bg <- c("yellow", "SpringGreen1", "mistyrose")
    Bg <- rep(Bg, length.out = length(block.MAX))
    
    ## "tapply" would be more efficient less readable here
    for (ib in 1:nlevels(block.MAX)) {
      
      ## z.MAX[as.integer(block.MAX) == ib] contains the wanted
      ## ordinates, but maybe  not in the right order.
      z.MAX.ib <- z.MAX[as.integer(block.MAX) == ib]
      ind.ib <- rev(order(z.MAX.ib))
      a <- ( N.pred[ib] + 1 ) / N.pred[ib]
      ww <- a * w.MAX[ib] / (1:r.MAX[ib]) 
      
      points(x = log(ww),
             y = z.MAX.ib[ind.ib],
             pch = 24, col = "red3", bg = Bg[ib], lwd = 2, cex = 1.1)
      
    }
  } 
  
  if (show.OTS && x$history.OTS$flag) {

    N.pred <- lambda * x$history.OTS$effDuration *
      (1 - x$funs$F.y(x = x$history.OTS$threshold, parm = est.y))
    w.OTS <- unlist(x$history.OTS$effDuration)
    z.OTS <- unlist(x$history.OTS$data)
    r.OTS <- as.numeric(x$history.OTS$r)
    block.OTS <- x$history.OTS$block
    threshold.OTS <- x$history.OTS$threshold
    
    ## Change by Yves 2012-01-21 old code next line
    ## N.pred[N.pred < r.OTS] <- r.OTS
    indc <- (N.pred < r.OTS)
    if (any(indc)) N.pred[indc] <- r.OTS[indc]
    
    ## identify blocks by background
    Bg <- c("yellow", "SpringGreen1", "mistyrose")
    Bg <- rep(Bg, length.out = length(block.OTS))

    for (ib in 1:nlevels(block.OTS)) {
      
      ## cat("XXX block OTS", ib, "r.OTS[ib] = ", r.OTS[ib], "\n")
      if (r.OTS[ib]>0) {
        
        z.OTS.ib <- z.OTS[as.integer(block.OTS) == ib]
        ind.ib <- rev(order(z.OTS.ib))
        a <- ( N.pred[ib] + 1 ) / N.pred[ib]
        ww <- a * w.OTS[ib] / (1:r.OTS[ib]) 
        
        points(x = log(ww), y = z.OTS.ib[ind.ib],
               pch = 22, col = "SpringGreen1", bg = Bg[ib], lwd = 2, cex = 1.1)
        
      } else {
        ## draw an horizontal segment.
        segments(x0 = par()$usr[1], x1 = log(w.OTS[ib]), 
                 y0 = threshold.OTS[ib], y1 = threshold.OTS[ib],
                 col = "SpringGreen1", lwd = 2)
        points(x = log(w.OTS[ib]), y = threshold.OTS[ib],
               col = "SpringGreen1", pch = 21)
      }
    }
    
  }

  if (!is.null(rl.mark)) {
    if ( mono || (is.null(col.mark)) ) col.mark <- rep("black", length(rl.mark))
    else  col.mark <- rep(col.mark, length.out = length(rl.mark))
    
    for (i in 1:length(rl.mark)) {
      abline(v = log(rl.mark[i]), col = col.mark[i])
      text(x = log(rl.mark[i]),
           y = ry[1] + (ry[2]-ry[1])/6,
           labels = labels.mark[i],
           col = col.mark[i],
           pos = 4)
    }
  }
  
}

