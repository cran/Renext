##******************************************************************
## Author: Yves Deville <deville.yves@alpestat.com> 
## 
## Fit aggregated POT renewal assuming GPD exccedances.   
## 
##****************************************************************** 

fGEV.MAX <- function(MAX.data,
                     MAX.effDuration,
                     MAX = NULL,
                     control = list(maxit = 300, fnscale = -1),
                     scaleData = TRUE,
                     cov = TRUE,
                     info.observed = TRUE,
                     trace = 0) {

    if (is.null(MAX)) {
        if (missing(MAX.data) || missing(MAX.effDuration)) {
            stop("when 'MAX' is not given,  'MAX.data' and 'MAX.effDuration' ",
                 "must be given")
            
        }
        MAX <- makeMAXdata(x = NULL,
                           data = MAX.data,
                           effDuration = MAX.effDuration) 
    }
    
    fixed.par.y <- NULL
    distname.y <- "GPD"

    ##========================================================================
    ## A very small value of 'delta' can be used. The smaller is
    ## 'delta', the smaller will be the rate 'lambda' which may be
    ## better for numeric precision. However, 'lambda' never gets very
    ## large as explained in 'Computing Details'.
    ## ========================================================================
    
    z.MAX <- unlist(MAX$data)
    ## delta <- mean(spacings(z.MAX)) / length(z.MAX)
    delta <- mean(spacings(z.MAX)) * 1e-6
    threshold <- min(z.MAX) - delta
    if (trace) cat("o Threshold (original scale) ", threshold, "\n\n")
    
    if (scaleData) {
        mexpon.OT <- floor(log(mean(z.MAX - threshold), base = 10))
        if (mexpon.OT >= 2) {
            scale.need <- TRUE
            scale.OT <- 10^mexpon.OT
        } else {
            scale.need = FALSE
            scale.OT <- 1.0
        }
        if (trace) cat("o Using scale factor ", scale.OT, "\n")
        zMod.MAX   <- (z.MAX - threshold) / scale.OT
        if (trace) {
            cat("Range of scaled exceedances\n")
            print(range(zMod.MAX))
            cat("\n")
        }
        
    } else {
        scale.OT <- 1.0
        zMod.MAX  <- z.MAX - threshold
    }
    
    par.ini <- parIni.MAX(MAX, threshold = threshold, distname.y = distname.y)
    par.ini["scale"] <- par.ini["scale"] / scale.OT
    
    myDist <-  checkDist(distname.y = distname.y,
                         scale.OT = scale.OT)
    
    funs <-  makeFuns(funname.y = myDist$funname.y,
                      parnames.y = myDist$parnames.y,
                      fixed.par.y = fixed.par.y,
                      trace = 0)
    
    ## z.MAX <- unlist(MAX$data)
    w.MAX <- MAX$effDuration
    block.MAX <- MAX$block
    r.MAX <- MAX$r
    nblock.MAX <- length(r.MAX)
    zrMod.MAX <- tapply(zMod.MAX, block.MAX, min)
    r.TOT <- sum(r.MAX)
    
    if ( (nblock.MAX > 1L) && (sd(w.MAX) > 0.0) ) {
        stop("GEV fit is at the time only possible for constant duration")
    }
    
    ##=========================================================================
    ## Function needed to compute 'lambda' at optimum 
    ##=========================================================================
    
    lambdaHatFun <- function(parms.y) {
        S.MAX <- (1.0 - funs$F.y(x = zrMod.MAX, parm = parms.y))
        w.TOT <- sum(w.MAX * S.MAX)
        r.TOT / w.TOT 
    }
    
    ##=========================================================================
    ## Concentrated version (with respect to lambda. Caution: lambda is not in
    ## 'parms' here!!!
    ##=========================================================================
    
    loglikC <- function(parms.y) {
        S.MAX <- (1.0 - funs$F.y(x = zrMod.MAX, parm = parms.y))
        w.TOT <- sum(w.MAX * S.MAX)
        lambda.hat <- r.TOT / w.TOT
        logL  <- r.TOT * log(lambda.hat)   
        logL <- logL + sum(funs$logf.y(parm = parms.y, x = zMod.MAX))
        logL
    }
    
    opt <- optim(par = par.ini[-1L],
                 fn = loglikC,
                 control = control)
    
    par.GPD <- c(lambda = lambdaHatFun(opt$par), opt$par)
        
    ## add a constant to match the logLik computed in
    ## the classical way
    logLik <- opt$value - r.TOT * log(scale.OT) - r.TOT
        
    if (trace) {
        cat("o estimated GPD parameters\n")
        print(par.GPD)
        cat("\n")
    }

    
    if (!cov) {
        ## transform to data scale
        par.GPD["scale"] <- par.GPD["scale"] * scale.OT
        par.GEV <- Ren2gev(par.GPD,
                           threshold = threshold,
                           jacobian = FALSE,
                           vcovRen = NULL)
        return(list(estimate = par.GEV[1L:3L],
                    opt = opt,
                    loglik = logLik))
    }
    
    ##=======================================================================
    ## Use exact derivatives to find the covariance matrix
    ## See the document "Renext Computing Details" for the details
    ##=======================================================================

    
    if (par.GPD["shape"] < -0.5) {
        ## transform to data scale
        par.GPD["scale"] <- par.GPD["scale"] * scale.OT
        par.GEV <- Ren2gev(par.GPD,
                           threshold = threshold,
                           jacobian = FALSE,
                           vcovRen = NULL)
        warning("estimated 'shape' is < -0.5. ML inference results not suitable")
        sds <- rep(NA, 3L)
        vcov <- matrix(NA, nrow = 3L, ncol = 3L)
        names(sds) <- rownames(vcov) <- colnames(vcov) <-
            c("loc", "scale", "shape")

    } else {

        if (info.observed) {
            
            ## compute the derivative at the 'zr' values (min in each block)
            DerMin <- parDeriv(par = par.GPD[-1L], x = zrMod.MAX, distname = "gpd",
                               sum = FALSE) 
            der2Smod <- sweep(x = DerMin$der2Surv, MARGIN = 1L,
                              STATS = -par.GPD["lambda"] * w.MAX, FUN = "*")
            der1S <- sweep(x = DerMin$derSurv, MARGIN = 1L,
                           STATS = -w.MAX, FUN = "*")
            vec <- apply(der1S, 2L, sum)
            names(vec) <- names(opt$par)
            mat2 <- apply(der2Smod, 2L:3L, sum)
            
            ## add the contribution of the log-density
            DerAll <- parDeriv(par = par.GPD[-1L], x = zMod.MAX, distname = "gpd",
                               sum = TRUE) 
            mat2 <- mat2 + DerAll$der2Logdens
            scal <- c("lambda" = -sum(r.MAX / par.GPD["lambda"] / par.GPD["lambda"]))
            info <- -rbind(c(scal, vec), cbind(vec, mat2))
            rownames(info) <- colnames(info)
            cov.Ren <- try(solve(info), silent = TRUE)
            
            if (inherits(cov.Ren, "try-error")) {
                warning("'info' could not be inverted")
                cov.Ren <- matrix(NA, nrow = 3L, ncol = 3L)
                rownames(cov.Ren) <- colnames(cov.Ren) <- colnames(info)
            }
            
            ## transform to data scale
            par.GPD["scale"] <- par.GPD["scale"] * scale.OT
            cov.Ren[ , "scale"] <- cov.Ren[ , "scale"] * scale.OT
            cov.Ren["scale", ] <- cov.Ren["scale", ] * scale.OT
          
            par.GEV <- Ren2gev(par.GPD,
                               threshold = threshold,
                               vcovRen = cov.Ren)
            
            vcov <- attr(par.GEV, "vcov")
            sds <- sqrt(diag(vcov))
            
        } else {
            
            if (!all(r.MAX == 1L)) {
                stop("'info.observed = FALSE' is only possible when ",
                     "all MAX blocks contain r = 1 observation")
            }
            
            ## transform to data scale
            par.GPD["scale"] <- par.GPD["scale"] * scale.OT
            par.GEV <- Ren2gev(par.GPD,
                               threshold = threshold,
                               jacobian = FALSE,
                               vcovRen = NULL)
            
            info <- infoGEV(par.GEV, n = r.TOT)
            vcov <- solve(info)
            sds <- sqrt(diag(vcov))
            
        }

    }
    
    list(estimate = par.GEV[1L:3L],
         opt = opt,
         loglik = logLik,
         sd = sds,
         cov = vcov)
    
}
