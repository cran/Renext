## Methode de renouvellement

rRenouv <- function(densfun.y = "exponential",
                    par.y = list(rate = 1),
                    densfun.N = "poisson",
                    par.N = list(lambda = 6),
                    threshold = 0,
                    aggreg = TRUE,
                    nb = 50,
                    labb = seq(to = 2009, by = 1, length = nb),
                    w = rep(1, nb)) {
  
  if (is.character(densfun.N)) {
    dists.N <- c("poisson", "negative binomial")
    
    if (!densfun.N %in% dists.N) stop(paste("si character, densfun.N peut etre:",
                                            paste(dists.N, collapse = ", ")))
  } else {
    if (!is (densfun.N, "function")) stop("densfun.N doit etre character ou function")
  }
  
  if (is.character(densfun.y)) {
    dists.y <- c("exponential", "weibull", "gpd")
    
    if (!densfun.y %in% dists.y) stop(paste("si character, densfun.y peut etre:",
                                            paste(dists.y, collapse = ", ")))
  } else {
    if (!is (densfun.y, "function")) stop("densfun.y doit etre character ou function")
  }
  
  ## simulate the counts
  
  if (aggreg) {
    
    if (densfun.N =="negative binomial") {
      if (is.null(par.N$gamma)) stop("par.N doit contenir un parametre gamma")
      if (is.null(par.N$prob)) stop("par.N doit contenir un parametre prob")
      
      N <- rnbinom(n = nb, size = par.N$gamma*w, prob = par.N$prob)
      
    } else if (densfun.N =="poisson") {
      if (is.null(par.N$lambda)) stop("par.N doit contenir un parametre lambda")
      N <- rpois(n = nb, lambda = par.N$lambda*w)
    }

    names(N) <- labb
    
    n.evt <- sum(N)
    
    if (is.character(densfun.y)) {
      
      dists <- c("exponential", "weibull", "gpd")
      
      if (!densfun.y %in% dists) stop(paste("si character, densfun.y peut etre:",
                                            paste(dists, collapse = ", ")))
      
      if (densfun.y =="exponential") {
        x <- rexp(n.evt, rate = par.y$rate) + threshold
      } else if (densfun.y =="gpd") { 
        x <- rgpd(n.evt, loc = threshold, shape = par.y$shape, scale = par.y$scale)
      } else if (densfun.y =="weibull") {
        x <- rweibull(n.evt, shape = par.y$shape, scale = par.y$scale) + threshold
      } else stop("distribution de x non disponible")
      
    }
    
    block <- rep(1:nb, times = N)

    if (!is.null(labb)) 
      block <- labb[block]

   ##  if (plot) {
   ##   df <- data.frame(block = block, x = x)
   ##   plotRenouv(data = df)
   ##  }
    
    res <- list(type = "aggreg",
                N = N,
                block = block,
                x = x)

    return(res)
    
  }

}


