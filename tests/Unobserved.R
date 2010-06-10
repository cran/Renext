library(Renext)

##=======================================================================
## In this test a sample of indepdnant excessances Over a Thershold
## is drawn at random from an exponential distribution with mean
## theta. The corresponding model is fitted.
##
## We the add the informaation that the (high) level X.U
## was never observed over aperiod of duration w.U. We control the
## theoretical relation between estimates
##
##      theta.hat.new = theta.hat.old
##                    - y.U * (lambda.hat.new/lamda.hat.old - 1)
## 
## where old is for OT only and new is for OT AND Unobserved.
##
##=========================================================================

set.seed(123)

n <- 20
theta <- 100
x <- 35 + theta*rexp(n)

x.U <- 100
w.U <- 100 + 100*runif(1)

fit.expon <- list()

opar <- par(mfrow = c(2, 1), mar = c(5, 5, 3, 1))

fit.expon[[1]] <- 
  fRenouv(x.OT = x, sumw.BOT = 50,
          distname.y = "exponential",
          threshold = 35,
          conf.pct = c(70, 95),
          prob.max = 0.99995,
          pred.period = c(10, 100, 1000),
          trace = 1,
          main = "\"exponential\" without x.U")

fit.expon[[2]] <- 
  fRenouv(x.OT = x, sumw.BOT = 50,
          x.U = x.U, w.U = w.U,                        ## Unobserved
          distname.y = "exponential",
          threshold = 35,
          conf.pct = c(70, 95),
          prob.max = 0.99995,
          pred.period = c(10, 100, 1000),
          control.H = list(fnscale = -1, trace = 3),
          trace = 1,
          main = "\"exponential\" plus unobserved x.U ")

par(opar)

estimates <- 
  array(NA,
        dim = c(2, length(fit.expon[[1]]$estimate)),
        dimnames = list(c("sans", "US"),
          names(fit.expon[[1]]$estimate)))

for (i in 1:2) estimates[i, ] <- fit.expon[[i]]$estimate

Test <- 1/estimates["US", "rate"] -
  1/estimates["sans", "rate"] -
  (x.U - 35) * (estimates["US", "lambda"]/estimates["sans", "lambda"] -1) 

## L'erruer doit être peite, disons ne pas dépasser 1%
RelError <- abs(Test)*estimates["sans", "rate"]

stopifnot(RelError < 1.0)
