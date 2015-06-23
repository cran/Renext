\name{parDeriv}
\alias{parDeriv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Derivation of probability functions with respect to the parameters
}
\description{
  Derivation of probability functions with respect to the parameters
  by using closed forms.
}
\usage{
parDeriv(par, x, distname, sum = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{

    Vector of parameter values.

  }
  \item{x}{

    Observations or data at which the derivatives are to be computed.

  }
  \item{distname}{
    
    Name of the distribution. See \bold{Details}.

  }
  \item{sum}{
    
    Logical. If \code{TRUE}, a summation over the element of \code{x} is
    carried. Otherwise, the first dimension of the result corresponds to
    the elements of \code{x}.
    
  }
}
\details{
  
  Only a few distributions are and will be available. For now, these are:
  the two-parameter Weibull \code{c("shape", "scale")}, the
  two-parameter Generalised Pareto, \code{c("scale", "shape")} and
  the two-parameter Lomax and maxlo distributions.

}
\value{
  
  A list of arrays containing the first and second order derivatives.
  \item{derLogdens, der2Logdens}{Derivatives of the log-density
  \eqn{\log f(x)}{log f(x)}.}  \item{derSurv, der2Surv}{Derivatives of
  the survival function \eqn{S(x)}.}  When \code{x} has length \eqn{n}
  and the distribution depends on \eqn{p} parameters, the arrays of
  first and second order derivatives have dimension \eqn{n \times
  p}{n*p} and \eqn{n \times p \times p}{n*p*p} when \code{sum} is
  \code{FALSE}. If \code{sum} is \code{TRUE} the summation drops the
  first dimension and the arrays are \eqn{p}{p} and \eqn{p \times
  p}{p*p}.
 
}
\references{
  See the \emph{Renext Computing Details} document.
}
\author{
  Yves Deville
}
\seealso{\code{\link{Maxlo}} and \code{\link{Lomax}}.}
\examples{
set.seed(1234)
distname <- "maxlo"
if (distname == "weibull") {
    logL <- function(par) {
        sum(dweibull(x, shape = par["shape"], scale = par["scale"], log = TRUE))
    }
    sumS <- function(par) {
        sum(pweibull(x, shape = par["shape"], scale = par["scale"],
                     lower.tail = FALSE))
    }
    pars <- c("shape" = rexp(1), "scale" = 1000 * rexp(1))
    x <- rweibull(n = 100, shape = pars["shape"], scale = pars["scale"])
    Der <- parDeriv(par = pars, x = x, distname = "weibull") 
} else if (distname == "gpd") {
    require(evd)
    logL <- function(par) {
        sum(dgpd(x, loc = 0, shape = par["shape"], scale = par["scale"],
                 log = TRUE))
    }
    sumS <- function(par) { 
        sum(pgpd(x, loc = 0, shape = par["shape"], scale = par["scale"],
                 lower.tail = FALSE))
    }
    pars <- c("scale" = 1000 * rexp(1),
              "shape" = runif(1, min = -0.4, max = 0.4))
    x <- rgpd(n = 100, loc = 0, shape = pars["shape"], scale = pars["scale"])
    Der <- parDeriv(par = pars, x = x, distname = "gpd")
} else if (distname == "lomax") {
    logL <- function(par) {
        sum(dlomax(x, shape = par["shape"], scale = par["scale"], log = TRUE))
    }
    sumS <- function(par) { 
        sum(plomax(x, shape = par["shape"], scale = par["scale"],
                   lower.tail = FALSE))
    }
    pars <- c( "shape" = 1 + rexp(1), "scale" = 1000 * rexp(1))
    x <- rlomax(n = 100, shape = pars["shape"], scale = pars["scale"])
    Der <- parDeriv(par = pars, x = x, distname = "lomax") 
} else if (distname == "maxlo") {
    logL <- function(par) {
        sum(dmaxlo(x, shape = par["shape"], scale = par["scale"], log = TRUE))
    }
    sumS <- function(par) { 
        sum(pmaxlo(x, shape = par["shape"], scale = par["scale"],
                   lower.tail = FALSE))
    }
    pars <- c( "shape" = 2.5 + runif(1), "scale" = 100 * rexp(1))
    x <- rmaxlo(n = 100, shape = pars["shape"], scale = pars["scale"])
    Der <- parDeriv(par = pars, x = x, distname = "maxlo") 
}

## check logdens
H <- numDeriv::hessian(func = logL, x = pars)
colnames(H) <- names(pars)
Grad <- numDeriv::grad(func = logL, x = pars)

cat("gradient for log density\n")
print(cbind(parDeriv = Der$derLogdens, num = Grad))

cat("hessian for log density\n")
print(cbind(exact = Der$der2Logdens, num = H))

## check survival
HS <- numDeriv::hessian(func = sumS, x = pars)
HS <- (HS + t(HS))/2
colnames(HS) <- names(pars)
GradS <- numDeriv::grad(func = sumS, x = pars)

cat("gradient for Survival\n")
print(cbind(parDeriv = Der$derSurv, num = GradS))

cat("hessian for Survival\n")
print(cbind(exact = Der$der2Surv, num = HS))
}
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
