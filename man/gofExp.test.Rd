\name{gofExp.test}
\alias{gofExp.test}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Goodness-of-fit test for exponential distribution}
\description{
  Bartlett's goodness-of-fit test for exponential distribution
}
\usage{
   gofExp.test(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{

    Sample with positive values.
    
  }
}
%%\details{}
\value{
  A list with elements
  \item{statistic}{

    Statistic.

  }
  \item{p.value}{

    Critical value.

  }
}
\references{See
  
  Yagouti A., Abi-Zeid I., Ouarda, T.B.M.J. and B. \enc{Bobée}{Bobee}
  (2001), Revue de processus ponctuels et \enc{synthèse}{synthese} de
  tests statistiques pour le choix d'un type de processus \emph{Revue
  des Sciences de l'Eau}, \bold{1}, pp. 323-361.

}
\author{Yves Deville}
%%\note{ }
\seealso{

  Among other goodness-of-fit tests \code{\link{ks.test}} in the
 \code{stats} package. See \code{\link{expplot}} for a graphical
 diagnostic.

}
\examples{
 ## a sample of size 30
 x <- rexp(30)
 res <- gofExp.test(x)

 ## ns samples: p.values should look as uniform on (0, 1)
 ns <- 100
 xmat <- matrix(rexp(30*ns), nrow = ns, ncol = 30)
 p.values <- apply(xmat, 1, function(x) gofExp.test(x)$p.value)
 plot(sort(p.values), type = "p", pch = 16)
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
