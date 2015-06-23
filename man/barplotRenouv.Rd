\name{barplotRenouv}
\alias{barplotRenouv}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Barplot for Renouv "Over Threshold" counts}
\description{
  Barplot for "Over Threshold" counts in time blocks (usually years)
}
\usage{
 barplotRenouv(data,
               blockname = colnames(data)[1],
               varname = colnames(data)[2],
               threshold = quantile(data[, varname], 0.2),
               na.block = NULL,
               plot = TRUE,
               main = NULL, xlab = NULL, ylab = NULL,
               mono = FALSE,
               prob.theo = 0.999,
               ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{

    A dataframe object containing the variables.

  }
  \item{blockname}{

    Name of the "block" variable (column in \code{data}). This variable
    should contain integers, or be of class "factor", but with integer
    values such as year numbers.

  }
  \item{varname}{

    Name of the variable (e.g. \code{"Surge"}).

  }
  \item{threshold}{

    Only obs for which the variable exceeds \code{threshold} will be
    taken into account.

  }
  \item{na.block}{

    Values of blocks containing missing values. See the Details section.

  }
  \item{plot}{

    If \code{FALSE} tests are computed without producing any plot.

  }
  \item{main}{

    Character for main title or \code{NULL} in which case a default main
    title is used.

  }
  \item{xlab}{

    Character for x axis label or \code{NULL} in which case a default
    lab is used.

  }
  \item{ylab}{

    Character for y axis or \code{NULL} in which case a default lab is
    used.

  }
  \item{mono}{

    If \code{FALSE} barplot will have colors, else greyscale will be
    used.

  }
  \item{prob.theo}{

    The total theoretical probability corresponding to the plotted
    (theoretical) bars.}  \item{\dots}{Further args to be passed to
    \code{barplot}.

  }
}
\details{
  
   Blocks described in the \code{na.block} are omitted in the
   determination of counts. The object given in the \code{na.block} is
   coerced to character and the same is done for values of \code{block}
   before comparing them to the \code{na.block} values.  If \code{block}
   variable is of class factor with levels representing years
   (e.g. 1980, 1981, etc.) missing blocks can be specified either as
   \code{c("1980", "1981")} or as numeric \code{c(1980, 1981)}.

   For the chi-square test, counts for neighbouring frequency classes are
   grouped in order to reach a minimum frequency of \code{5} in each
   group. E.g. if we expect respectively \code{1.0}, \code{3.8} and
   \code{7.0} blocks with frequency \code{0}, \code{1} and \code{2} for
   events, the three counts are grouped in one group with frequency
   \code{1.0+3.8+7.0=11.8}. Note that this strategy of grouping is not
   unique and is likely to weaken the power of the test. Before
   grouping, the higher class theoretical probability is computed as the
   probability to obtain a count equal to or greater than the max value.
}
\value{
  A list with the following objects.
  \item{freq}{

    frequency table (matrix) giving observed and theoretical (Poisson)
    frequencies as well as a group number for the chi-square test.

  }
  \item{overdispersion}{

    the overdispersion coefficient (variance/mean ratio).

  }
  \item{disp.test}{

    a list giving results of the (over)dispersion
    test. See the reference Yagouti and al. in the \bold{References}
    section.

  }
  \item{chisq.test}{

    a list giving results for the chis-square test of goodness-of-fit to
    the Poisson distribution.

  }
  \item{tests}{a matrix with the two tests displayed in two rows.}
  For both tests, the statistic follows a chi-square distribution under the
  null hypothesis . The list of results contains the statistic
  \code{statistic}, the number of degrees of freedom \code{df} and 
  the \eqn{p}{p}-value \code{p.value}.
}
\references{

  See Yagouti A., Abi-Zeid I., Ouarda, T.B.M.J. and
  B. \enc{Bobée}{Bobee} (2001), Revue de processus ponctuels et
  \enc{synthèse}{synthese} de tests statistiques pour le choix d'un type
  de processus \emph{Revue des Sciences de l'Eau}, \bold{1},
  pp. 323-361.
  
}
\author{Yves Deville}
\note{
  
  The two tests: (over-)dispersion and chi-square have \emph{one-sided} (upper
  tail) \eqn{p}{p}-value. In other words, we do not intend to reject when
  statistics take "abnormally small" values, but only when abnormally
  large values are met.

%%  The 'expected' and 'observed' counts in the \code{freq} object do
%%  have slightly different totals since the Poisson distribution
%%  takes values that are greater than the observed ones.
   
}
\seealso{\code{\link{plot.Rendata}}}
\examples{
## na.block influence for Brest data
opar <- par(mfrow = c(2, 2))

bp1 <- barplotRenouv(data = Brest.years, threshold = 30,
         main = "missing periods ignored")
bp2 <- barplotRenouv(data = Brest.years, threshold = 30,
         na.block = 1992, main = "1992 missing")
bp3 <- barplotRenouv(data = Brest.years, threshold = 30,
         na.block = 1991:1993, main ="1991:1993 missing")
bp4 <- barplotRenouv(data = Brest.years, threshold = 30,
         na.block = Brest.years.missing, main = "all missing periods")

par(opar)

## threshold influence
opar <- par(mfrow = c(2,2))

thresh <- c(30, 35, 40, 50)

for (i in 1:length(thresh)) {
  bp  <- barplotRenouv(data = Brest.years, threshold = thresh[i],
                   na.block = Brest.years.missing,
                   main = paste("threshold =", thresh[i], "cm at Brest"))
}
par(opar)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
