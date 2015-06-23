\name{gof.date}
\alias{gof.date}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Goodness-of-fit for the distribution of dates}
\description{
  Goodness-of-fit diagnostics for the distribution of event dates
  in a (assumed) Poisson process
}
\usage{
    gof.date(date,
             start = NULL,
             end = NULL,
             plot = TRUE,
             main = NULL,
             skip = NULL,
             plot.type = "skip")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{date}{

    Object of class \code{POSIXct} (or that can be coerced to this
    class) giving the dates to be tested. Must be in strictly increasing
    order.

  }
  \item{start}{

    The beginning of the interval, a \code{POSIXct} object. If
    \code{NULL}, the first event in \code{date} is used.

  }
  \item{end}{

    Object of class \code{POSIXct} the end of the interval. If
    \code{NULL}, the last event in \code{date} is used.

  }
  \item{plot}{

    Should a plot be shown?

  }
  \item{main}{

    Character giving the main title of the plot. The default \code{NULL}
    stands for a default main describing the period.

  }
  \item{skip}{

    Optional data.frame with columns \code{start} and \code{end}
    indicating start and end of skipped periods. The two columns need to
    be coerced to POSIXct objects. They can be POSIXct or character with
    POSIX datetime format.

  }
  \item{plot.type}{

    Character indicating the type of plot to produce when a \code{skip}
    data.frame is given. With \code{plot.type = "skip"} the plot shows
    missing periods as greyed rectangles and the displays the results of
    a Kolmogorov-Smirnov (KS) test performed on the events. For the
    \code{"omit"} case the missing periods are collapsed into vertical
    lines on the plot and the displayed results are for an "effective"
    KS test of uniformity performed omitting the missing periods.

  }
}
\details{

  In the homogeneous Poisson process, events occur on a time interval in
  a uniform fashion. More precisely, for a given time interval the
  distribution of the event dates conditional to their number \eqn{n}{n} is the
  distribution of the order statistics of a sample of size \eqn{n}{n} of
  the uniform distribution on this interval.
  
  When the interval has limits taken at events the uniformity statement
  remains true, but for \emph{inner} events. This behaviour is met when
  \code{start} and \code{end} are not given and taken as the first and
  last events in \code{date}.

}
\value{
  A list
  \item{effKS.statistic, KS.statistic}{

    Kolmogorov-Smirnov global test statistic for uniformity (bilateral
    test) omitting slipped periods or not.

  }
  \item{effKS.pvalue, KS.pavalue}{

    Critical probability in the KS test omitting skipped periods or not.

  }
  \item{effnevt, nevt}{

    Number of events omitting skipped periods or not.

  }
  \item{effduration, duration}{

    Effective duration i.e. total duration of non-skipped periods. In
    years, omitting skipped periods or not.

  }
  \item{effrate, rate}{

    Occurrence rate in \bold{number of events by year}, omitting skipped
      periods or not.

  }
  \item{effduration, duation}{

    Total duration in \bold{years}, omitting missing periods or not.

  }
  \item{noskip}{

    Data.frame object giving indications on the periods that are NOT
    skipped over (hence usually non-missing periods). These are :
    \code{start}, \code{end} (POSIX), \code{duration} (in years)
    \code{rate} (in number of events by year) and Kolmogorov test
    statistic and p-value. This data.frame is only available when a
    suitable \code{skip} has been given.

  }

  When the number of events corresponding to the indications of args is
  \code{0}, the function returns \code{NULL} with a warning. When the
  number of events is less than \code{6} a warning is shown.

}
%%\references{ }
\author{

  Yves Deville

}
\note{
  
  In practical contexts missing periods are often met in the
  datasets. The diagnostic should therefore be applied on \emph{every
  period with no missing data}. Even if the event dates seem reasonably
  uniform, it is a good idea to check that the rates do not differ
  significantly over intervals.

  When some events are missing and no suitable information is given via
  the \code{skip} argument, the global \code{rate}, \code{KS.statistic}
  and \code{KS.pvalue} are of little interest. Yet the graph might be
  instructive.

}
\section{Warning}{
  
 When skipped periods exist the number of events, duration, rate the
 global KS test must be computed by omitting the skipped periods in the
 duration and retaining only valid interevents. The indication given in
 \code{nevt} \code{rate} and \code{duration} should be used only when no
 skipped period exist (\code{skip = NULL} on input) and replaced by
 \code{effnevt}, \code{effrate} and \code{effduration} otherwise.

}
\seealso{

  \code{\link{interevt}} function for the determination of interevents
  ans subsequent diagnostics.

}
\examples{
## Use "Brest" dataset
## simple plot. Kolmogorov-Smirnov is not useful
gof1 <- gof.date(date = Brest$OTdata$date)

## consider missing periods. Much better!
gof2 <- gof.date(date = Brest$OTdata$date,
         skip = Brest$OTmissing,
         start = Brest$OTinfo$start,
         end = Brest$OTinfo$end)

print(gof2$noskip)

## Second type of graph
gof3 <- gof.date(date = Brest$OTdata$date,
         skip = Brest$OTmissing,
         start = Brest$OTinfo$start,
         end = Brest$OTinfo$end,
         plot.type = "omit")

## non-skipped periods at Brest
ns <- skip2noskip(skip = Brest$OTmissing,
                 start = Brest$OTinfo$start,
                 end = Brest$OTinfo$end)

## say 9 plots/diagnostics
oldpar <- par(mar = c(3, 4, 3, 2), mfcol = c(3, 3))

for (i in 1:9) {
  GOF <- gof.date(date = Brest$OTdata$date,
           start = ns$start[i],
           end = ns$end[i])
}

par(oldpar)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
