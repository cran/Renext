##=======================================================================
## Author: Y. Deville
#"
## Takes a data.frame with 'start' and 'end' columns and
## return the POSIXct 'bounds' object
##
##    c(start[1], end[1], start[2], end[2], ..., )
##
## If check is TRUE and 'bounds' are not coherent, an error occurs
##
## At the time check must be TRUE!!!
##
##=======================================================================

timeints2bounds <- function(timeint,
                            check = TRUE) {  
  ## convert noskip to a vector start[1], end[1], start[2], end[2], ...
  nt <- nrow(timeint)
  if (nt > 0) {
    inds <- as.vector(t(matrix(1:(2*nt), nrow = nt, ncol = 2)))
    bounds <- c(timeint$start, timeint$end)[inds]
    if (check && is.unsorted(bounds, strictly = TRUE))
      stop("expecting non-overlapping time intervals in strictly ascending order")
  } else bounds <- NULL
  
  bounds
  
}

##=====================================================================
## Author: Y. Deville
##
## skip2noskip returns a data.frame with 'start' and 'end' columns
## for the NON skipped periods from the same type of data.frame
## but for skipped periods.
##
## This function was re-written at version 0.6-0
##
##=====================================================================

skip2noskip <- function(skip = NULL,
                        start = NULL,
                        end =  NULL) {
  
  if ( is.null(skip) || nrow(skip) == 0) {
    if ( is.null(start) || is.null(end) )
      stop("when 'skip' is NULL 'start' and 'end' both must be non-null")
    noskip <- data.frame(start = I(as.POSIXct(start)),
                         end = I(as.POSIXct(end)))
  } else {

    ns <- nrow(skip)
    
    if (is.null(start)) start <- as.POSIXct(skip$start[1])
    else start <- as.POSIXct(start)
    
    if (is.null(end)) end <- as.POSIXct(skip$end[ns])
    else end <- as.POSIXct(end)    

    if (start >= end) stop("'start' must be < 'end'")
    
    ## bounds[1] begins a skipped time interval
    ## pers give the period just after the bound
    bounds <- timeints2bounds(skip)  
    periods <- 1:(2*ns)
    
    ## keep the skipped interval in the period
    ind <- (bounds > start) & (bounds < end)  

    if (any(ind)) {
      
      bounds.ns <- bounds[ind]
      periods.ns <- periods[ind]
      per1 <- periods.ns[1]
      per2 <- periods.ns[length(periods.ns)]
      
      ## if 'start' is in a non skip period, add it
      ## thus the first elt in 'periods.ns' is always even
      ## that is: alway begins a non-skip period
      if (per1 %% 2 == 1) {
        bounds.ns <- c(start, bounds.ns)
        periods.ns <- c(per1-1, periods.ns)
      }
      
      ## if end is in a non skip period, dd it
      ## thus the last elt 'periods.ns' is always odd
      if (per2 %% 2 == 0) {
        bounds.ns <- c(bounds.ns, end)
        periods.ns <- c(periods.ns, per2+1)
      }
      
      ## starts are know with odd periods, and ends
      ## with even periods
      ind <- seq(from = 1, to = length(bounds.ns), by = 2)

      noskip <- data.frame(start = as.POSIXct(bounds.ns[ind]),
                           end = as.POSIXct(bounds.ns[ind+1]),
                           period = periods.ns[ind] %/% 2 )
      
    } else {
      noskip <- data.frame(start = as.POSIXct(start),
                           end = as.POSIXct(end),
                           period = 1)
    }
  }
  
  noskip

}

##========================================================================
## Author: Y. Deville
##
## Computes intervents IN DAYS corresponding to the
## given dates of evts. 
##
## If 'skip' or 'noskip' is precised, only successive dates
## in the same no skip period  are retained.
##
## This function wa rewritten at version 0.6-0 and the output was
## redefined at this circumstance.
##
##=========================================================================

interevt <- function(date,
                     skip = NULL,
                     noskip = NULL) {
  
  date <- as.POSIXct(date)
  
  if (!is.null(skip) && !is.null(noskip))
    stop("only one of the two args 'skip' and 'noskip' args can be given")
  
  ## easier to work with 'noskip' than 'skip'...
  if(!is.null(skip))
    noskip <- skip2noskip(skip = skip,
                          start = date[1]  ,
                          end = date[length(date)]  )
  
  if (is.null(noskip)) {
    ndate <- length(date)
    per  <- rep(1, ndate-1) 
    interevt <-   as.numeric(diff(date, units = "days"))
    t.start <- date[-ndate]
    t.end <- date[-1]

    ie <- data.frame(period = per,
                     start =  as.POSIXct(t.start),
                     end = as.POSIXct(t.end),
                     duration = interevt)
    
    res <- list(interevt = ie,
                noskip = NULL)
    
  } else {

    chgs <- timeints2bounds(noskip)

    if (any(diff(chgs) < 0))
      stop("non-consistent data in 'noskip' data.frame bounds must be in ascending order")
    
    ## period is a period number
    periods <- findInterval(x = as.numeric(date),
                            vec = as.numeric(chgs),
                            rightmost.closed = TRUE,
                            all.inside = FALSE)
    ip <- unique(periods)
    interevt <- numeric(0)
    per <- integer(0)
    Flag <- FALSE
    ticks <- character(0)
    daysfrom <- numeric(0) 
    cumd <- 0
    
    for (i in ip) {
      ## i should be odd
      if (i %% 2 == 1) {

        ind <- periods == i

        if (sum(ind) > 1) {
          dind <- date[ind]
          dt <- diff(dind, units = "days")
          interevt <- c(interevt, dt)
          ## Check this rule
          per <- c(per, ((periods[ind])[-1]+1) %/% 2)

          nind <- length(dind)
          y1 <- as.integer(format(dind[1], "%Y"))
          y2 <- as.integer(format(dind[nind], "%Y"))

          if ( y1 < y2) {
            ys <- seq(from = as.POSIXct(sprintf("%4d-01-01 12:00", y1+1)),
                      to = dind[nind],
                      by = "year")
            ds <- as.numeric(difftime(ys, dind[1], units = "days"))
            ticks <- c(ticks, format(ys, "%Y"))
            daysfrom <- c(daysfrom, cumd + ds)
          }
          
          if (Flag) {
            t.start <- c(t.start, dind[-nind])
            t.end <- c(t.end, dind[-1])
          } else {
            t.start <- dind[-nind]
            t.end <- dind[-1]
            Flag <- TRUE
          }
          cumd <- cumd + as.numeric(difftime(dind[nind], dind[1], units = "days"))
        }
      } 
    }

    ie <- data.frame(period = per,
                     start =  as.POSIXct(t.start),
                     end = as.POSIXct(t.end),
                     duration = interevt)
    
    res <- list(axis = list(ticks = ticks, daysfrom = daysfrom),
                interevt = ie,
                noskip = noskip)

  }

  return(res)
  
}


