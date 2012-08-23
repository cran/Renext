##============================================================
## Given a POSIXct object 'x', and a period data.frame,
## returns a logical vector ot the same length as 'x'
## and telling if the element of x falls within a period.
##
## NOTE: the periods must be clean: non-overlapping and
## in ascending order.
##
## The intervals are assumed to be <= ... < as in the
## 'findInterval' function
##============================================================

inInts <- function(x,
                   periods,
                   check = TRUE) {
  
  chgs <- timeints2bounds(periods[ , c("start", "end")], check = check)  
  period <- findInterval(x = as.numeric(x),
                         vec = as.numeric(chgs),
                         rightmost.closed = TRUE,
                         all.inside = FALSE)

  inInt <-  as.logical(period %% 2)
 
}
  
##============================================================
## Build non-overlapping time intervals in ascrending order
## The intervals 'data' is data.frame with columns 'start',
## 'end' and 'comment'; it must have one row or more.
##============================================================

cleanInt <- function(data,
                     start = NULL,
                     end = NULL,
                     trace = 0) {

  ## limited control
  if ( !("start" %in% colnames(data))
      || !("end" %in% colnames(data)) )
    stop("'data' must contain columns 'start' and 'end'")

  if ( !(is(data$start, "POSIXct"))
      || !(is(data$end, "POSIXct")) )
    stop("columns 'start' and 'end' must be of class \"POSIXct\"")

  ind <- data$start > data$end
  if ( any(ind) ) {
    warning("'data' contains rows with start > end")
    data <- data[!ind, ]
  } 
  
  inds <- order(data$start)


  if (nrow(data) > 1 ) {
  
    if (trace) cat("simplifying periods...\n")
    
    startCol <- data[inds, ]$start
    endCol <- data[inds, ]$end
    commentCol <- data[inds, ]$comment
    
    start.old <- startCol[1]
    end.old <- endCol[1]
    comment.old <- commentCol[1]
    
    n <- 0
    
    for (i in 2:length(startCol)){
      
      if (startCol[i] > end.old) {
        if (n==0) {
          start.clean <- start.old
          end.clean <- end.old
          comment.clean <- comment.old
        } else {
          start.clean <- c(start.clean, start.old)
          end.clean <- c(end.clean, end.old)
          comment.clean <- c(comment.clean, comment.old)
        }
        n <- n+1
        start.old <- startCol[i]
        end.old <- endCol[i]
        comment.old <- commentCol[i]
      } else {
        end.old <- endCol[i]
      }
      
    }

    ## build a dataframe. It can not have zero rows
    ## at this stage, but car is needed if 'start' or
    ## 'end' formal is given
    
    start.clean <- c(start.clean, start.old)
    end.clean <- c(end.clean, end.old)
    comment.clean <- c(comment.clean, comment.old)
    
    daf <- data.frame(start = as.POSIXct(start.clean),
                      end = as.POSIXct(end.clean),
                      comment = I(comment.clean))
    
    if (as.numeric(trace) >= 1) {
      print(daf)
    }
    
  } else {
     daf <- data
  }
  
  ## Cut at begining or end if necessary
  if (!is.null(start)) {
    start <- as.POSIXct(start)
    if (start > daf$start[1]) {
      ## period is a period number
      chgs <- timeints2bounds(daf[ , c("start", "end")])
      period <- findInterval(x = as.numeric(start),
                             vec = as.numeric(chgs),
                             rightmost.closed = TRUE,
                             all.inside = FALSE)
      
      ## start is in a gap... or not
      inGap <- period %% 2
      numGap <- period %/% 2 + 1
      if (trace)
        cat("restrict with 'start' = ", format(start), "inGap = ", inGap, "numGap = ", numGap, "\n")
      daf <- daf[ 1:nrow(daf) >= numGap, ]
      if ( (nrow(daf) > 0) && inGap) daf$start[1] <- start
    }
  }
  
  if (as.numeric(trace) >= 1) {
    print(daf)
  }

  if (nrow(daf) && (!is.null(end))) {
    end <- as.POSIXct(end)
    if (end < daf$end[nrow(daf)]) {
      ## period is a period number
      chgs <- timeints2bounds(daf[ , c("start", "end")])
      period <- findInterval(x = as.numeric(end),
                             vec = as.numeric(chgs),
                             rightmost.closed = TRUE,
                             all.inside = FALSE)
      ## end is in a gap... or not
      inGap <- period %% 2
      numGap <- period %/% 2 + 1 + inGap
      if (trace)
        cat("restrict with 'end' = ", format(end), "inGap = ", inGap, "numGap = ", numGap, "\n")
      daf <- daf[ 1:nrow(daf) < numGap, ]
      if ( (nrow(daf) > 0) && inGap) daf$end[nrow(daf)] <- end
    }
  }
  
  daf
  
}

##============================================================
## Returns an <event></events> node as a data.frame or NULL
## object if no <event/> chld is found.
## Note that 'date' is a POSIXct object, and that the format
## is assumed to be understandable by R.
## Other formats (such as locale specific) should be
## used within files and not within the XML file directly.
##============================================================

parse.eventsNode <- function(node,
                             varName = "Var",
                             trace = 1) {
  
  event.nodes   <- xmlElementsByTagName(node, "event")
  if (length(event.nodes)>0) {  
    date <- NULL
    Var <- NULL
    comment  <- NULL

    ## 'date' is kept as character here
    for (i in 1:length(event.nodes)) {
      node.i <- event.nodes[[i]]
      date <- c(date, xmlGetAttr(node.i, "date"))
      Var <- c(Var, as.numeric(xmlValue(node.i)))
      comment.i <- xmlGetAttr(node.i, "comment")
      if (length(comment.i) == 0) comment.i <- ""
      comment <- c(comment, comment.i)
    }
    
    ## 'date' is kept as character here. Are there true
    ## datetime in it?
    date2 <- as.POSIXct(rep(NA, length(date)))
    ind <- grep("[0-9]{1,4}-[0-1][0-9]-[0-3][0-9]($|( ([0-1][0-9]|[2][0-3]):([0-5][0-9])))",
                date, perl = TRUE)
    ## use loop since format could differ(?)
    if (length(ind)) {
      for (i in ind) date2[i] <- as.POSIXct(date[i])
    }

    Data <- data.frame(date = as.POSIXct(date2),
                       Var = as.numeric(Var),
                       comment = I(comment))
    colnames(Data)[2] <- varName
    Data
  } else {
    NULL
  }
  
}

##============================================================
## Returns an <period></periods> node as a data.frame or NULL
## object.
## Note that 'start' and 'end' are POSIXct objects and that
## their format is assumed to be understandable by R.
## Other formats (such as locale specific) should be
## used within files and not within the XML file directly.
##============================================================

parse.periodsNode <- function(node,
                              trace = 1) {
  
  period.nodes   <- xmlElementsByTagName(node, "period")
 ## year.nodes   <- xmlElementsByTagName(node, "year")
  if (length(period.nodes)>0) {
    start <- NULL
    end <- NULL
    comment  <- NULL
    
    for (i in 1:length(period.nodes)) {
      node.i <- period.nodes[[i]]
      start <- c(start, xmlGetAttr(node.i, "start"))
      end <- c(end, xmlGetAttr(node.i, "end"))
      comment.i <- xmlGetAttr(node.i, "comment")
      if (length(comment.i) == 0) comment.i <- ""
      comment <- c(comment, comment.i)
    }
    
    ## start <- as.POSIXct(start)
    ## end <- as.POSIXct(end)
    Prov <- data.frame(start = as.POSIXct(start),
                       end = as.POSIXct(end),
                       comment = as.character(comment))
    
    Data <- cleanInt(data = Prov,
                     trace = trace)
  
  } else{
    NULL
  } 
  
}

##=============================================================
## Read an 'events' file with description given in 'node'.
## The varname must be given since it is an attribute of an
## ancestor of node.
##=============================================================

read.eventsFile <- function(node,
                            dir,
                            varName = "Var",
                            trace = 1) {
  
  path <- file.path(dir, xmlGetAttr(node, "path"))
  
  if (trace) cat("Reading the data file\n", path,"...\n")
  nbCols <- as.integer(xmlGetAttr(node, "nbCols"))
  dtCol <- as.integer(xmlGetAttr(node, "dtCol"))
  varCol <- as.integer(xmlGetAttr(node, "varCol"))
  commentCol <- as.integer(xmlGetAttr(node, "commentCol"))
  if (length(commentCol) == 0) commentCol <- 0
  
  colClasses <- rep("NULL", nbCols)
  colClasses[dtCol] <- "character"
  colClasses[varCol] <- "numeric"
  if (commentCol > 0) colClasses[commentCol] <- "character"

  colNames <- paste("V", 1:nbCols)
  dtLab <- xmlGetAttr(node, "dtLab")
  colNames[dtCol] <- dtLab
  colNames[varCol] <- varName
  if (commentCol > 0) colNames[commentCol] <- "comment"
  
  readData <- read.csv(file = path,
                       header = FALSE,
                       sep = xmlGetAttr(node, "sep"),
                       skip = as.integer(xmlGetAttr(node, "skip")),
                       colClasses = colClasses,
                       col.names= colNames)
    
  dtFormat <- xmlGetAttr(node, "dtFormat")
  
  ## Now re-arrange datetime
  readData[ , dtCol] <- as.POSIXct(strptime(readData[ , dtCol],
                                            format = dtFormat))
  
  if (commentCol == 0) 
    readData <- data.frame(readData,
                           comment = I(rep("", nrow(readData))))
  
  if (trace) {
    cat("events data read from file\n")
    if (nrow(readData) > 8L){
      print(head(readData, n = 4L))
      cat(" ... <more lines>\n")
      print(tail(readData, n = 4L))
    } else print(readData)
    
  }
    
  readData
  
}

##===================================================================
## Read a 'periods' file with desctiption in 'node'. The file path is
## in the 'path' attribute of 'node', and must be located in the
## directory given in 'dir'.
##
## Note that 'start' and 'end' must be in the same format
##===================================================================

read.periodsFile <- function(node,
                             dir,
                             trace = 1) {
  
  path <- file.path(dir, xmlGetAttr(node, "path"))
  
  if (trace) cat("Reading the data file\n", path,"...\n")
  nbCols <- as.integer(xmlGetAttr(node, "nbCols"))
  startCol <- as.integer(xmlGetAttr(node, "startCol"))
  endCol <- as.integer(xmlGetAttr(node, "endCol"))
  commentCol <- as.integer(xmlGetAttr(node, "commentCol"))
  if (length(commentCol) == 0) commentCol <- 0

  dtFormat <- xmlGetAttr(node, "dtFormat")
  
  colClasses <- rep("NULL", nbCols)
  colClasses[startCol] <- "character"
  colClasses[endCol] <- "character"
  if (commentCol > 0) colClasses[commentCol] <- "character"

  colNames <- paste("V", 1:nbCols)
  dtLAb <- xmlGetAttr(node, "dtLab")
  colNames[startCol] <- "start"
  colNames[endCol] <- "end"
  if (commentCol > 0) colNames[commentCol] <- "comment"
  
  readData <- read.csv(file = path,
                       header = FALSE,
                       sep = xmlGetAttr(node, "sep"),
                       skip = xmlGetAttr(node, "skip"),
                       colClasses = colClasses,
                       as.is = TRUE,
                       col.names= colNames)
  
  ## Now re-arrange datetimes
  readData[ , "start"] <- as.POSIXct(strptime(readData[ , "start"], format = dtFormat))
  readData[ , "end"] <- as.POSIXct(strptime(readData[ , "end"], format = dtFormat))
  
  if (trace) {
    cat("events data read from file\n")
    if (nrow(readData) > 8L){
      print(head(readData, n = 4L))
      cat(" ... <more lines>\n")
      print(tail(readData, n = 4L))
    } else print(readData)
    
  }

  readData
  
}

##====================================================================
## Read or Parse 'events' as specified in 'node' which should have
## children only of type "events" or "file". Anyway, other children
## will not be taken into consideration. The node given in 'node'
## formal is typically of type "data".
##
## The varname must be given since it is an attribute of an ancestor
## of node and not of the node itself.
##====================================================================

readOrParse.events <- function(node,
                               dir,
                               varName = "Var",
                               trace = 1) {

  events <- list()
  nevents <- 0
  ## file node???
  file.nodes   <- xmlElementsByTagName(node, "file")
  
  if (length(file.nodes) > 1) 
    stop("at the time a 'data' can not contain more than one 'file' child node")
  
  if (length(file.nodes) > 0) {
    
    for (i in  length(file.nodes) ) {
      nevents <- nevents + 1
      events.i <- read.eventsFile(node = file.nodes[[i]],
                                  dir = dir,
                                  varName = varName,
                                  trace = trace)
      if (nevents > 1) events <- rbind(events, events.i)
      else  events <- events.i
    }
  }
  
  ## events node???
  events.nodes   <- xmlElementsByTagName(node, "events")
  
  if (length(events.nodes) > 0) {
    for (i in  length(events.nodes) ) {
      nevents <- nevents + 1
      events.i <- parse.eventsNode(node = events.nodes[[i]],
                                   varName = varName,
                                   trace = trace)
      if (nevents > 1) events <- rbind(events, events.i)
      else  events <- events.i
    }
  }

  ## returns a data.frame or NULL
  if (nevents) events
  else NULL
}

##=================================================================
## Read or Parse 'periods' as specified in 'node' which should be
## have children only of type "periods" or "file". Anyway, other
## children will not be taken into consideration. The node given
## in 'node' is typically of type "data".
##=================================================================

readOrParse.periods <- function(node,
                                dir,
                                trace = 1) {
  periods <- list()
  nperiods <- 0
  ## file node???
  file.nodes   <- xmlElementsByTagName(node, "file")
  
  if (length(file.nodes) > 1) 
    stop("at the time a 'data' can not contain more than one 'file' child node")
  
  if (length(file.nodes) > 0) {
    for (i in  length(file.nodes) ) {
      nperiods <- nperiods + 1
      periods.i <- read.periodsFile(node = file.nodes[[i]],
                                    dir = dir,
                                    trace = trace)
      if (nperiods > 1) periods <- rbind(periods, periods.i)
      else  periods <- periods.i
    }
  }
  
  ## periods node???
  periods.nodes   <- xmlElementsByTagName(node, "periods")
  
  if (length(periods.nodes) > 0) {
    for (i in  length(periods.nodes) ) {
      nperiods <- nperiods + 1
      periods.i <- parse.periodsNode(node = periods.nodes[[i]],
                                               trace = trace)
      if (nperiods > 1) periods <- rbind(periods, periods.i)
      else  periods <- periods.i
    }
  }
  
  if (nperiods) periods
  else NULL

}

##==============================================================
## The readXML funtion reads heterogenous data from an XML
## repository
## 
##                  
##===============================================================

readXML <- function(name,
                    dir,
                    index = "index.xml",
                    trace = 0) {

  ##----------------------------------------------------
  ## check that the file exists and can be read
  ##----------------------------------------------------
  if ((file.access(names = dir, mode = 0)!=0) ||
      (!file.info(dir)$isdir)) stop("dir =", dir, "is not an existing directory")
  
  file <- file.path(dir, index)
  
  if (trace) cat("Reading 'index' file\n", file, "...\n")

  if (file.access(names = file, mode = 4)!=0) {
    if (file.access(names = file, mode = 0)!=0) stop("'index' file not found")
    else stop("not allowed to open file")
  }
  
  arbre   <- xmlTreeParse(file)
  racine  <- xmlRoot(arbre)
  
  dataset.node   <- xmlElementsByTagName(racine, "dataset")
  dataset.names <- sapply(dataset.node, xmlGetAttr, "name")
  
  if (trace) {
    cat("datasets declared\n")
    print(dataset.names)
  }
  
  m <- match(dataset.names, name)

  if (trace) {
    cat("matching...\n")
    print(m)
  }

  nm <- sum(m == 1, na.rm = TRUE)
  
  if (nm==0) stop("no correct dataset")
  else if (nm>1) stop("several datasets match 'name'")
  
  ind <- (1:length(dataset.node))[!is.na(m) & m==1]
  dataset.node <- dataset.node[[ind]]

  ##=============================================
  ## Information about the dataset
  ##=============================================

  info <- list()
  for (attr  in c("name", "shortLab", "longLab", "varName",
                  "varShortLab", "varUnit")) {
    info[[attr]] <- xmlGetAttr(dataset.node, attr)
  }

  describe.node   <- xmlElementsByTagName(dataset.node, "describe")
  if (length(describe.node) == 0) describe.node <- NULL
  else describe.node   <- as(describe.node[[1]], "character")

  ##====================================================
  ## "Over Threshold" data
  ## The dataset must contain EXACTLY ONE OTdata node.
  ##====================================================
  
  OTdata.nodes   <- xmlElementsByTagName(dataset.node, "OTdata")
  if ( (length(OTdata.nodes) > 1) || (length(OTdata.nodes) == 0 ) )
    stop("a 'dataset' node must contain exactly one 'OTdata' child node")
  
  OTinfo <- list()
  for (attr in c("start", "end")) 
    OTinfo[[attr]] <- as.POSIXct(xmlGetAttr(OTdata.nodes[[1]], attr))
 
  for (attr in c("effDuration", "threshold"))
    OTinfo[[attr]] <- as.numeric(xmlGetAttr(OTdata.nodes[[1]], attr))

  data.nodes   <- xmlElementsByTagName(OTdata.nodes[[1]], "data")
  if ( length(data.nodes) != 1 ) 
    stop("an 'OTdata' node must contain exactly one 'data' child node")

  OTdata <- readOrParse.events(node = data.nodes[[1]],
                               dir = dir, varName = info$varName,
                               trace = trace)

  ##print(OTdata)
  
  if (any(OTdata$date < OTinfo$start)) {
    warning("'OTdata' contain events with date < OTinfo$start. They are removed.")
    OTdata <- OTdata[OTdata$date >= OTinfo$start, ]
  }
    
  if (any(OTdata[ , "date"] > OTinfo$end)) {
    warning("'OTdata' contain events with date > OTinfo$end. They are removed. ")
    OTdata <- OTdata[OTdata$date <= OTinfo$end, ]
  }
    
  if (any(OTdata[ , info$varName] < OTinfo$threshold))
      warning("'OTdata' column '", info$varName, "' contains a value < info$threshold")
  
  ##=============================================================
  ## Missing periods for OTdata
  ## We must find 'start' and 'end' for the successive periods,
  ## and caution is needed for successive periods possibly
  ## collapsing.
  ##=============================================================
  
  missing.nodes   <- xmlElementsByTagName(OTdata.nodes[[1]], "missing")
  
  if( length(missing.nodes) > 1)  
    stop("an \"OTdata\" node can not have more than one child  \"missing\"")

  if( length(missing.nodes) == 1) {

    OTmissing <- readOrParse.periods(node = missing.nodes[[1]],
                                     dir = dir,
                                     trace = trace)
    
    OTmissing <- cleanInt(data= OTmissing,
                          start = OTinfo$start,
                          end = OTinfo$end)
    
    ## warn if some OT events fall into an OTmissing period
    inGaps <- inInts(x = OTdata$date,
                     periods = OTmissing)
    
    if (any(inGaps)) {
      warning(sum(inGaps),"'OTdata' events are in 'OTmissing' periods. They are removed.")
      OTdata <- OTdata[inGaps, ]
    }

    ## Check effective duration
    noskip <- skip2noskip(skip = OTmissing,
                        start = OTinfo$start,
                        end = OTinfo$end)
    
    checkDuration <-
      sum(as.numeric(difftime(noskip$end, noskip$start, units = "day") / 365.25))

    if (abs(checkDuration - OTinfo$effDuration) > 0.1)
      warning("'effDuration' attributes does not agrre with computation")
    
    if (trace) cat("checkDuration", checkDuration, "effDuration", OTinfo$effDuration, "\n")
    
  } else {
    checkDuration <-
       sum(as.numeric(difftime(OTinfo$end, OTinfo$start, units = "day") / 365.25))

    if (abs(checkDuration - OTinfo$effDuration) > 0.1)
      warning("'effDuration' attributes does not agrre with computation")
    
    if (trace) cat("checkDuration", checkDuration, "effDuration", OTinfo$effDuration, "\n")
    
    OTmissing <- NULL
  }
  
  ##=============================================================
  ## Prepare results list
  ##=============================================================

  MyList <- list(info = info,
                 describe = describe.node,
                 OTinfo = OTinfo,
                 OTdata = OTdata,
                 OTmissing = OTmissing)

  
  ##=============================================================
  ## MAX data nodes
  ##=============================================================
  
  MAXdata.nodes   <- xmlElementsByTagName(dataset.node, "MAXdata")
  MAXinfo <- list()
  
  if (length(MAXdata.nodes)) {
    
    if (trace) cat("Processing", length(MAXdata.nodes), "'MAX' (historical) data\n")
    
    MAXdata <- list()
  
    for (i in 1:length(MAXdata.nodes)){

      node <-  MAXdata.nodes[[i]]
      
      if (i == 1)  MAXinfo[["block"]] <- 1
      else MAXinfo[["block"]] <- c(MAXinfo[["block"]], i)
      
      for (attr  in c("start", "end")) {
        if (i == 1)  MAXinfo[[attr]] <- xmlGetAttr(node, attr)
        else MAXinfo[[attr]] <- c(MAXinfo[[attr]], xmlGetAttr(node, attr))
      }

      ##--------------------------------------------------------
      ## data (events)
      ## Note that there can be no events in some case.
      ## Then there is a 'data' node but with no 'events' nor
      ## 'file' child.
      ##--------------------------------------------------------
    
      data.nodes   <- xmlElementsByTagName(node, "data")
      
      if(length(data.nodes) != 1)  
        stop("a \"MAXdata\" node must have exactly one child \"data\"")

      MAXdata.i <- readOrParse.events(data.nodes[[1]],
                                      dir = dir, varName = info$varName,
                                      trace = trace)
      
      ri <- 0
      if (!is.null(MAXdata.i)) {
        ri <- nrow(MAXdata.i)
        MAXdata.i <- data.frame(block = rep(i, times = ri),
                                MAXdata.i)
      } 
      
      if (i == 1) MAXdata <- MAXdata.i
      else if (ri) MAXdata <- rbind(MAXdata, MAXdata.i)
     
      if (i == 1)  MAXinfo[["r"]] <- ri
      else MAXinfo[["r"]] <- c(MAXinfo[["r"]], ri)
      
      if (is.null(MAXdata.i))
        warning("'MAXdata' ", i, " contains no events!")
      if ( !is.null(MAXdata.i) &&
          any(MAXdata.i$date[!is.na(MAXdata.i$date)] < as.POSIXct(MAXinfo$start[i])))
        warning("'MAXdata' ", i, " contains an event before 'start' given in 'MAXinfo'!")
      if ( !is.null(MAXdata.i) &&
          any(MAXdata.i$date[!is.na(MAXdata.i$date)] > as.POSIXct(MAXinfo$end[i]))) {
        warning("'MAXdata' ", i, " contains an event after 'end' given in 'MAXinfo'!")
      }
      
    }
    
    ## MAXinfo <- as.data.frame(MAXinfo)
    MAXinfo$start <- as.POSIXct(MAXinfo$start)
    MAXinfo$end <- as.POSIXct(MAXinfo$end)
    MAXinfo$r <- as.integer(MAXinfo$r)
    ## Only if missing
    MAXinfo <-
      data.frame(start = as.POSIXct(MAXinfo$start),
                 end = as.POSIXct(MAXinfo$end),
                 duration = round(as.numeric(difftime(MAXinfo$end, MAXinfo$start, units= "day") /365), 
                   digits = 2))
    
    MyList$MAXinfo <- MAXinfo
    MyList$MAXdata <- MAXdata
    
  } else {
    if (trace) cat("No MAX (historical) data\n")
  }

  ##========================================================================
  ## OTS data nodes
  ##========================================================================
  
  OTSdata.nodes   <- xmlElementsByTagName(dataset.node, "OTSdata")
  
  if (length(OTSdata.nodes)) {
    
    if (trace) cat("Processing", length(OTSdata.nodes), "'OTS' (historical) data\n")
    
    OTSinfo <- list()
    
    for (i in 1:length(OTSdata.nodes)){

      node <-  OTSdata.nodes[[i]]
      
      if (i == 1)  OTSinfo[["block"]] <- 1
      else OTSinfo[["block"]] <- c(OTSinfo[["block"]], i)
      
      for (attr  in c("start", "end", "threshold")) {
        if (i == 1)  OTSinfo[[attr]] <- xmlGetAttr(node, attr)
        else OTSinfo[[attr]] <- c(OTSinfo[[attr]], xmlGetAttr(node, attr))
      }
      
      ##------------------------------------------------------
      ## data (events)
      ## Note that there can be no events in some case.
      ## Then there is a 'data' node but with no 'events'
      ## nor 'file' child.
      ##------------------------------------------------------
      data.nodes   <- xmlElementsByTagName(node, "data")
      
      if( length(data.nodes) != 1)  
        stop("a \"MAXdata\" node must have exactly one child \"data\"")

      OTSdata.i <- readOrParse.events(data.nodes[[1]],
                                      dir = dir, varName = info$varName,
                                      trace = trace)
      ri <- 0
      if (!is.null(OTSdata.i)) {
        ri <- nrow(OTSdata.i)
        OTSdata.i <- data.frame(block = rep(i, times = ri),
                                OTSdata.i)
      } 
      
      if (i == 1) OTSdata <- OTSdata.i
      else if (ri) OTSdata <- rbind(OTSdata, OTSdata.i)
     
      if (i == 1)  OTSinfo[["r"]] <- ri
      else OTSinfo[["r"]] <- c(OTSinfo[["r"]], ri)

      if ( !is.null(OTSdata.i) && any(OTSdata.i$date < as.POSIXct(OTSinfo$start[i])))
        warning("'OTSdata' ", i, " contains an event before 'start' given in 'OTSinfo'")
      if ( !is.null(OTSdata.i) && any(OTSdata.i$date > as.POSIXct(OTSinfo$end[i]))) {
        warning("'OTSdata' ", i, " contains an event after 'end' given in 'OTSinfo'")
      }

    }

    OTSinfo$start <- as.POSIXct(OTSinfo$start)
    OTSinfo$end <- as.POSIXct(OTSinfo$end)
    
    OTSinfo <- data.frame(start = as.POSIXct(OTSinfo$start),
                          end = as.POSIXct(OTSinfo$end),
                          duration = round(as.numeric(difftime(OTSinfo$end, OTSinfo$start, units = "day") /365),
                            digits = 2),
                          threshold = as.numeric(OTSinfo$threshold),
                          r = as.integer(OTSinfo$r))

    if (any(OTSinfo$threshold < OTinfo$threshold) )
        warning("'OTSdata' specified with a threshold < OTinfo$threshold")
    
    MyList$OTSinfo <- OTSinfo
    MyList$OTSdata <- OTSdata
    
  } else {
    if (trace) cat("No OTS (historical) data\n")
  }
 

  attr(MyList, "class") <- "Rendata"
  MyList
  
}

