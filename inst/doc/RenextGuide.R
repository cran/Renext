### R code from vignette source 'RenextGuide.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(prompt = "> ", continue = "  ", width = 80)
Renext.Version <- packageVersion("Renext")


###################################################
### code chunk number 2: RenextGuide.Rnw:478-480
###################################################
library(Renext)
names(Brest)


###################################################
### code chunk number 3: RenextGuide.Rnw:484-486
###################################################
head(Brest$OTdata, n = 4)
str(Brest$OTinfo)


###################################################
### code chunk number 4: <
###################################################
End <- Brest$OTinfo$end; Start <- Brest$OTinfo$start
Dur <- as.numeric(difftime(End,  Start, units = "days"))/365.25
Dur
Dur - as.numeric(Brest$OTinfo$effDuration) 


###################################################
### code chunk number 5: summaryBrest
###################################################
class(Brest)
summary(Brest)


###################################################
### code chunk number 6: RenplotBrest
###################################################
plot(Brest)


###################################################
### code chunk number 7: missing
###################################################
head(Brest$OTmissing, n = 4)


###################################################
### code chunk number 8: labelgaronneMAX
###################################################
names(Garonne)
Garonne$MAXinfo
head(Garonne$MAXdata, n = 4)


###################################################
### code chunk number 9: RenplotGaronne
###################################################
plot(Garonne)


###################################################
### code chunk number 10: exppBrest
###################################################
expplot(x = Brest$OTdata$Surge, main = "expplot for \"Brest\"")


###################################################
### code chunk number 11: weibpBrest
###################################################
weibplot(x = Brest$OTdata$Surge-30, main = "weibplot for \"Brest\" (surge - 30)")


###################################################
### code chunk number 12: ExpPlot
###################################################
library(evd); set.seed(136)
X <- rgumbel(400); X <- X[X > 0.6]           ## X is truncated Gumbel
n <- length(X); 
Z <- sort(X); F <- (1:n)/(n+1)               ## distribution function
y.exp <- -log(1-F); y.gum <- -log(-log(F))   
plot(Z, y.exp, col = "red3", main = "exponential plot")


###################################################
### code chunk number 13: GumPlot
###################################################
plot(Z, y.gum, col = "SteelBlue3", main = "Gumbel plot")


###################################################
### code chunk number 14: spGaronne
###################################################
plot(Flow ~ date, data = Garonne$OTdata, type = "h", main = "Flows > 2500 m3/s")


###################################################
### code chunk number 15: RenextGuide.Rnw:750-751
###################################################
subset(Garonne$OTdata, date >= as.POSIXct("1945-01-01") & date <= as.POSIXct("1950-01-01"))


###################################################
### code chunk number 16: gdGaronne
###################################################
gof.date(date = Garonne$OTdata$date)


###################################################
### code chunk number 17: ieGaronne
###################################################
ie <- interevt(date = Garonne$OTdata$date)
names(ie)
d <- ie$interevt$duration
expplot(d, main = "Exponential plot for interevents")
bt <- gofExp.test(d) 
bt


###################################################
### code chunk number 18: gdBrest
###################################################
gof.Brest  <- gof.date(date = Brest$OTdata$date, skip = Brest$OTmissing,
                       start = Brest$OTinfo$start, end = Brest$OTinfo$end)
print(names(gof.Brest))


###################################################
### code chunk number 19: gdBrest
###################################################
head(gof.Brest$noskip, n = 2)


###################################################
### code chunk number 20: gdBrest2
###################################################
gof.Brest2  <- gof.date(date = Brest$OTdata$date, 
                        skip = Brest$OTmissing, plot.type = "omit",
                        start = Brest$OTinfo$start, end = Brest$OTinfo$end)


###################################################
### code chunk number 21: bp40
###################################################
data(Brest.years); data(Brest.years.missing)
bp40 <- barplotRenouv(data = Brest.years, threshold = 40,
           na.block = Brest.years.missing, main = "threshold = 40 cm")


###################################################
### code chunk number 22: bp50
###################################################
bp50 <- barplotRenouv(data = Brest.years, threshold = 50,
           na.block = Brest.years.missing, main = "threshold = 50 cm")


###################################################
### code chunk number 23: RenextGuide.Rnw:1034-1036
###################################################
bp40$tests
bp50$tests


###################################################
### code chunk number 24: RenextGuide.Rnw:1041-1042
###################################################
bp50$freq


###################################################
### code chunk number 25: feGaronne
###################################################
fit.exp <- Renouv(x = Garonne$OTdata$Flow,
                  effDuration = 65, threshold = 2500,
                  distname.y = "exponential",
                  main = "exponential")
fit.exp$estimate


###################################################
### code chunk number 26: namefitexp
###################################################
names(fit.exp)
class(fit.exp)


###################################################
### code chunk number 27: fwGaronne
###################################################
fit.weibull <- Renouv(x = Garonne$OTdata$Flow,
                      effDuration = 65, threshold = 2500,
                      distname.y = "weibull",
                      main = "Weibull")
fit.weibull$estimate
fit.weibull$sigma


###################################################
### code chunk number 28: rlGaronneOpt1
###################################################
plot(fit.weibull,  Tlim = c(1, 100), main = "return periods from 0 to 100 years")


###################################################
### code chunk number 29: rlGaronneOpt2
###################################################
plot(fit.weibull,
     Tlim = c(1, 100), ylim = c(3000, 10000),
     pct.conf = 95,
     main = "return levels and 95% limits")


###################################################
### code chunk number 30: feGaronneH
###################################################
Garonne$MAXdata$Flow


###################################################
### code chunk number 31: feGaronneH
###################################################
fit.exp.H <- Renouv(x = Garonne$OTdata$Flow,
                    effDuration = 65, threshold = 2500,
                    MAX.data=  list(Garonne$MAXdata$Flow),
                    MAX.effDuration = Garonne$MAXinfo$duration,
                    distname.y = "exponential",
                    main = "Garonne data, \"exponential\" with MAXdata")


###################################################
### code chunk number 32: fwGaronneH
###################################################
fit.weib.H <- Renouv(x = Garonne$OTdata$Flow,
                     effDuration = 65, threshold = 2500,
                     MAX.data=  list(Garonne$MAXdata$Flow),
                     MAX.effDuration = Garonne$MAXinfo$duration,
                     distname.y = "weibull",
                     main = "Garonne data, \"Weibull\" with MAXdata")


###################################################
### code chunk number 33: nonorth
###################################################
fit.exp.H$corr


###################################################
### code chunk number 34: fwGaronneObj
###################################################
fitWithObj <- Renouv(x = Garonne)


###################################################
### code chunk number 35: fwGaronneObj1
###################################################
fitWithObj1 <- Renouv(x = Garonne, threshold = 3000)


###################################################
### code chunk number 36: fixweibGaronneH
###################################################
fit.weib.fixed.H <- 
  Renouv(x = Garonne$OTdata$Flow,
         effDuration = 65, threshold = 2500,
         MAX.data = list(Garonne$MAXdata$Flow),
         MAX.effDuration = Garonne$MAXinfo$duration,
         distname.y = "weibull",
         fixed.par.y = list(shape = 1.2),
         start.par.y = list(scale = 2000),
         trace = 0,
         main = "Garonne data, \"Weibull\" with MAXdata and fixed shape")


###################################################
### code chunk number 37: nonorth
###################################################
fit.weib.fixed.H$estimate


###################################################
### code chunk number 38: fixSLTWGaronneH
###################################################
fit.SLTW.H <- 
  Renouv(x = Garonne$OTdata$Flow,
         effDuration = 65, threshold = 2500,
         MAX.data = list(Garonne$MAXdata$Flow),
         MAX.effDuration = Garonne$MAXinfo$duration,
         distname.y = "SLTW",
         fixed.par.y = list(delta = 2000, shape = 1.2),
         start.par.y = list(scale = 2000),
         main = "Garonne data, \"SLTW\" with MAXdata, delta and shape fixed")


###################################################
### code chunk number 39: nonorth
###################################################
fit.SLTW.H$cov


