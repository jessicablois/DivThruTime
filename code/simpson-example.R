### Simpson blog posts about examining trends through time ####

# First post: http://www.fromthebottomoftheheap.net/2011/06/11/global-warming-since-1995-now-significant/ ####
URL <- url("http://www.cru.uea.ac.uk/cru/data/temperature/HadCRUT3v-gl.dat")
gtemp <- read.table(URL, fill = TRUE)

## Don't need the even rows --- perhaps do as case weights?
gtemp <- gtemp[-seq(2, nrow(gtemp), by = 2), ]

## set the Year as rownames
rownames(gtemp) <- gtemp[,1]

## Add colnames
colnames(gtemp) <- c("Year", month.abb, "Annual")

## Data for 2011 incomplete so work only with 1850-2010 data series
gtemp <- gtemp[-nrow(gtemp), ]
head(gtemp)

par(mfrow=c(1,1))
ylab <- expression(Temperature~Anomaly~(1961-1990)~degree*C)
plot(Annual ~ Year, data = gtemp, type = "o", ylab = ylab,
main = "Global mean temperature anomaly 1850-2010")

ylab <- expression(Temperature~Anomaly~(1961-1990)~degree*C)
plot(Annual ~ Year, data = gtemp, type = "o", ylab = ylab,
main = "Global mean temperature anomaly 1850-2010")

grecent <- subset(gtemp, subset = Year >= 1995, select = c(Year, Annual))

## plot
plot(Annual ~ Year, data = grecent, type = "o", ylab = ylab,
main = "Global mean temperature anomaly 1995-2010")

## fit a linear model through these recent data
gm1 <- lm(Annual ~ Year, data = grecent)
gm0 <- lm(Annual ~ 1, data = grecent)
anova(gm0, gm1)
coef(gm1)

acf(resid(gm1), main = "Residuals of linear model")

layout(matrix(1:4, ncol = 2))
plot(gm1)
layout(1)


require(nlme)
gg0 <- gls(Annual ~ 1, data = grecent, method = "ML")
gg1 <- gls(Annual ~ Year, data = grecent, method = "ML")
anova(gg0, gg1)

gg1 <- update(gg1, method = "REML")
gg2 <- gls(Annual ~ Year, data = grecent,
           correlation = corARMA(form = ~ Year, p = 1), method = "REML")
gg3 <- gls(Annual ~ Year, data = grecent,
           correlation = corARMA(form = ~ Year, p = 2), method = "REML")

anova(gg1, gg2, gg3)
anova(gg1, gg3)

confint(gg1)

acf(resid(gg3, type = "normalized"))
summary(gg3)
intervals(gg3)

pred <- data.frame(Year = 1995:2010)
pred <- transform(pred, yhat = predict(gg3, newdata = pred))
with(pred, yhat)


layout(matrix(1:2, ncol = 2))
## plot full data
plot(Annual ~ Year, data = gtemp, type = "o", ylab = ylab, cex.main=0.75,
     main = "Global mean temperature anomaly 1850-2010")
lines(yhat ~ Year, data = pred, lwd = 2, col = "red")
## plot the 1995-2010 close-up
plot(Annual ~ Year, data = grecent, type = "o", ylab = ylab, cex.main=0.75,
     main = "Global mean temperature anomaly 1995-2010")
lines(yhat ~ Year, data = pred, lwd = 2, col = "red")
layout(1)


# Second post: http://www.fromthebottomoftheheap.net/2011/06/12/additive-modelling-and-the-hadcrut3v-global-mean-temperature-series/ ####
## load the packages and code we need
require(mgcv)
require(nlme)
## load custom functions
tmp <- tempfile()
download.file("https://github.com/gavinsimpson/random_code/raw/master/derivFun.R", tmp)
source(tmp)
tmp <- tempfile()
download.file("https://github.com/gavinsimpson/random_code/raw/master/tsDiagGamm.R", tmp)
source(tmp)

## Global temperatures
URL <- url("http://www.cru.uea.ac.uk/cru/data/temperature/HadCRUT3v-gl.dat")
gtemp <- read.table(URL, fill = TRUE)
## Don't need the even rows
gtemp <- gtemp[-seq(2, nrow(gtemp), by = 2), ]
## set the Year as rownames
rownames(gtemp) <- gtemp[,1]
## Add colnames
colnames(gtemp) <- c("Year", month.abb, "Annual")
## Data for 2011 incomplete so work only with 1850-2010 data series
gtemp <- gtemp[-nrow(gtemp), ]
## Plot the data
ylab <- expression(Temperature~Anomaly~(1961-1990)~degree*C)
plot(Annual ~ Year, data = gtemp, type = "o", ylab = ylab)

## Fit a smoother for Year to the data
m1 <- gamm(Annual ~ s(Year, k = 20), data = gtemp)
summary(m1$gam)

## look at autocorrelation in residuals:
acf(resid(m1$lme, type = "normalized"))
## ...wait... look at the plot, only then do...
pacf(resid(m1$lme, type = "normalized"))
## seems like like sharp cut-off in ACF and PACF - AR terms probably best

## ...so fit the AR1
m2 <- gamm(Annual ~ s(Year, k = 30), data = gtemp,
           correlation = corARMA(form = ~ Year, p = 1))
## ...and fit the AR2
m3 <- gamm(Annual ~ s(Year, k = 30), data = gtemp,
           correlation = corARMA(form = ~ Year, p = 2))

anova(m1$lme, m2$lme, m3$lme)

plot(m2$gam, residuals = TRUE, pch = 19, cex = 0.75)
with(gtemp, tsDiagGamm(m2, timevar = Year, observed = Annual))

plot(Annual ~ Year, data = gtemp, type = "p", ylab = ylab)
pdat <- with(gtemp,
             data.frame(Year = seq(min(Year), max(Year),
                                   length = 200)))
p1 <- predict(m1$gam, newdata = pdat)
p2 <- predict(m2$gam, newdata = pdat)
lines(p1 ~ Year, data = pdat, col = "red")
lines(p2 ~ Year, data = pdat, col = "blue")
legend("topleft",
       legend = c("Uncorrelated Errors","AR(1) Errors"),
       bty = "n", col = c("red","blue"), lty = 1)

m2.d <- DerivMod(m2, n = 200)
plot(m2.d, sizer = TRUE, alpha = 0.01)


plot(Annual ~ Year, data = gtemp, type = "p", ylab = ylab)
lines(p2 ~ Year, data = pdat)
CI <- confint(m2.d, alpha = 0.01)
S <- signifD(p2, m2.d$Year$deriv, CI$Year$upper, CI$Year$lower,
             eval = 0)
lines(S$incr ~ Year, data = pdat, lwd = 3, col = "blue")
lines(S$decr ~ Year, data = pdat, lwd = 3, col = "red")
