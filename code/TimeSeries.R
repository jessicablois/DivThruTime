### Times series analysis ####

load("workspaces/richness.RData")
source("code/Part0-GlobalParam.R")

richTS<- ts(richness, start=-21000, end=0, deltat=500)

library(forecast)

#coeff<- as.list(NULL)
par(mfrow=c(3,3), mar=c(3,3,3,0))
for (i in 1:9){
  t<- auto.arima(richTS[,i], seasonal=F)
  plot(t$x)
  lines(fitted(t), col="red")
  acf(richTS[,i], na.action=na.pass)
  acf(residuals(t), na.action=na.pass)
  }



richTS[,1][which(!is.na(richTS[,1]))]


acf(richTS[,i], na.action=na.pass)


## GAMM example based on Simpson post ####
times<- seq(-21000, 0, by=500)
richNew<- cbind(times, richness)
rownames(richNew) <- richNew[,1]

richNew<- as.data.frame(richNew)

i=2
richTest<- richNew[,c(1,i)]
richTest<- richTest[-(na.action(na.omit(richTest[,2]))),]
colnames(richTest)[2]<- "siteDat"

## Plot the data
ylab <- expression(Genus~Richness~(21000-0))
plot(siteDat ~ times, data = richTest, type = "o", ylab = ylab, 
     main = colnames(richNew)[i])


## Fit a smoother for Year to the data
m1 <- gamm(siteDat ~ s(times, k = 10), data = richTest)
summary(m1$gam)

## look at autocorrelation in residuals:
acf(resid(m1$lme, type = "normalized"))
pacf(resid(m1$lme, type = "normalized"))

## ...so fit the AR1
m2 <- gamm(siteDat ~ s(times, k = 10), data = richTest,
           correlation = corARMA(form = ~ times, p = 1))
## ...and fit the AR2
m3 <- gamm(siteDat ~ s(times, k = 10), data = richTest,
           correlation = corARMA(form = ~ times, p = 2))

anova(m1$lme, m2$lme, m3$lme)

plot(m1$gam, residuals = TRUE, pch = 19, cex = 0.75)
with(richTest, tsDiagGamm(m1, timevar = times, observed = siteDat))
