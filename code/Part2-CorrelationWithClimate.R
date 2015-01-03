#### Part 2: COMPOSITIONAL CHANGE VERSUS CLIMATE CHANGE ####
load("workspaces/richness.RData")
source("code/Part0-GlobalParam.R")

# Did more compositional change occurred between periods with more rapid climate change?
library(raster)

#### site level richness and climate change matrix #### 

#read in climate velocity stacks.
# the temporalGrad layer has the magnitudes of climate change
climateDir<- "/Volumes/bloisgroup/bloislab/Data/Climate/Paleo/CCSM3_500/With_PaleoShorelines"
var<- "tmax_year_ave"

siteClimChanges<- matrix(nrow=nrow(siteRichChanges), ncol=ncol(siteRichChanges))

for (k in 1:ncol(siteRichChanges)){
  t1<- as.numeric(colnames(siteRichChanges)[k])
  t2<- t1-1000
  
  #read in climate data and extract values from shared sites
  climVeloc<- stack(paste(climateDir, "/Climate Velocity/", var, "-", t1, "-", t2, ".tif", sep=""))
  names(climVeloc)<- c("temporalGrad", "spatialGrad", "Velocity", "BRNG")
  
  siteClimChanges[,k]<- extract(climVeloc$temporalGrad, siteLocs[match(rownames(siteRichChanges), sites)])
  
}

siteClimChanges[which(is.na(siteRichChanges))]<- NaN
siteClimChanges<- siteClimChanges*1000 #(convert to total amount of clim change, not per year)

#### extract means ####
richMean<- colMeans(siteRichChanges, na.rm=T)
climMean<- colMeans(siteClimChanges, na.rm=T)




#### Mean plotting and models ####
plot(richMean~climMean, pch=16, xlab="Temperature change")
points(richMean[1:11]~climMean[1:11], col="red", pch=16)
points(richMean[12:21]~climMean[12:21], col="blue", pch=16)

climateChangeModel<- summary(lm(richMean~climMean))  # no relationship.  

#Pleistocene change [1:19]
plot(richMean[12:21]~climMean[12:21], pch=16, xlab="Temperature change")
summary(lm(richMean[12:21]~climMean[12:21]))
abline(lm(richMean[12:21]~climMean[12:21]))

#Holocene [20:42]
plot(richMean[1:11]~climMean[1:11], pch=16, xlab="Temperature change")
summary(lm(richMean[1:11]~climMean[1:11]))
abline(lm(richMean[1:11]~climMean[1:11]))

#### Site level plotting ####
richTemp<- as.vector(siteRichChanges)
climTemp<- as.vector(siteClimChanges)

plot(richTemp~climTemp, pch=16, xlab="Site temperature change", ylab="Site richness change")
abline(lm(richTemp~climTemp))
summary(lm(richTemp~climTemp))

#Plot site level and means together
pdf(file="figures/Correlation-Richness-ClimChange-3Panels.pdf", height=2, width=6)

  par(mfrow=c(1,3), mar=c(4,4,0,0)+0.1)
  
  plot(richTemp~climTemp, pch=16, xlab="Site temperature change", ylab="Site richness change")
  #abline(lm(richTemp~climTemp))
  summary(lm(richTemp~climTemp))
  
  plot(richMean[12:21]~climMean[12:21], pch=16, xlim=c(-1,8), ylim=c(-1, 1),
       xlab="Mean temperature change", ylab="Mean richness change")
  abline(lm(richMean[12:21]~climMean[12:21]))
  pModel<- lm(richMean[12:21]~climMean[12:21])
  
  plot(richMean[1:11]~climMean[1:11], pch=16, xlim=c(-1,8), ylim=c(-1, 1),
       xlab="Mean temperature change", ylab="Mean richness change")
  abline(lm(richMean[1:11]~climMean[1:11]))
  hModel<- lm(richMean[1:11]~climMean[1:11])

dev.off()


## climate change at sites ##

posRichChange<- siteRichChanges[sigPos,]
posClimChange<- siteClimChanges[sigPos,]
plot(posRichChange~posClimChange, pch=16)
summary(lm(as.vector(posRichChange)~as.vector(posClimChange)))


negRichChange<- siteRichChanges[sigNeg,]
negClimChange<- siteClimChanges[sigNeg,]
plot(negRichChange~negClimChange, pch=16)
summary(lm(as.vector(negRichChange)~as.vector(negClimChange)))

par(mfrow=c(1,1))
plot(richTemp~climTemp, pch=16, xlab="Site temperature change", ylab="Site richness change")
points(posRichChange~posClimChange, pch=16, col="blue")
points(negRichChange~negClimChange, pch=16, col="red")
