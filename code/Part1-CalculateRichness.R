#### Part 1: BASIC PATTERNS ####
source("code/Part0-GlobalParam.R")
source("code/DivThruTimeFunctions.R")

# Do you want to restrict analysis to sites with X number of samples?
restrict <- "yes"
if (restrict =="yes"){
  minSamp <- 6
}

# Calculate richness at each site
richness<- matrix(ncol=length(allTimes), nrow=length(sites))
colnames(richness)<- allTimes
rownames(richness)<- sites

for (i in 1:length(sites)){
  sitePath<- files[match(sites[i], gsub(".gdm.data.csv", "", files))]
  dat<- read.csv(paste(pollenDir, sitePath, sep="")) #read data
  datTimes<- dat[,1] # pull out time periods
  
  if (length(datTimes)<=1){
    richness[i,]<- NA     
  }else{
    matchedTimes<- allTimes[na.omit(match(datTimes, allTimes))]
    
    minTime<- min(matchedTimes, na.rm=T)
    maxTime<- max(matchedTimes, na.rm=T)
    
    richness[i,]<- calcSiteRichness(dat, minTime, maxTime, pollenThreshold, interval) 
  }
}

#richness<- richness[,ncol(richness):1]  # switch richness dataframe to run from past to present
richness<- t(richness)

#### restrict analysis to sites with at least 6 samples ####
if (restrict == "yes"){
  pa <- richness
  pa[which(pa > 0)] <- 1
  noSamp <- colSums(pa, na.rm=T)
  richness<- richness[,which(noSamp >= minSamp)]
}

richnessMeans<- rowMeans(richness, na.rm=T)
richnessVar<- apply(richness, 1, var, na.rm=T)

# Calculate sample size
sampSize<- apply(richness, 1, sampleSize)

#### at which sites is diversity increasing or decreasing? ####
timeNeg<- -(as.numeric(rownames(richness)))

linearSlopes<- vector(length=ncol(richness)) # each site has a slope
pvals<- vector(length=ncol(richness))  # each site has a p value

for (i in 1:ncol(richness)){
  if (length(which(!is.na(richness[,i])))>2){ #if there are more than two points, fit a model and store the slope
    lin.mod <- lm(richness[,i]~timeNeg)
    plot(richness[,i]~timeNeg, xlim=c(-21000,0), pch=16, ylim=c(0, max(richness, na.rm=T)))
    abline(lin.mod)
    linearSlopes[i]<- (lin.mod$coefficients[2])  
    pvals[i]<- summary(lin.mod)$coefficients[8]
  }else{
    linearSlopes[i]<- NaN
    pvals[i]<- NaN
  }
}

significants<- which(pvals<=0.05)
positives<- which(linearSlopes>0) #these are sites that show INCREASING RICHNESS through time
negatives<- which(linearSlopes< 0) #these are sites that show DECREASING RICHNESS through time
neutrals<- which(linearSlopes==0)

sigPos<- intersect(significants, positives) #these are sites that show INCREASING RICHNESS through time
sigNeg<- intersect(significants, negatives) #these are sites that show DECREASING RICHNESS through time
nonSig<- seq(1:length(linearSlopes))[-union(sigPos, sigNeg)]

posRichnessMeans<- rowMeans(richness[,sigPos], na.rm=T)
negRichnessMeans<- rowMeans(richness[,sigNeg], na.rm=T)
nonSigRichnessMeans<- rowMeans(richness[,nonSig], na.rm=T)

richZeros <- richness
richZeros[which(is.na(richZeros))]<- 0
siteRichChanges<- -(richZeros[2:nrow(richZeros),]-richZeros[1:(nrow(richZeros)-1),])
propRichChanges<- siteRichChanges/richness[2:(nrow(richZeros)), ]

propRichMeans<- rowMeans(propRichChanges, na.rm=T)

timeStep<- -(as.numeric(rownames(siteRichChanges)))

# Determine number of sites showing >= 25% change in richness
noOverallSites <- vector(length=length(timeStep))
noPosSites <- vector(length=length(timeStep))
noNegSites <- vector(length=length(timeStep))
propPosSites <- vector(length=length(timeStep))
propNegSites <- vector(length=length(timeStep))

for (j in 1:nrow(propRichChanges)){
  noOverallSites[j] <- length(which(!is.na(propRichChanges[j, ])))
  noPosSites[j] <- length(which(propRichChanges[j, ] >= .25))
  propPosSites[j] <- noPosSites[j] / noOverallSites[j]
  noNegSites[j] <- length(which(propRichChanges[j, ] <= -.25))
  propNegSites[j] <- noNegSites[j] / noOverallSites[j]
}

## Of sites with significant change, what was the average magnitude of change?

save.image(file="workspaces/richness.RData")


#### Plotting ####
# load("workspaces/richness.RData")

#### Genus richness- overall and at individual sites ####

# Plot the mean genus richness across all sites and sample size 

pdf(file="figures/Correlation-Richness-NoSites.pdf", height=4, width=6)
plot(richnessMeans ~ sampSize, pch=16, xlab="Number of sites", ylab="Mean richness")
test1<- cor.test(richnessMeans, sampSize) #is richness correlated with sample size?
legend("bottomright", bty="n", paste("cor=",round(test1$estimate, 4), "; p=", round(test1$p.value, 4), sep=""), cex=0.75)
dev.off()

pdf(file="figures/RichnessSampSizeThruTime-all.pdf", height=4, width=6)
par(mar=c(4,4,4,4)+0.1)
plot(richnessMeans~timeNeg, type="l", xlim=c(-21000,0), xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
par(new=T)
plot(sampSize~timeNeg, pch=16, col="darkgray", ylim=c(0, max(sampSize)), xlim=c(-21000,0), axes=F, xlab="", ylab="", cex=0.75)
axis(4, at=seq(0, max(sampSize), by=100), col.axis="darkgray")
mtext("Number of Sites", 4, line=3, col="darkgray")
dev.off()

# Plot richness trajectory at each site, all together

pdf(file="figures/RichnessThruTime-all-withLines.pdf", height=4, width=6)
par(mfrow=c(1,1))
plot(richnessMeans ~ timeNeg, 
     xlim=c(-21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (i in 1:ncol(richness)){
  lines(richness[,i]~timeNeg, col="gray") 
}
lines(richnessMeans~timeNeg, col="black", lwd=2)
dev.off()

# Plot each individual site- the increasing vs decreasing vs non sig sites 

pdf(file="figures/RichnessThruTime-INCREASING.pdf", height=9, width=10)
par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
for (k in 1:length(sigPos)){
  plot(richness[,sigPos[k]]~timeNeg, 
       type="l", 
       ylim=c(0, max(richness, na.rm=T)), xlim=c(-21000,0),
       xlab="", ylab="Site Richness")
  mtext(colnames(richness)[sigPos[k]], side=3, line=-1.5, adj=1, cex=0.75)
}
dev.off()

pdf(file="figures/RichnessThruTime-DECREASING.pdf", height=9, width=10)
par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
for (k in 1:length(sigNeg)){
  plot(richness[,sigNeg[k]]~timeNeg, 
       type="l", 
       ylim=c(0, max(richness, na.rm=T)), xlim=c(-21000,0),
       xlab="", ylab="Site Richness")
  mtext(colnames(richness)[sigNeg[k]], side=3, line=-1.5, adj=1, cex=0.75)
}
dev.off()

pdf(file="figures/RichnessThruTime-nonSig.pdf", height=9, width=10)
par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
for (k in 1:length(nonSig)){
  plot(richness[,nonSig[k]]~timeNeg, 
       type="l", 
       ylim=c(0, max(richness, na.rm=T)),
       xlab="", ylab="Site Richness")
  mtext(colnames(richness)[nonSig[k]], side=3, line=-1.5, adj=1, cex=0.75)
}
dev.off()

# Plot richness trajectory at each site- same plot with sig pos vs sig neg colored differently

pdf(file="figures/RichnessThruTime-all-withPosNegLines.pdf", height=4, width=6)
par(mfrow=c(1,1))
plot(richnessMeans ~ timeNeg, 
     xlim=c(-21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(nonSig)){
  lines(richness[,nonSig[k]]~timeNeg, col="gray")
}
for (k in 1:length(sigPos)){
  lines(richness[,sigPos[k]]~timeNeg, col="red")
}
for (k in 1:length(sigNeg)){
  lines(richness[,sigNeg[k]]~timeNeg, col="blue")
}
lines(richnessMeans~timeNeg, col="black", lwd=2)
dev.off()

# Plot richness trajectory at each site- different plots for sig pos vs sig neg

pdf(file="figures/RichnessThruTime-threePanels-withPosNegLines.pdf", height=4, width=10)
par(mfrow=c(1,3))
plot(richnessMeans ~ timeNeg, 
     xlim=c(-21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(nonSig)){
  lines(richness[,nonSig[k]]~timeNeg, col="gray")
}
lines(nonSigRichnessMeans~timeNeg, col="black", lwd=2)
legend("topleft", legend="Non-significant sites", bty="n")

warmPalette<- colorRampPalette(c("red", "orange"))(43)
plot(richnessMeans ~ timeNeg, 
     xlim=c(-21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(sigPos)){
  lines(richness[,sigPos[k]]~ timeNeg, col=warmPalette[max(which(!is.na(richness[,sigPos[k]])))-1]) #this sets the color according to the time the time series begins
}
lines(posRichnessMeans~timeNeg, col="black", lwd=2)
legend("topleft", legend="Positive sites", bty="n")

coolPalette<- colorRampPalette(c("green", "blue"))(43)
plot(richnessMeans ~ timeNeg, 
     xlim=c(-21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(sigNeg)){
  lines(richness[,sigNeg[k]]~ timeNeg, col=coolPalette[max(which(!is.na(richness[,sigNeg[k]])))-1])
}
lines(negRichnessMeans~timeNeg, col="black", lwd=2)
legend("topleft", legend="Negative sites", bty="n")

dev.off()

#### Plot proportional richness change- overall and at individual sites ####

# Plot each individual line plus mean proportional change

pdf(file="figures/PropRichnessThruTime-all-withLines.pdf", height=4, width=6)
par(mfrow=c(1,1))
plot(propRichMeans ~ timeStep, 
     xlim=c(-21000,0), ylim=c(min(propRichChanges, na.rm=T), max(propRichChanges, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Proportional Richness Change")
for (i in 1:ncol(propRichChanges)){
  lines(propRichChanges[,i]~timeStep, col=rainbow(ncol(propRichChanges))[i]) 
}
lines(propRichMeans~timeStep, col="black", lwd=2)
dev.off()


# Plot only those sites with >25% change

pdf(file="figures/PropRichnessThruTime-pos and neg.pdf", height=6, width=4)
par(mfrow=c(2,1))
plot(noPosSites~timeStep,
     xlim=c(-21000,0), ylim=c(0, max(noOverallSites)), 
     pch=16, col="red", 
     xlab="Time step starting (yr BP)", ylab="Number of sites with >25% change")
points(noNegSites~timeStep,
       xlim=c(-21000,0), ylim=c(0, max(noPosSites, noNegSites)), 
       pch=16, col="blue")
points(noOverallSites~timeStep, pch=16, col="gray")

plot(propPosSites~timeStep,
     xlim=c(-21000,0), ylim=c(0, max(propPosSites, propNegSites)), 
     type="l", col="red", 
     xlab="Time step starting (yr BP)", ylab="Proportion of sites with >25% change")
points(propNegSites~timeStep,
       xlim=c(-21000,0), ylim=c(0, max(propPosSites, propNegSites)), 
       type="l", col="blue")
dev.off()

# Is there a relationship between times with lots of positive and negative change?
c<- cor.test(propNegSites, propPosSites)
pdf(file="figures/Correlation-propNeg-propPos.pdf", height=4, width=6)
  plot(propNegSites~propPosSites, pch=16)
  abline(lm(propNegSites~propPosSites))
  legend("topright", legend=paste("t=", round(c$statistic, 2), "; df=", c$parameter, "; p=", round(c$p.value, 2), sep=""), bty="n")
dev.off()

# plot the density of different rates of change with time periods
# note: density plot isn't working, or at least, density through time is 0.

pdf(file="figures/PropRichnessChangeThruTime-pos and neg.pdf", height=4, width=6)
  proportionDF<- matrix(ncol=2, nrow=0)
  for (i in 1:ncol(propRichChanges)){
    temp<- as.data.frame(cbind(as.numeric(rownames(propRichChanges)), propRichChanges[,i]))
    proportionDF <- rbind(proportionDF, temp)
  }
  colnames(proportionDF)<- c("TimeStep", "propChange")   
  proportionDF<- proportionDF[-which(is.na(proportionDF$propChange)), ]  
  plot(propChange~TimeStep, data=proportionDF, pch=16, cex=0.75, col=rgb(0,0,0,0.25))
dev.off()

#### When you control for length of record, is there a pattern to richness change? ####
startTimes <- apply(richness, 2, function(x) max(as.numeric(names(which(!is.na(x))))))

holoSites <- which(startTimes <= 11000)
pleistSites <- which(startTimes > 11000)

hist(linearSlopes[pleistSites])
hist(linearSlopes[holoSites])

hist(linearSlopes[holoSites], xlim=c(-.003, .003), ylim=c(0, 150), col=rgb(0,0,1, alpha=0.75), xlab="", ylab="", main="")
par(new=T)
hist(linearSlopes[pleistSites], xlim=c(-.003, .003), ylim=c(0, 150), col=rgb(1,1,0, alpha=0.5), axes=F, xlab="", ylab="", main="")

sigPosPleist <- intersect(sigPos, pleistSites)
sigNegPleist <- intersect(sigNeg, pleistSites)

sigPosHolo <- intersect(sigPos, holoSites)
sigNegHolo <- intersect(sigNeg, holoSites)

length(sigPosPleist)/length(pleistSites)
length(sigNegPleist)/length(pleistSites)
length(sigPosHolo)/length(holoSites)
length(sigNegHolo)/length(holoSites)

## Results ####
# > length(sigPosPleist)/length(pleistSites)
# [1] 0.4327485
# > length(sigNegPleist)/length(pleistSites)
# [1] 0.1169591
# > length(sigPosHolo)/length(holoSites)
# [1] 0.2169118
# > length(sigNegHolo)/length(holoSites)
# [1] 0.1323529
