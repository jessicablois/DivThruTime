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

richness<- richness[,ncol(richness):1]  # switch richness dataframe to run from past to present
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
timeKYR<- rev(allTimes/1000)
richMat<- cbind(timeKYR, richness)

linearSlopes<- vector(length=ncol(richness)) # each site has a slope
pvals<- vector(length=ncol(richness))  # each site has a p value

for (i in 2:ncol(richMat)){
  if (length(which(!is.na(richMat[,i])))>2){ #if there are more than two points, fit a model and store the slope
    lin.mod <- lm(richMat[,i]~richMat[,1])
    plot(richMat[,i]~richMat[,1], xlim=c(21,0), pch=16, ylim=c(0, max(richMat[,-1], na.rm=T)))
    abline(lin.mod)
    linearSlopes[i-1]<- -(lin.mod$coefficients[2])  #need to make this the opposite of original to have time running forward
    pvals[i-1]<- summary(lin.mod)$coefficients[8]
  }else{
    linearSlopes[i-1]<- NaN
    pvals[i-1]<- NaN
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
siteRichChanges<- -(richZeros[1:(nrow(richZeros)-1),]-richZeros[2:nrow(richZeros),])
propRichChanges<- siteRichChanges/richness[1:(nrow(richZeros)-1), ]

propRichMeans<- rowMeans(propRichChanges, na.rm=T)

timeStep<- as.numeric(rownames(siteRichChanges))

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
plot(richnessMeans~timeKYR, type="l", xlim=c(21,0), xlab="Time slice (kyr BP)", ylab="Mean Genus Richness")
par(new=T)
plot(sampSize~timeKYR, pch=16, col="darkgray", ylim=c(0, max(sampSize)), xlim=c(21,0), axes=F, xlab="", ylab="", cex=0.75)
axis(4, at=seq(0, max(sampSize), by=100), col.axis="darkgray")
mtext("Number of Sites", 4, line=3, col="darkgray")
dev.off()

# Plot richness trajectory at each site, all together

pdf(file="figures/RichnessThruTime-all-withLines.pdf", height=4, width=6)
par(mfrow=c(1,1))
plot(richnessMeans ~ timeKYR, 
     xlim=c(21,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (kyr BP)", ylab="Mean Genus Richness")
for (i in 1:ncol(richness)){
  lines(richness[,i]~timeKYR, col=rainbow(nrow(richness))[i]) 
}
lines(richnessMeans~timeKYR, col="black", lwd=2)
dev.off()

# Plot each individual site- the increasing vs decreasing vs non sig sites 

pdf(file="figures/RichnessThruTime-DECREASING.pdf", height=9, width=10)
par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
for (k in 1:length(sigPos)){
  plot(richness[,sigPos[k]]~timeKYR, 
       type="l", 
       ylim=c(0, max(richness, na.rm=T)), xlim=c(21,0),
       xlab="", ylab="Site Richness")
  mtext(colnames(richness)[sigPos[k]], side=3, line=-1.5, adj=1, cex=0.75)
}
dev.off()

pdf(file="figures/RichnessThruTime-INCREASING.pdf", height=9, width=10)
par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
for (k in 1:length(sigNeg)){
  plot(richness[,sigNeg[k]]~timeKYR, 
       type="l", 
       ylim=c(0, max(richness, na.rm=T)), xlim=c(21,0),
       xlab="", ylab="Site Richness")
  mtext(colnames(richness)[sigNeg[k]], side=3, line=-1.5, adj=1, cex=0.75)
}
dev.off()

pdf(file="figures/RichnessThruTime-nonSig.pdf", height=9, width=10)
par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
for (k in 1:length(nonSig)){
  plot(richness[,nonSig[k]]~timeKYR, 
       type="l", 
       ylim=c(0, max(richness, na.rm=T)),
       xlab="", ylab="Site Richness")
  mtext(colnames(richness)[nonSig[k]], side=3, line=-1.5, adj=1, cex=0.75)
}
dev.off()

# Plot richness trajectory at each site- same plot with sig pos vs sig neg colored differently

pdf(file="figures/RichnessThruTime-all-withPosNegLines.pdf", height=4, width=6)
par(mfrow=c(1,1))
plot(richnessMeans ~ timeKYR, 
     xlim=c(21,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(nonSig)){
  lines(richness[,nonSig[k]]~timeKYR, col="gray")
}
for (k in 1:length(sigPos)){
  lines(richness[,sigPos[k]]~timeKYR, col="blue")
}
for (k in 1:length(sigNeg)){
  lines(richness[,sigNeg[k]]~timeKYR, col="red")
}
lines(richnessMeans~timeKYR, col="black", lwd=2)
dev.off()

# Plot richness trajectory at each site- different plots for sig pos vs sig neg

pdf(file="figures/RichnessThruTime-threePanels-withPosNegLines.pdf", height=4, width=6)
par(mfrow=c(1,3))
plot(richnessMeans ~ timeKYR, 
     xlim=c(21,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(nonSig)){
  lines(richness[,nonSig[k]]~timeKYR, col="gray")
}
lines(nonSigRichnessMeans~timeKYR, col="black", lwd=2)

coolPalette<- colorRampPalette(c("green", "blue"))(43)
plot(richnessMeans ~ timeKYR, 
     xlim=c(21,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(sigPos)){
  lines(richness[,sigPos[k]]~timeKYR, col=coolPalette[max(which(!is.na(richness[,sigPos[k]])))-1])
}
lines(posRichnessMeans~timeKYR, col="black", lwd=2)

warmPalette<- colorRampPalette(c("red", "orange"))(43)
plot(richnessMeans ~ timeKYR, 
     xlim=c(21,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(sigNeg)){
  lines(richness[,sigNeg[k]]~timeKYR, col=warmPalette[max(which(!is.na(richness[,sigNeg[k]])))-1])
}
lines(negRichnessMeans~timeKYR, col="black", lwd=2)

dev.off()

#### Plot proportional richness change- overall and at individual sites ####

# Plot each individual line plus mean proportional change

pdf(file="figures/PropRichnessThruTime-all-withLines.pdf", height=4, width=6)
par(mfrow=c(1,1))
plot(propRichMeans ~ timeStep, 
     xlim=c(21000,0), ylim=c(min(propRichChanges, na.rm=T), max(propRichChanges, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Proportional Richness Change")
for (i in 1:ncol(propRichChanges)){
  lines(propRichChanges[,i]~timeStep, col=rainbow(nrow(propRichChanges))[i]) 
}
lines(propRichMeans~timeStep, col="black", lwd=2)
dev.off()

# Plot only those sites with >25% change

pdf(file="figures/PropRichnessThruTime-pos and neg.pdf", height=4, width=6)
par(mfrow=c(2,1))
plot(noPosSites~timeStep,
     xlim=c(21000,0), ylim=c(0, max(noOverallSites)), 
     pch=16, col="red", 
     xlab="Time step starting (yr BP)", ylab="Number of sites with >25% change")
points(noNegSites~timeStep,
       xlim=c(21000,0), ylim=c(0, max(noPosSites, noNegSites)), 
       pch=16, col="blue")
points(noOverallSites~timeStep, pch=16, col="gray")

plot(propPosSites~timeStep,
     xlim=c(21000,0), ylim=c(0, max(propPosSites, propNegSites)), 
     type="l", col="red", 
     xlab="Time step starting (yr BP)", ylab="Proportion of sites with >25% change")
points(propNegSites~timeStep,
       xlim=c(21000,0), ylim=c(0, max(propPosSites, propNegSites)), 
       type="l", col="blue")

# Is there a relationship between times with lots of positive and negative change?
cor.test(propNegSites, propPosSites)
plot(propNegSites~propPosSites, pch=16)
abline(lm(propNegSites~propPosSites))

# plot the density of different rates of change with time periods
# note: density plot isn't working, or at least, density through time is 0.

proportionDF<- matrix(ncol=2, nrow=NULL)
for (i in 1:ncol(propRichChanges)){
  temp<- as.data.frame(cbind(as.numeric(rownames(propRichChanges)), propRichChanges[,i]))
  proportionDF <- rbind(proportionDF, temp)
}
colnames(proportionDF)<- c("TimeStep", "propChange")   
proportionDF<- proportionDF[-which(is.na(proportionDF$propChange)), ]  
plot(propChange~TimeStep, data=proportionDF, pch=16)





