
load("workspaces/richness.RData")
source("code/Part0-GlobalParam.R")

library(rgeos)
library(maptools)

#### Spatially plot the positives, negatives, and neutrals ####

#these are sites with no changes in diversity through time
spatialNeutrals<- latlongs[match(colnames(richness)[nonSig], latlongs$Handle),]
coordinates(spatialNeutrals)<- ~Longitude+Latitude
neutralCentroid<- gCentroid(spatialNeutrals)
plot(spatialNeutrals, col="lightgray")
points(neutralCentroid, col="darkgray", pch=16, cex=1.5)

#these are sites with decreases in diversity through time
spatialNegatives<- latlongs[match(colnames(richness)[sigNeg], latlongs$Handle),]
coordinates(spatialNegatives)<- ~Longitude+Latitude
negativeCentroid<- gCentroid(spatialNegatives)
plot(spatialNegatives, col="blue", add=T)
points(negativeCentroid, col="darkblue", pch=16, cex=1.5)

#these are sites with increases in diversity through time
spatialPositives<- latlongs[match(colnames(richness)[sigPos], latlongs$Handle),]
coordinates(spatialPositives)<- ~Longitude+Latitude
positiveCentroid<- gCentroid(spatialPositives)
plot(spatialPositives, col="red", add=T)
points(positiveCentroid, col="darkred", pch=16, cex=1.5)


warmPalette<- colorRampPalette(c("red", "yellow"))(length(timeNeg)-1)
warmPalette<- c("gray", rev(warmPalette))
coolPalette<- colorRampPalette(c("green", "blue"))(length(timeNeg)-1)
coolPalette<- c("gray", coolPalette)
#reds= orange=time series starts younger, time series starts older 
#greens= time series starts younger, blues=time series starts older

negCols<- apply(richness[,sigNeg], 2, function(x) max(which(!is.na(x))))
posCols<- apply(richness[,sigPos], 2, function(x) max(which(!is.na(x))))

data(wrld_simpl)
pdf(file="figures/Map-RichnessThruTime.pdf", height=9, width=10)
par(mfrow=c(1,1))
plot(wrld_simpl[grep("United States", wrld_simpl@data$NAME),], xlim=c(-160, -60), ylim=c(15, 80), axes=FALSE, col='light yellow')
plot(wrld_simpl[grep("Canada", wrld_simpl@data$NAME),], axes=FALSE, col='light yellow', add=T)
plot(wrld_simpl[grep("Mexico", wrld_simpl@data$NAME),], axes=FALSE, col='light yellow', add=T)

plot(spatialNeutrals, pch=3, cex=0.75, col="lightgray", add=T)
plot(spatialNegatives, pch=17, col=coolPalette[negCols], add=T)
plot(spatialPositives, pch=15, col=warmPalette[posCols], add=T)

l <- legend("bottomleft", legend=rep(NA, 16), 
       ncol=2, pch=c(rep(15, 8), rep(17, 8)), bty="n",
       col=c(warmPalette[seq(43,1, -6)], coolPalette[seq(43,1, -6)]))
text(l$text$x+4, l$text$y, c(abs(timeNeg[seq(43,1, -6)]), rep(NA, 8)), pos = 4, cex=0.5)

dev.off()

pdf(file="figures/Map-RichnessThruTime-3panels.pdf", height=9, width=10)
par(mfrow=c(1,3))

#plot neutrals
plot(wrld_simpl[grep("United States", wrld_simpl@data$NAME),], xlim=c(-160, -60), ylim=c(15, 80), axes=FALSE, col='light yellow')
plot(wrld_simpl[grep("Canada", wrld_simpl@data$NAME),], axes=FALSE, col='light yellow', add=T)
plot(wrld_simpl[grep("Mexico", wrld_simpl@data$NAME),], axes=FALSE, col='light yellow', add=T)
plot(spatialNeutrals, pch=3, cex=0.75, col="lightgray", add=T)

l <- legend("bottomleft", legend=rep(NA, 16), 
            ncol=2, pch=c(rep(15, 8), rep(17, 8)), bty="n",
            col=c(warmPalette[seq(43,1, -6)], coolPalette[seq(43,1, -6)]))
text(l$text$x+4, l$text$y, c(abs(timeNeg[seq(43,1, -6)]), rep(NA, 8)), pos = 4, cex=0.5)

#plot negatives
plot(wrld_simpl[grep("United States", wrld_simpl@data$NAME),], xlim=c(-160, -60), ylim=c(15, 80), axes=FALSE, col='light yellow')
plot(wrld_simpl[grep("Canada", wrld_simpl@data$NAME),], axes=FALSE, col='light yellow', add=T)
plot(wrld_simpl[grep("Mexico", wrld_simpl@data$NAME),], axes=FALSE, col='light yellow', add=T)
plot(spatialNegatives, pch=17, col=coolPalette[negCols], add=T)

#plot positives
plot(wrld_simpl[grep("United States", wrld_simpl@data$NAME),], xlim=c(-160, -60), ylim=c(15, 80), axes=FALSE, col='light yellow')
plot(wrld_simpl[grep("Canada", wrld_simpl@data$NAME),], axes=FALSE, col='light yellow', add=T)
plot(wrld_simpl[grep("Mexico", wrld_simpl@data$NAME),], axes=FALSE, col='light yellow', add=T)
plot(spatialPositives, pch=15, col=warmPalette[posCols], add=T)

dev.off()


#### plot spatial patterns of richness change- latitudinal and longitudinal gradients ####

par(mfrow=c(2,4), mar=c(5,4,4,1))
times<- seq(21000, 0, -3000)
for (j in 1:length(times)){
  specificLocs<- siteLocs[match(colnames(richness), sites)]
  plot(richness[which(rownames(richness)==times[j]),]~specificLocs@coords[,2], 
       pch=16, xlab="Latitude", ylab="Site richness") #plot richness change as a function of latitude
  abline(lm(richness[which(rownames(richness)==times[j]),]~specificLocs@coords[,2]))
  legend("topright", legend=times[j], bty="n")
  print(times[j])
  print(summary(lm(richness[which(rownames(richness)==times[j]),]~specificLocs@coords[,2])))  
}

par(mfrow=c(2,4), mar=c(5,4,4,1))
times<- seq(21000, 0, -3000)
for (j in 1:length(times)){
  specificLocs<- siteLocs[match(colnames(richness), sites)]
  plot(richness[which(rownames(richness)==times[j]),]~specificLocs@coords[,1], 
       pch=16, xlab="Longitude", ylab="Site richness") #plot richness change as a function of latitude
  abline(lm(richness[which(rownames(richness)==times[j]),]~specificLocs@coords[,1]))
  legend("topright", legend=times[j], bty="n")
  print(times[j])
  print(summary(lm(richness[which(rownames(richness)==times[j]),]~specificLocs@coords[,1])))  
}


pdf(file="figures/LatLongGradient.pdf", width=15, height=10)
par(mfrow=c(2,4), mar=c(5,4,4,1))
times<- seq(21000, 0, -3000)
richSeq <- seq(min(richness, na.rm=T), max(richness, na.rm=T), 1)
cexOverall<- seq(min(richness, na.rm=T), max(richness, na.rm=T), 1)/12
for (j in 1:length(times)){
  locs<- latlongs[match(colnames(richness), sites),] 
  coordinates(locs)<- ~Longitude+Latitude
  locs<- locs[which(!is.na(richness[which(rownames(richness)==times[j]),])),]
  tempDat<- richness[which(rownames(richness)==times[j]), which(!is.na(richness[which(rownames(richness)==times[j]),]))]
  tempDat<- as.data.frame(tempDat)
  colnames(tempDat)<- "richness"
  tempDF<- SpatialPointsDataFrame(coords=locs@coords, data=tempDat)
  cexPts<- cexOverall[tempDF@data$richness]
  plot(wrld_simpl[grep("United States", wrld_simpl@data$NAME),], 
       xlim=c(-160, -60), ylim=c(15, 80), axes=FALSE, main=times[j])
  plot(wrld_simpl[grep("Canada", wrld_simpl@data$NAME),], axes=FALSE, add=T)
  plot(wrld_simpl[grep("Mexico", wrld_simpl@data$NAME),], axes=FALSE, add=T)
  plot(locs, col="darkgreen", pch=16, cex=cexPts, add=T)
}
temp<- seq(2, 24, 2)
legend("bottom", pch=16, pt.cex=cexOverall[temp], legend=richSeq[temp], horiz=T)
dev.off()


par(mfrow=c(5,5))
for (j in 1:nrow(siteRichChanges)){
  specificLocs<- siteLocs[match(colnames(siteRichChanges), sites)]
  plot(siteRichChanges[j,]~specificLocs@coords[,2], pch=16) #plot richness change as a function of latitude
  summary(lm(siteRichChanges[j,]~specificLocs@coords[,2]))  
}

pdf(file="figures/Map-RichnessThruTime-All.pdf", height=9, width=10)
par(mfrow=c(2,2), mar=c(4,4,4,4)+0.1)
plot(richnessMeans ~ timeNeg, 
     xlim=c(-21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness",
     main="All sites")
for (i in 1:ncol(richness)){
  lines(richness[,i]~timeNeg, col="gray") 
}
lines(richnessMeans~timeNeg, col="black", lwd=2)

plot(richnessMeans ~ timeNeg, 
     xlim=c(-21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness",
     main="Neutral sites")
for (k in 1:length(nonSig)){
  lines(richness[,nonSig[k]]~timeNeg, col="gray")
}
lines(nonSigRichnessMeans~timeNeg, col="black", lwd=2)

plot(richnessMeans ~ timeNeg, 
     xlim=c(-21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness",
     main="Increasing sites")
for (k in 1:length(sigPos)){
  lines(richness[,sigPos[k]]~timeNeg, col=warmPalette[max(which(!is.na(richness[,sigPos[k]])))-1])
}
lines(posRichnessMeans~timeNeg, col="black", lwd=2)

plot(richnessMeans ~ timeNeg, 
     xlim=c(-21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness",
     main="Decreasing sites")
for (k in 1:length(sigNeg)){
  lines(richness[,sigNeg[k]]~timeNeg, col=coolPalette[max(which(!is.na(richness[,sigNeg[k]])))-1])
}
lines(negRichnessMeans~timeNeg, col="black", lwd=2)
dev.off()


#### Regional analyses ####
# superimpose present-day biome map to delimit regions, then examine richness changes within them
# determine whether same set of species are involved in majority of colonizations or extirpations
# Need to be in Albers

#read in biome layer
biomes<- readShapePoly("~/Dropbox/Research-Wisconsin/GIS Data/TerrestrialEcoregions/NorthAmerica_ecoregions.shp", proj4string=originalCRS) 

#reproject everything to AlbersCRS
biomesSp<- spTransform(biomes, albersCRS)

#convert biomes (column 6) to factor for plotting
biomesSp@data[,6]<- as.factor(biomesSp@data[,6])
biome.list.full<- levels(biomesSp$BIOME)
biomeCols<- terrain.colors(length(biome.list.full))

biomeCols[14]<- "gray"  #swap colors to make 4,5,6,8,9 more distinguishable, 98=gray
biomeCols[c(1,5)]<- biomeCols[c(5,1)]  #swap colors to make 4,5,6,8,9 more distinguishable, 98=gray
biomeCols[c(12,4)]<- biomeCols[c(4,12)]  #swap colors to make 4,5,6,8,9 more distinguishable, 98=gray

jpeg(file="figures/PresentDayBiomesAtPoints.jpg", height=960, width=960, units="px")
  par(mfrow=c(1,1))
  plot(biomesSp, col=biomeCols[match(biomesSp@data$BIOME, biome.list.full)])
  points(siteLocs, col="red", pch=16, cex=0.75)
dev.off()


#plot biomes with correct colors
plot(biomesSp, lwd=0.05)
for (i in 1:length(biome.list.full)){
  plot(biomesSp[which(biomesSp$BIOME==biome.list.full[i]),], add=T, col=biomeCols[i], lwd=0.05)
}
legend("topleft", biome.list.full, pch=19, col=biomeCols)
points(siteLocs, col="red", pch=16, cex=0.5)


biomeAtLocs<- over(siteLocs, biomesSp)
biomesinSamp<- as.numeric(as.character(biomeAtLocs$BIOME))
biome.list.sites<- as.numeric(na.omit(unique(biomesinSamp)))

if (restrict == "yes"){  #need to trim out sites with low sampling, the same as in Part 1a.
  biomesinSamp<- biomesinSamp[which(noSamp >= minSamp)] #values should be carried over in richness workspace
  biome.list.sites<- as.numeric(na.omit(unique(biomesinSamp)))
}

biomeMeanRichness <- matrix(ncol=length(biome.list.sites), nrow=nrow(richness))
for (i in 1:length(biome.list.sites)){
  if (!is.vector(richness[,which(biomesinSamp==biome.list.sites[i])])){
    biomeMeanRichness[,i] <- rowMeans(richness[,which(biomesinSamp==biome.list.sites[i])], na.rm=T)
  }
  if (is.vector(richness[,which(biomesinSamp==biome.list.sites[i])])){
    biomeMeanRichness[,i] <- richness[,which(biomesinSamp==biome.list.sites[i])]
  }
}
colnames(biomeMeanRichness) <- biome.list.sites
rownames(biomeMeanRichness) <- rownames(richness)

pdf(file="figures/RichnessByBiome.pdf", height=6, width=10)
  plotCols<- biomeCols[match(biome.list.sites, biome.list.full)]
  biomeSampSize <- vector(length=length(biome.list.sites))
  for (i in 1:length(biome.list.sites)-1){
    if(i==5){biomeSampSize[i]<- 1}else{
    biomeSampSize[i] <- ncol(richness[,which(biomesinSamp==biome.list.sites[i])])
    }
  }
  
  #biomeSampSize[7]<- 1
  
  par(mar=c(4,4,4,4)+0.1)
  plot(biomeMeanRichness[,1]~timeNeg, type="n", xlim=c(-21000,0), ylim=c(0, max(biomeMeanRichness, na.rm=T)),
       xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
  for (i in 1:length(biome.list.sites)){
    lines(biomeMeanRichness[,i]~timeNeg, pch=16, lwd=2, col=plotCols[i])
  }
  lines(richnessMeans ~ timeNeg, lwd=2.5, col="black")
  legend("bottomleft", bty="n", paste("biome=",biome.list.sites, " (", biomeSampSize, " sites)", sep=""), col=plotCols, lwd=2, cex=0.75)
dev.off()


save(list=c("spatialNeutrals", "spatialNegatives", "spatialPositives", 
            "biomesinSamp", "biome.list.sites", "biomeMeanRichness", "biomeSampSize"), 
     file="workspaces/spatial.RData")


#### plot spatial patterns of richness change ####
# this just plots richness as a function of latitude or longitude- not useful

# par(mfrow=c(5,5))
# for (j in 1:nrow(siteRichChanges)){
#   specificLocs<- siteLocs[match(colnames(siteRichChanges), sites)]
#   plot(siteRichChanges[j,]~specificLocs@coords[,2], pch=16) #plot richness change as a function of latitude
#   summary(lm(siteRichChanges[j,]~specificLocs@coords[,2]))  
# }

