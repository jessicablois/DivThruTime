
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
plot(spatialNegatives, col="red", add=T)
points(negativeCentroid, col="darkred", pch=16, cex=1.5)

#these are sites with increases in diversity through time
spatialPositives<- latlongs[match(colnames(richness)[sigPos], latlongs$Handle),]
coordinates(spatialPositives)<- ~Longitude+Latitude
positiveCentroid<- gCentroid(spatialPositives)
plot(spatialPositives, col="blue", add=T)
points(positiveCentroid, col="darkblue", pch=16, cex=1.5)


warmPalette<- colorRampPalette(c("red", "orange"))(length(timeNeg))
negCols<- apply(richness[,sigNeg], 2, function(x) max(which(!is.na(x))))
#reds= time series starts younger, orange=time series starts older

coolPalette<- colorRampPalette(c("green", "blue"))(length(timeNeg))
posCols<- apply(richness[,sigPos], 2, function(x) max(which(!is.na(x))))
#greens= time series starts younger, blues=time series starts older


data(wrld_simpl)
pdf(file="figures/Map-RichnessThruTime.pdf", height=9, width=10)
par(mfrow=c(1,1))
plot(wrld_simpl[grep("United States", wrld_simpl@data$NAME),], xlim=c(-160, -60), ylim=c(15, 80), axes=FALSE, col='light yellow')
plot(wrld_simpl[grep("Canada", wrld_simpl@data$NAME),], axes=FALSE, col='light yellow', add=T)
plot(wrld_simpl[grep("Mexico", wrld_simpl@data$NAME),], axes=FALSE, col='light yellow', add=T)

plot(spatialNeutrals, pch=3, cex=0.75, col="lightgray", add=T)
plot(spatialNegatives, pch=15, col=warmPalette[negCols], add=T)
plot(spatialPositives, pch=17, col=coolPalette[posCols], add=T)

legend("bottom", pch=c(rep(15, 8), 3, rep(17, 8)), col=c(warmPalette[seq(43,1, -6)], "gray", (coolPalette[seq(1,43, 6)])), 
       legend=timeNeg[c(seq(43,1, -6), 0, seq(1,43, 6))], cex=0.75, horiz=T)

dev.off()


#### plot spatial patterns of richness change ####

par(mfrow=c(5,5))
for (j in 1:nrow(siteRichChanges)){
  specificLocs<- siteLocs[match(colnames(siteRichChanges), sites)]
  plot(siteRichChanges[j,]~specificLocs@coords[,2], pch=16) #plot richness change as a function of latitude
  summary(lm(siteRichChanges[j,]~specificLocs@coords[,2]))  
}


#### Regional analyses ####
# superimpose present-day biome map to delimit regions, then examine richness changes within them
# determine whether same set of species are involved in majority of colonizations or extirpations
# Need to be in Albers
library(maptools)

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

pdf(file="figures/RichnessByBiome.pdf", height=4, width=6)
  plotCols<- biomeCols[match(biome.list.sites, biome.list.full)]
  biomeSampSize <- vector(length=length(biome.list.sites))
  for (i in 1:length(biome.list.sites)-1){
    biomeSampSize[i] <- ncol(richness[,which(biomesinSamp==biome.list.sites[i])])
  }
  biomeSampSize[7]<- 1
  
  par(mar=c(4,4,4,4)+0.1)
  plot(biomeMeanRichness[,1]~timeNeg, type="n", xlim=c(-21000,0), ylim=c(0, max(biomeMeanRichness, na.rm=T)),
       xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
  for (i in 1:length(biome.list.sites)){
    lines(biomeMeanRichness[,i]~timeNeg, pch=16, col=plotCols[i])
  }
  lines(richnessMeans ~ timeNeg, lwd=1.5, col="black")
  legend("bottomleft", bty="n", paste("biome=",biome.list.sites, " (", biomeSampSize, " sites)", sep=""), col=plotCols, lwd=1, cex=0.75)
dev.off()
