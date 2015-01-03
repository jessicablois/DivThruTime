library(rgeos)

## Calculate centroids of each of the 20 taxa at each time.
pollenDir<- "~/Dropbox/Research/Community Paleomodels/projects/pollen/"

times<- seq(0, 21000, by=1000)

albersCRS<- CRS("+proj=aea +lat_1=20 +lat2=60 +lat0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

latlongs<- read.delim(paste(pollenDir, "data/All.site.data-withagemodel-finalv2.txt", sep=""), header=T)
siteLocs<- SpatialPoints(latlongs[,c('Longitude','Latitude')], proj4string = CRS("+proj=longlat +datum=WGS84"))
siteLocs<- spTransform(siteLocs, albersCRS)

#times<- seq(minTime, maxTime, by=1000)
nTax<- 20
taxa<- c("Abies", "Acer", "Alnus", "Ambrosia.type", "Artemisia", "Betula", "Carya", 
         "Corylus", "Fagus", "Fraxinus", "Juglans", "Larix", "Ostrya.Carpinus", "Picea",
         "Pinus", "Platanus", "Populus", "Quercus", "Salix", "Ulmus")

#set pollen threshold
dat0<- read.csv("~/Dropbox/Research/Community Paleomodels/projects/pollen/output/all data by time/PollenAbund_Q_0bp.csv", header=T, row.names=1)
pollenMax<- apply(dat0, 2, max)
pollenThreshold<- pollenMax*0.05

centroids<- as.data.frame(matrix(nrow=nTax, ncol=3))
colnames(centroids)<- c("taxa", "long", "lat")
centroidsAll<- list(centroids)

for (i in 1:length(times)){
  dat<- read.csv(paste(pollenDir, "output/all data by time/PollenAbund_", times[i], "bp.csv", sep=""))
  datRev<- as.matrix(dat[,-1]) # remove times and convert to matrix
  matchOrder<- match(colnames(datRev), names(pollenThreshold))
 
  #convert to presence-absence
  for (k in 1:nrow(datRev)){
    datRev[k,which(datRev[k,] >= pollenThreshold[matchOrder])]<- 1
    datRev[k,which(datRev[k,] < pollenThreshold[matchOrder])]<- 0   
  }
  
  #calculate which sites a taxon is located at
  sites<- as.character(dat[,1])
  for (j in 1:length(taxa)){
    focal<- taxa[j]
    col<- match(focal, colnames(datRev))
    taxaSites<- sites[which(datRev[,col]==1)]
  
    if (length(taxaSites)==0){
      centroids[j,1]<- taxa[j]
      centroids[j,2:3]<- c("NaN", "NaN")     
    }else{
    tempCentroid<-gCentroid(siteLocs[match(taxaSites, latlongs$Handle)])
    plot(siteLocs[match(taxaSites, latlongs$Handle)])
    plot(tempCentroid, col="red", add=T, pch=16)
    mtext(side=3, taxa[j])
    
    centroids[j,1]<- taxa[j]
    centroids[j,2:3]<- coordinates(tempCentroid)
    }
  }
  
  centroidsAll[[i]]<- centroids  
}

save.image(file="workspaces/centroids.RData")


