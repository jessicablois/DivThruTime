#### Spatial analyses ####
pollenDir<- "~/Dropbox/Research/Community Paleomodels/projects/pollen/output/all data by site/"
files<- list.files(path=pollenDir, pattern="gdm.data")
sites<- gsub(".gdm.data.csv", "", files)
times<- seq(0, 21000, by=500)

#### Find names of unique species ####
uniqueSpp<- vector(length=0)
for (i in 1:length(sites)){
  site<- sites[i]
  sitePath<- files[match(site, gsub(".gdm.data.csv", "", files))]
  dat<- read.csv(paste(pollenDir, sitePath, sep="")) #read data
  uniqueSpp<- unique(c(uniqueSpp, colnames(dat[,-1])))
}
uniqueSpp<- uniqueSpp[order(uniqueSpp)]

#### Prepare layers of each species at each time ####
library(sp)
siteMeta<-read.delim("~/Dropbox/Research/Community Paleomodels/projects/pollen/data/All.site.data-withagemodel-finalv2.txt", sep="\t", header=T) 
siteLocs<- siteMeta[match(sites, siteMeta$Handle), c('Handle', 'lon_alb_km', "lat_alb_km")]

for (j in 1:length(times)){
  datByTime<- read.csv(file=paste("~/Dropbox/Research/Community Paleomodels/projects/pollen/output/all data by time/PollenAbund_", times[i], "bp.csv", sep=""), header=T, row.names=1)
  datTrev<- as.matrix(datByTime)
  datTrev[which(datTrev>0)]<- 1

  for (i in length(colnames(datTrev))){
   species<- colnames(datTrev)[i] 
   sitesWithSpp<- names(which(datTrev[,i]>0))
   tempLocs<- siteLocs[match(sitesWithSpp, siteLocs$Handle),]
   tempLocs<- SpatialPoints(tempLocs[,2:3])
   
   temp2<- siteLocs[match(sitesWithSpp, siteLocs$Handle),]
   coordinates(temp2)<- ~lon_alb_km+lat_alb_km
   tempCentroid<- gCentroid(temp2)
   plot(temp2)
   points(tempCentroid, col="red")
  }
}  
  
speciesBySite<- matrix(nrow=length(uniqueSpp), ncol=length(sites))
rownames(speciesBySite)<- uniqueSpp
colnames(speciesBySite)<- sites

for (i in 1:length(uniqueSpp)){
  
  site<- sites[i]
  sitePath<- files[match(site, gsub(".gdm.data.csv", "", files))]
  dat<- read.csv(paste(pollenDir, sitePath, sep="")) #read data
  colIndex<- match(dat$sample.ages, colnames(speciesByTime))

  datPA<- as.matrix(dat[,-1])
  rowIndex<- match(colnames(datPA), rownames(speciesByTime))
  
  datPA[which(datPA>0)]<- 1
  
  dim(speciesByTime[rowIndex, colIndex])
  }
}  
  
site<- sites[1]
sitePath<- files[match(site, gsub(".gdm.data.csv", "", files))]
dat<- read.csv(paste(pollenDir, sitePath, sep="")) #read data

datTimes<- dat[,1] # pull out time periods
datRev<- as.matrix(dat[,-1]) # remove times and convert to matrix
datRev[which(datRev>0)]<- 1 #replace abundance with presence-absence
datRichness<- rowSums(datRev) #calculate genus richness for each time (row)

temp<- cbind(datTimes, datRichness)
