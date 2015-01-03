#### Part3: Calc taxon additions and subtractions ####
source("code/Part0-GlobalParam.R")

changes<- as.data.frame(matrix(ncol=nTax+2, nrow=length(allTimes)-1))
changes[,1]<- allTimes[length(allTimes):2]
changes[,2]<- allTimes[(length(allTimes)-1):1]
colnames(changes)[1:2]<- c("Btime", "Etime")
colnames(changes)[3:ncol(changes)]<- taxa
changesAll<- list(changes)

for (i in 1:length(sites)){
  sitePath<- files[match(sites[i], gsub(".gdm.data.csv", "", files))]
  dat<- read.csv(paste(pollenDir, sitePath, sep="")) #read data
  
  dat<- dat[na.omit(match(allTimes, dat$sample.ages)),]
  datTimes<- dat[,1] # pull out time periods
  
  if (length(datTimes)<=1){
    changesAll[[i]]<- NA     
  }else{ 
    
    #replace abundance with presence-absence 
    datRev<- as.matrix(dat[,-1]) # remove times and convert to matrix
    matchOrder<- match(colnames(datRev), names(pollenThreshold))
 
    for (k in 1:nrow(datRev)){  
      datRev[k,which(datRev[k,] >= pollenThreshold[matchOrder])]<- 1
      datRev[k,which(datRev[k,] < pollenThreshold[matchOrder])]<- 0    
    }
       
    for (j in length(datTimes):2){
      changesTemp<- datRev[j-1,]-datRev[j,]  #0 means no change, -1 means species loss, +1 means species gain
      changesTemp<- changesTemp[-which(is.na(match(names(changesTemp), taxa)))] # remove taxa not in original 20
      changesTemp<- changesTemp[order(match(names(changesTemp), taxa))]
      changes[which(changes[,1]==datTimes[j]), na.omit(match(names(changesTemp), colnames(changes)))]<- changesTemp
    } 
    changesAll[[i]]<-changes 
  }
  
}

#### compare changes with centroids for each taxon and time ####
library(geosphere)

#0 means no change, -1 means species loss, +1 means species gain

load("workspaces/centroids.RData")

bearings<- as.data.frame(matrix(ncol=nTax+2, nrow=length(allTimes)-1))
bearings[,1]<- allTimes[length(allTimes):2]
bearings[,2]<- allTimes[(length(allTimes)-1):1]
colnames(bearings)[1:2]<- c("Btime", "Etime")
colnames(bearings)[3:ncol(bearings)]<- taxa
bearingsAll<- list(bearings)

for (i in 1:length(sites)){
 siteCoords<- siteLocs[match(sites[i], latlongs$Handle)]
   
 changes<- changesAll[[i]]
 if (all(is.na(changes))){
   bearingsAll[[i]]<-NaN     
 }else{
 
 goodRows<- which(!is.na(rowSums(changes[,-c(1:2)])))
 
 for (j in 1:length(goodRows)){
   Btime<- changes[goodRows[j],1] 
   Etime<- changes[goodRows[j],2]
   allChanges<- changes[goodRows[j],3:22]
 
  k<- which(allTimes==Btime)
  centroids<- centroidsAll[[k]]
  centroids[,2]<- as.numeric(centroids[,2])
  centroids[,3]<- as.numeric(centroids[,3])
  if (any(is.na(centroids[,3]))){
  centroids<- centroids[-which(is.na(centroids[,2])),]
  }
  
  taxaCoords<- SpatialPoints(centroids[,2:3], proj4string=albersCRS)
  #taxaCoords<- SpatialPoints(centroids[match(colnames(allChanges), centroids$taxa),2:3], proj4string=albersCRS)

  taxaCoords<- spTransform(taxaCoords, originalCRS)
  siteCoords<- spTransform(siteCoords, originalCRS)
  bearings[goodRows[j], match(centroids$taxa, colnames(bearings))]<- bearing(taxaCoords, siteCoords)
    
 }
 bearingsAll[[i]]<-bearings 
  
 }
}
 
rowMeans(bearingsAll[[1]])

bearingMeans<- matrix(ncol=3, nrow=length(sites))
colnames(bearingMeans)<- c("Additions", "Subtractions", "Neutrals")
for (i in 1:length(sites)){
  changes<- changesAll[[i]]
  bearingsVector<- as.vector(as.matrix(bearingsAll[[i]]))
  
  bearingMeans[i,1]<- mean(bearingsVector[which(changes==1)], na.rm=T)
  bearingMeans[i,2]<- mean(bearingsVector[which(changes==-1)], na.rm=T)
  bearingMeans[i,3]<- mean(bearingsVector[which(changes==0)], na.rm=T)
    
}

pdf(file="figures/rose-diagrams.pdf", 
par(mfrow=c(1,3))
x<- circular(bearingMeans[,1], type="directions", units="degrees", template="geographics") 
rose.diag(x, bins = 18, main = 'Additions')
x<- circular(bearingMeans[,2], type="directions", units="degrees", template="geographics")
rose.diag(x, bins = 18, main = 'Subtractions')
x<- circular(bearingMeans[,3], type="directions", units="degrees", template="geographics")
rose.diag(x, bins = 18, main = 'Neutrals')
dev.off()



 