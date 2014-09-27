#### Diversity Thru Time Functions ####

#### calcSiteRichness ####
calcSiteRichness<- function(site, minTime, maxTime){
  sitePath<- files[match(sites[i], gsub(".gdm.data.csv", "", files))]
  times<- seq(minTime, maxTime, by=500)
  
  dat<- read.csv(paste(pollenDir, sitePath, sep="")) #read data
  datTimes<- dat[,1] # pull out time periods
  
  datRev<- as.matrix(dat[,-1]) # remove times and convert to matrix
  datRev[which(datRev>0)]<- 1 #replace abundance with presence-absence
  datRichness<- rowSums(datRev) #calculate genus richness for each time (row)
  
  fullRich<- as.vector(rep(NaN, length(times)), mode="numeric")
  fullRich[match(datTimes, times)]<- datRichness  
  return(fullRich)
}

#### sampleSize ####
sampleSize<- function(x) length(which(!is.na(x)))