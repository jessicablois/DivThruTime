#### Diversity Thru Time Functions ####

#### calcSiteRichness ####
calcSiteRichness<- function(dat, minTime, maxTime, pollenThreshold, interval){
  interval<- 1000
  allTimes<- seq(0, 21000, by=interval)
  siteTimes<- seq(minTime, maxTime, by=interval)
  
  datRev<- as.matrix(dat[,-1]) # remove times and convert to matrix
  matchOrder<- match(colnames(datRev), names(pollenThreshold))
  
  for (k in 1:nrow(datRev)){
    
    datRev[k,which(datRev[k,] >= pollenThreshold[matchOrder])]<- 1
    datRev[k,which(datRev[k,] < pollenThreshold[matchOrder])]<- 0
    
  }
  
  datRichness<- rowSums(datRev) #calculate genus richness for each time (row)
  datRichness<- datRichness[match(siteTimes, dat[,1])]
  
  fullRich<- as.vector(rep(NaN, length(allTimes)), mode="numeric")
  fullRich[match(siteTimes, allTimes)]<- datRichness  
  
  return(fullRich)
}


#### sampleSize ####
sampleSize<- function(x) length(which(!is.na(x)))




#### OLD: calcSiteRichness ####
# calcSiteRichness<- function(site, minTime, maxTime){
#   sitePath<- files[match(sites[i], gsub(".gdm.data.csv", "", files))]
#   times<- seq(minTime, maxTime, by=500)
#   
#   dat<- read.csv(paste(pollenDir, sitePath, sep="")) #read data
#   datTimes<- dat[,1] # pull out time periods
#   
#   datRev<- as.matrix(dat[,-1]) # remove times and convert to matrix
#   datRev[which(datRev>0)]<- 1 #replace abundance with presence-absence
#   datRichness<- rowSums(datRev) #calculate genus richness for each time (row)
#   
#   fullRich<- as.vector(rep(NaN, length(times)), mode="numeric")
#   fullRich[match(datTimes, times)]<- datRichness  
#   return(fullRich)
# }

