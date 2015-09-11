#### Plant macrofossils ####

# Tower Lake as a test core

library("neotoma")
library("vegan")

datasets <- get_dataset(9738)

core_data <- get_download(datasets)

  i=1
  macro1chron <- core_data[[i]]$sample.meta
  counts <- core_data[[i]]$counts
  macro1counts <- counts[match(rownames(counts), macro1chron$sample.id), ]
  macro1all <- cbind(macro1chron, macro1counts)
  
  i=2
  pollen1chron <- core_data[[i]]$sample.meta
  counts <- core_data[[i]]$counts
  pollen1counts <- counts[match(rownames(counts), pollen1chron$sample.id), ]
  pollen1all <- cbind(pollen1chron, pollen1counts)
  
  i=3
  macro2chron <- core_data[[i]]$sample.meta
  counts <- core_data[[i]]$counts
  macro2counts <- counts[match(rownames(counts), macro2chron$sample.id), ]
  macro2all <- cbind(macro2chron, macro2counts)
  
  i=4
  pollen2chron <- core_data[[i]]$sample.meta
  counts <- core_data[[i]]$counts
  pollen2counts <- counts[match(rownames(counts), pollen2chron$sample.id), ]
  pollen2all <- cbind(pollen2chron, pollen2counts)

  i=5
  macro3chron <- core_data[[i]]$sample.meta
  counts <- core_data[[i]]$counts
  macro3counts <- counts[match(rownames(counts), macro3chron$sample.id), ]
  macro3all <- cbind(macro3chron, macro3counts)
  
  i=6
  pollen3chron <- core_data[[i]]$sample.meta
  counts <- core_data[[i]]$counts
  pollen3counts <- counts[match(rownames(counts), pollen3chron$sample.id), ]
  pollen3all <- cbind(pollen3chron, pollen3counts)
  

 # check datasets
  range(macro1chron$age)
  range(pollen1chron$age)
  range(macro2chron$age)
  range(pollen2chron$age)
  range(macro3chron$age)
  range(pollen3chron$age)
  
#   macro1
#   pollen1
#   macro2
#   pollen2
#   macro3
#   pollen3
#   
#   order is: 2, 3, 1, from youngest to oldest
  
  
# track changes in diversity through time.  keep track of the number of additions, number of extirpations, and standing richness, for each time step
  # From Brown et al. 2001: "We measured colonization by counting the number of
  # taxa present in one time step that were not present in the previous time
  # step. We determined extinction by counting the number of taxa that were
  # present in one time step that were absent in the subsequent time step." 
  # --> This will be referred to in this code as actual colonization and extinctions
  # --> macro1div$adds, macro1div$exts
  
  # "Therefore, we plotted the actual time series of taxonomic richness 
  # --> e.g., macro1div$richness) 
  # as well as the cumulative colonizations (the number of species that would have
  # accumulated if no extinctions had offset the observed colonizations), and the
  # cumulative extinctions (the number of species that would have remained if no
  # colonizations had offset the observed ex- tinctions). In calculating the
  # cumulative colonization and extinction curves, only a single event of
  # colonization or extinction, respectively, was permitted for each species. It
  # frequently happened in the real data sets that a taxon went extinct and
  # subsequently recolonized or vice versa, but in these analyses we counted only
  # the first colonization or the first extinction event for each taxon."
  # --> this will be referred to as cumulative colonizations and extinctions
 
# calculate actual colonization and extinction, as well as standing richness   #### 
macro1div <- matrix(ncol=3, nrow=nrow(macro1all))
colnames(macro1div) <- c("adds", "exts", "richness")
rownames(macro1div)<- macro1all$depth
macro1div <- as.data.frame(macro1div)  

macro1pa <- decostand(macro1counts, method="pa")
macro1div[, 'richness'] <- rowSums(macro1pa)

j <- 1
macro1div[j, 'adds'] <- macro1div[j, 'exts'] <- 0

for (j in 2:nrow(macro1pa)){
  balance <- macro1pa[j-1,] - macro1pa[j,]
  
  #additions (time runs forward)
  if (length(which(balance == 1))>0){
    macro1div[j, 'adds'] <- length(which(balance == 1))
  }else{
    macro1div[j, 'adds'] <- 0
  }
  if (length(which(balance == -1))>0){
    macro1div[j, 'exts'] <- length(which(balance == -1))
  }else{
    macro1div[j, 'exts'] <- 0
  }
}

actualcol <- vector(mode="numeric", length=nrow(macro1div))
actualext <- vector(mode="numeric", length=nrow(macro1div))

startingpoint <- macro1div$richness[length(macro1div$richness)]
actualcol[nrow(macro1div)]<- startingpoint
actualext[nrow(macro1div)]<- startingpoint

for (i in (nrow(macro1div)-1):1){
  actualcol[i] <- actualcol[i+1] + macro1div[i+1,'adds']
  actualext[i] <- actualext[i+1] - macro1div[i+1,'exts']
}

# plot actual numbers
plot(macro1div$richness ~ macro1chron$age, xlim=c(10000, 0), ylim=c(-190, 190), pch=16, type="b", col="black")
points(actualcol ~ macro1chron$age, pch=15, type="b", col="red")
points(actualext ~ macro1chron$age, pch=17, type="b", col="blue")

# calculate cumulative colonizations and extinctions  ####
macro1pacumul <- matrix(data=0, nrow=nrow(macro1pa), ncol=ncol(macro1pa))

# +1 = extinctions (e.g., 1 (present) - 0 (absent))
# -1 = colonization (e.g., 0 (absent) - 1 (present))

# calculate col and ext for whole matrix, then take first and last occurrence of +1 and -1
macro1cum_temp <- macro1pa[2:nrow(macro1pa),] - macro1pa[1:nrow(macro1pa)-1,]

for (j in 1:ncol(macro1cum_temp)){
  exts<- which(macro1cum_temp[,j]==1)
  colons<- which(macro1cum_temp[,j]==-1)
  
  firstext<- max(exts)
  firstcol<- max(colons) 
  
  macro1pacumul[firstcol,j] <- -1
  macro1pacumul[firstext,j] <- 1
}


macro1cumul <- matrix(data=0, ncol=2, nrow=nrow(macro1all))
colnames(macro1cumul) <- c("cumuladds", "cumulexts")
macro1cumul<- as.data.frame(macro1cumul)

for (j in 1:nrow(macro1pacumul)){
  macro1cumul[j,"cumulexts"] <- length(which(macro1pacumul[j,] == 1))
  macro1cumul[j,"cumuladds"] <- length(which(macro1pacumul[j,] == -1))
}

colonizations <- vector(mode="numeric", length=nrow(macro1cumul))
extinctions <- vector(mode="numeric", length=nrow(macro1cumul))

startingpoint <- macro1div$richness[length(macro1div$richness)]
colonizations[nrow(macro1cumul)]<- startingpoint
extinctions[nrow(macro1cumul)]<- startingpoint

for (i in (nrow(macro1cumul)-1):1){
  colonizations[i] <- colonizations[i+1] + macro1cumul[i,'cumuladds']
  extinctions[i] <- extinctions[i+1] - macro1cumul[i,'cumulexts']
}

# plot cumulative numbers
plot(macro1div$richness ~ macro1chron$age, xlim=c(10000, 0), ylim=c(-15, 30), pch=16, type="b", col="black")
points(colonizations ~ macro1chron$age, pch=15, type="b", col="red")
points(extinctions ~ macro1chron$age, pch=17, type="b", col="blue")



par(mfrow=c(2, 1))
# plot actual numbers
plot(macro1div$richness ~ macro1chron$age, xlim=c(10000, 0), ylim=c(-190, 190), 
     pch=16, type="b", col="black", main="Actual col and ext")
points(actualcol ~ macro1chron$age, pch=15, type="b", col="red")
points(actualext ~ macro1chron$age, pch=17, type="b", col="blue")

# plot cumulative numbers
plot(macro1div$richness ~ macro1chron$age, xlim=c(10000, 0), ylim=c(-15, 30), 
     pch=16, type="b", col="black", main="Cumulative col and ext")
points(colonizations ~ macro1chron$age, pch=15, type="b", col="red")
points(extinctions ~ macro1chron$age, pch=17, type="b", col="blue")

#### JESSICA NOTES: START HERE #### OK, this doesn't seem to add up. Perhaps for
#the cumulative plot, I need to re-do standing richness, with the assumptions 
#stated (e.g., that only first col and first ext matter?)


