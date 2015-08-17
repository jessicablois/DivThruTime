#### Plant macrofossils ####

# Tower Lake as a test core

library("neotoma")
library("vegan")

datasets <- get_dataset(9738)

core_data <- get_download(datasets)

# for (i in 1:length(core_data)){
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
  
  # }
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

# plot
plot(macro1div$richness ~ macro1chron$age, xlim=c(10000, 0), ylim=c(0, 15), type="l")


# but this is not cumulative additions
# JESSICA START HERE NEXT TIME
