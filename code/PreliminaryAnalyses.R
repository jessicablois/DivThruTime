#### Preliminary Code ####
source("code/DivThruTimeFunctions.R")

# Note also that this is calling the original data
# Do I want to use the pollen thresholds established by Kaitlin and Diego?

# pollenDir<- "/Volumes/bloisgroup/bloislab/Data/Biological" #fix this later to connect to server instead of local computer
pollenDir<- "~/Dropbox/Research/Community Paleomodels/projects/pollen/output/all data by site/"
files<- list.files(path=pollenDir, pattern="gdm.data")
sites<- gsub(".gdm.data.csv", "", files)
times<- seq(0, 21000, by=500)

richness<- matrix(ncol=length(times), nrow=length(sites))
colnames(richness)<- times
rownames(richness)<- sites

for (i in 1:length(sites)){
  richness[i,]<- calcSiteRichness(sites[i], 0, 21000)
}

# Calculate sample size
sampSize<- apply(richness, 2, sampleSize)

# Plot the mean genus richness across all sites and overlay sample size
pdf(file="figures/RichnessSampSizeThruTime-all.pdf", height=4, width=6)
  par(mar=c(4,4,4,4)+0.1)
  plot(colMeans(richness, na.rm=T)~times, type="l", xlim=c(0, 21000), xlab="Time slice (kyr BP)", ylab="Mean Genus Richness")
  par(new=T)
  plot(sampSize~times, pch=16, col="red", ylim=c(0, max(sampSize)), xlim=c(0, 21000), axes=F, xlab="", ylab="")
  axis(4, at=seq(0, max(sampSize), by=100))
  mtext("Number of Sites", 4, line=3, col="red")
dev.off()