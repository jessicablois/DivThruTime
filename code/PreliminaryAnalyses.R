#### Preliminary Code ####

# pollenDir<- "/Volumes/bloisgroup/bloislab/Data/Biological" #fix this to connect to server
pollenDir<- "~/Dropbox/Research/Community Paleomodels/projects/pollen/output/all data by site/"

files<- list.files(path=pollenDir, pattern="gdm.data")
sites<- gsub(".gdm.data.csv", "", files)

times<- seq(0, 21000, by=500)

dat<- read.csv(paste(pollenDir, files[i], sep=""))