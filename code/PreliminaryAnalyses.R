#### Preliminary Code ####

# pollenDir<- "/Volumes/bloisgroup/bloislab/Data/Biological" #fix this to connect to server
pollenDir<- "~/Dropbox/Research/Community Paleomodels/projects/pollen/output/all data by site/"

files<- list.files(path=pollenDir, pattern="rel.genus.sum")
sites<- gsub(".rel.genus.sum.csv", "", files)
