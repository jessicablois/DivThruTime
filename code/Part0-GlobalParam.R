#### Part0: Set Global Parameters ####

#### Arrange pollen data so that presence-absence corresponds to the 5% threshold ####
dat0<- read.csv("~/Dropbox/Research/Community Paleomodels/projects/pollen/output/all data by time/PollenAbund_Q_0bp.csv", header=T, row.names=1)
pollenMax<- apply(dat0, 2, max)
pollenThreshold<- pollenMax*0.05

#### set times ####
interval<- 1000
allTimes<- seq(0, 21000, by=interval)

#### set sites list ####
pollenDir<- "~/Dropbox/Research/Community Paleomodels/projects/pollen/output/all data by site/"
files<- list.files(path=pollenDir, pattern="gdm.data")
sites<- gsub(".gdm.data.csv", "", files)

#### set taxa ####
nTax<- 20
taxa<- c("Abies", "Acer", "Alnus", "Ambrosia.type", "Artemisia", "Betula", "Carya", 
         "Corylus", "Fagus", "Fraxinus", "Juglans", "Larix", "Ostrya.Carpinus", "Picea",
         "Pinus", "Platanus", "Populus", "Quercus", "Salix", "Ulmus")

#### set site locations ####
library(rgdal)

albersCRS<- CRS("+proj=aea +lat_1=20 +lat2=60 +lat0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
originalCRS<- CRS("+proj=longlat +datum=WGS84")

latlongs<- read.delim("~/Dropbox/Research/Community Paleomodels/projects/pollen/data/All.site.data-withagemodel-finalv2.txt", header=T)
latlongs<- latlongs[match(sites, latlongs$Handle),]
siteLocs<- SpatialPoints(latlongs[,c('Longitude','Latitude')], proj4string = CRS("+proj=longlat +datum=WGS84"))
siteLocs<- spTransform(siteLocs, albersCRS)




