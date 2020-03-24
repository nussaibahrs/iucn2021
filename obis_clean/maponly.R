library(sp) # for world map
library(rgdal) # also   
library(rgeos)
library(raster)

load("/home/nussaibah/Dropbox/nuss/iucn_traits/obis_clean/obis_cor.RData")
obis.corals <- obis.cleaned
x2 <- obis.corals$decimalLongitude
y2 <- obis.corals$decimalLatitude

r <- raster() 
ras <- rasterize(cbind(x2, y2), r)

test.points <- data.frame(lon=x2, lat=y2)

URL <- "http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/physical/ne_110m_ocean.zip"
fil <- basename(URL)
if (!file.exists(fil)) download.file(URL, fil)
fils <- unzip(fil)
oceans <- readOGR(grep("shp$", fils, value=TRUE), "ne_110m_ocean",
                  stringsAsFactors=FALSE, verbose=FALSE)


coordinates(test.points) <- ~lon+lat
proj4string(test.points) <- CRS(proj4string(oceans))
# Do some cleaning
# Points on land
proper <- over(test.points, oceans) # NAs indicate points on land

remove <- which(is.na(proper$scalerank))


oc <- obis.corals[-remove,]   # obis clean
# round coordinates to 4 digits and make unique
oc$decimalLongitude <- round(oc$decimalLongitude, 4)
oc$decimalLatitude <- round(oc$decimalLatitude, 4)  
oc <- unique(oc)

test <- subset(oc, decimalLatitude < -33)  
test <- subset(oc, decimalLatitude < -22 & decimalLongitude< -60 &  decimalLongitude > -120)  
test <- subset(oc, decimalLatitude > 35) 
test <- subset(oc, decimalLatitude > 32 & decimalLongitude < 0 &  decimalLongitude > -30) 
test <- subset(oc, decimalLatitude > 32 & decimalLongitude < -30 &  decimalLongitude > -80) 
test <- subset(oc, decimalLatitude > 27 & decimalLongitude < -100 &  decimalLongitude > -130) 

# Remove non-reef corals
oc <- subset(oc, scientificName!="Astrangia poculata")
oc <- subset(oc, scientificName!="Ceratotrochus magnaghii")
oc <- subset(oc, scientificName!="Cladopsammia rolandi")
oc <- subset(oc, scientificName!="Crispatotrochus galapagensis")
oc <- subset(oc, scientificName!="Heterocyathus aequicostatus")
oc <- subset(oc, scientificName!="Madracis pharensis")
oc <- subset(oc, scientificName!="Madracis asperula")  
oc <- subset(oc, scientificName!="Oculina diffusa")
oc <- subset(oc, scientificName!="Oculina varicosa")
oc <- subset(oc, scientificName!="Polycyathus isabela")

x2 <- oc$decimalLongitude
y2 <- oc$decimalLatitude

# load file to remove strange coordinates
rem.oc <- read.csv("/home/nussaibah/Dropbox/nuss/iucn_traits/obis_clean/obis.cleaning.csv", header=T)
rem.oc <- subset(rem.oc, select=c(round.long, round.lat))
rem.oc <- unique(rem.oc)
rem.oc <- paste(rem.oc$round.long, rem.oc$round.lat)
coord <- paste(x2, y2)  
x <- which(coord %in% rem.oc)

oc <- oc[-x,]

# Get values to plot
x2 <- oc$decimalLongitude
y2 <- oc$decimalLatitude

save(obis.cleaned, file="obis_cor.RData")
