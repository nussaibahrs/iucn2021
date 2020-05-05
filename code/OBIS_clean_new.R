# original download from obis for scleractinia
load("data/original/2019-11-13_obis_scleractinia.RData")

library(sp) # for world map
library(rgdal) # also   
library(rgeos)
library(raster)


# Load world map
w <- chronosphere::fetch("NaturalEarth", "land")

# taxonomic cleaning -----------------------------------------------------

scler <- obis.corals
nrow(scler)

#rename columns
colnames(scler)[grep("decimal", colnames(scler))] <- c("latitude", "longitude") #coordinates
scler$tname <- scler$scientificName #tname = taxon name
coral <- subset(scler, !is.na(latitude)==F | !is.na(longitude))

# Isolate genus name
gen <- sub(" .*", "", coral$tname)
gen <- table(gen)

# Check Maeandrina, Montastraea, Lophohelia, Tubastrea, Tubinaria
coral$tname <- sub("Maeandrina", "Meandrina", coral$tname)
coral$tname <- sub("Montastraea", "Montastrea", coral$tname)
coral$tname <- sub("Lophohelia", "Lophelia", coral$tname)
coral$tname <- sub("Caulastraea", "Caulastrea", coral$tname)
coral$tname <- sub("Tubinaria", "Turbinaria", coral$tname)

# remove white spaces
coral$tname <- trimws(coral$tname)

# # Revisions based on Wallace (1999) and Veron (2000) and WoRMS (2010) #not run
# syno <- read.csv("c:/daten/PBDB/coral-synonyms.csv", sep=";", header=T) # carefully checked until 10 occurrences
# 
# 
# # Takes 1 minute
# for (i in 1:nrow(syno)) coral$tname <- sub(syno$OBIS.name[i], syno$valid.name[i], coral$tname, fixed=T)


# Get information on symbiotic status and clean taxonomy
cor.sp <- readxl::read_excel("data/coral_species_new.xlsx") # WK this step could be ignored - follow CTD throughout

# update symbio status
ctdb <- read.csv("data/ctdb_resolved.csv", stringsAsFactors = FALSE)
symbio <- unique(subset(ctdb, trait_name=="Zooxanthellate", select = c("valid", "value")))

#check for duplicates
x <- table(symbio$valid)
dup <- names(x[x > 1])

symbio[symbio$valid %in% dup,]$value <- "both"
symbio <- unique(symbio)

#rename levels
symbio$value <- as.factor(symbio$value)
levels(symbio$value) <- c("az", "ap", "z")

symbio <- symbio[symbio$value == "z",]

cor.sp <- subset(cor.sp, is.na(valid))
cor.sp <- unique(cor.sp)

Species <- cor.sp$scientificName
cor.sp <- cbind(cor.sp, Species)

# Drop taxa not symbiotic and not in ctdb
# remove anything within brackets to be able to match with other dbs
coral$tname <- sub(" \\(.*\\)","", coral$tname) 
coral2 <- merge(coral, symbio, by.x="tname", by.y="valid")

######

# Look for individual species names
View(coral[grep("harttii", coral$tname),])

View(coral2[grep("Mussismilia h", coral2$tname),])
View(coral[grep("Mussismilia h", coral$tname),])
nrow(coral2)

# Consider removing all rounded coordinates (currently not implemented)
# coral2 <- subset(coral2, latitude == round(latitude) & longitude== round(longitude)))

# make bathymetry numeric
coral2$depth <- as.character(coral2$depth)
coral2$depth <- as.numeric(coral2$depth)
# replace greater than 10000 m with NA
coral2$depth[coral2$depth>10000] <- NA

# Solve Problem with bathymetry
# American data might be in feet rather than meters
#    St. Croix, USVI Benthic Composition and Monitoring Data (2002 - Present)
x <- which(coral2$resname==
             "St. Croix, USVI Benthic Composition and Monitoring Data (2002 - Present)")
coral2$depth[x] <- coral2$depth[x]*0.3048       

#    Biogeographic Characterization of Benthic Habitat Communities within the Flower Garden Banks National Marine Sanctuary (2006 - Present)
# Maximum depth 35 m according to website
x <- which(coral2$resname==
             "Biogeographic Characterization of Benthic Habitat Communities within the Flower Garden Banks National Marine Sanctuary (2006 - Present)")
coral2$depth[x] <- coral2$depth[x]*0.3048

#    St. John, USVI Benthic Composition and Monitoring Data (2002 - Present)
# Maximum depth 20 m according to Harborne et al. (2006, Ecology)
x <- which(coral2$resname==
             "St. John, USVI Benthic Composition and Monitoring Data (2002 - Present)")
coral2$depth[x] <- coral2$depth[x]*0.3048

# La Parguera, Puerto Rico Benthic Composition and Monitoring Data (2002 - Present)
# Maximum depth 100 ft according to website
x <- which(coral2$resname==
             "La Parguera, Puerto Rico Benthic Composition and Monitoring Data (2002 - Present)")
coral2$depth[x] <- coral2$depth[x]*0.3048
# Coordinates cleaning ----------------------------------------------------

# round coordinates to 4 digits
coral2$latitude <- round(coral2$latitude, 4)
coral2$longitude <- round(coral2$longitude, 4)  

# exclude high latitude collection with strange data (lat, long mixup?
coral2 <- subset(coral2, latitude!= -55.5333)

# omit deep stuff, currently not implemented
# x <- which(coral2$symbio!="az" & coral2$depth>100)
#   coral3 <- coral2[-x,] 


library(sp)

# delete last remaining errors
obis <- subset(coral2, latitude > -50)

# exclude strange record south of Australia
obis <- subset(obis, latitude!= -40.5292) 

#  exclude occurrences on land
obis <- subset(obis, longitude!= -82.225 & longitude!= -86.000 &
                 longitude!= -79.9989 & longitude!= -49 & longitude!= 134)        

# Shift Madagaskar occurrences from land to sea
obis$longitude[obis$longitude==46.9] <- 44.15

# Shift Kenya occurrences from land to sea
obis$longitude[obis$latitude== -0.0237] <- 40.97
obis$latitude[obis$latitude== -0.0237] <- -2.27
obis$longitude[obis$latitude== -0.0236] <- 41.02
obis$latitude[obis$latitude== -0.0236] <- -2.1

# Delete A. cervicornis from Indian Ocean
del <- which(obis$scientificName=="Acropora cervicornis" & obis$longitude> 10)
obis <- obis[-del,]

########################################################################
# WK This causes issues due to NAs in the datasetName 
# Ditch here. You may consider using the "dataset_id" of these rapid assessments to manually remove them later. 
# However, for the current analyses I don't think it hurts to leave the rapid assessments in

# Omit Virgin Islands National Park Coral Transplant Study with just three species reported
#obis <- subset(obis, datasetName!="Virgin Islands National Park Coral Transplant Study")

#nrow(obis) # occurrences after cleaning
#nrow(unique(subset(obis, select=c("latitude", "longitude"))))
#length(levels(factor(obis$scientificName)))

# filter our rapid assessments
#obis <- subset (obis, datasetName!="CRED Rapid Ecological Assessments of Coral Population in the Pacific Ocean" &
#                  datasetName!="Atlantic and Gulf Rapid Reef Assessment - Benthic" &
#                  datasetName!="CRED Rapid Ecological Assessments of Coral Health and Disease in the Pacific Ocean 2005-2008" &
#                  datasetName!="CRED Rapid Ecological Assessment of Benthic Habitat Cover in the Pacific Ocean 2005-2010" &
#                  datasetName!="CRED Rapid Ecological Assessments of Coral Population in the Pacific Ocean 2007-2010" &
#                  datasetName!="CRED Rapid Ecological Assessment of Invertebrate in the Pacific Ocean")

#x <- subset(obis, datasetName=="St. Croix, USVI Benthic Composition and Monitoring Data (2002 - Present)" | 
#              datasetName=="St. John, USVI Benthic Composition and Monitoring Data (2002 - Present)" |
#              datasetName=="La Parguera, Puerto Rico Benthic Composition and Monitoring Data (2002 - Present)")

# Filter out Carribean monitoring, because the depress rarefied diversity # revised 14.12.2011
#obis <- subset (obis, datasetName!="St. Croix, USVI Benthic Composition and Monitoring Data (2002 - Present)" &
#                  datasetName!="St. John, USVI Benthic Composition and Monitoring Data (2002 - Present)" &
#                  datasetName!="La Parguera, Puerto Rico Benthic Composition and Monitoring Data (2002 - Present)")


#obis <- rbind(obis, x)

nrow(obis) # occurrences after cleaning
nrow(unique(subset(obis, select=c("latitude", "longitude"))))
length(levels(factor(obis$scientificName)))

# check occurrences falling on land 
# WK Your code takes forewer - I even got a sys
# Changed here to my original code
URL <- "http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/physical/ne_110m_ocean.zip"
fil <- basename(URL)
if (!file.exists(fil)) download.file(URL, fil)
fils <- unzip(fil)
oceans <- readOGR(grep("shp$", fils, value=TRUE), "ne_110m_ocean",
                  stringsAsFactors=FALSE, verbose=FALSE)


coordinates(obis) <- ~ longitude + latitude
proj4string(obis) <- CRS(proj4string(oceans))
# Do some cleaning
# Points on land
proper <- over(obis, oceans) # NAs indicate points on land

remove <- which(is.na(proper$scalerank))

oc <- obis[-remove,]   # obis clean
oc <- data.frame(oc) # WK back to ordinary data frame

#w <- chronosphere::fetch("NaturalEarth", "land")[1,] #continents only 
#dat <- obis 
#dat$row_n <- 1:nrow(dat)
#coordinates(dat) <- ~ longitude + latitude
#proj4string(dat) <- proj4string(w)

#land_points <- over(w, dat)

#oc <- obis[-as.numeric(na.omit(land_points$row_n)),]
#oc <- unique(oc)

# Remove non-reef corals if ever remained
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

# load file to remove strange coordinates
rem.oc <- read.csv("obis.cleaning.csv", header=T)
rem.oc <- subset(rem.oc, select=c(round.long, round.lat))
rem.oc <- unique(rem.oc)
rem.oc <- paste(rem.oc$round.long, rem.oc$round.lat)

x2 <- oc$longitude
y2 <- oc$latitude

coord <- paste(x2, y2)  
x <- which(coord %in% rem.oc)

if (length(x) > 0) oc <- oc[-x,]



##############

# download only those in iucn from obis
iucn.cnidaria <- read.csv("data/iucn_resolved.csv", stringsAsFactors = FALSE)
iucn.cnidaria$valid <- sub(" \\(.*\\)", "", iucn.cnidaria$valid)
obis.corals <- subset(oc, scientificName %in% iucn.cnidaria$valid)

length(unique(obis.corals$scientificName))


obis.corals <- subset(obis.corals, select=c("scientificName", "longitude", "latitude"))
colnames(obis.corals)[c(2,3)] <- c("decimalLongitude", "decimalLatitude")
obis.corals <- unique(obis.corals) # make unique at 4 digits

x2 <- obis.corals$decimalLongitude
y2 <- obis.corals$decimalLatitude

# WK Test plot

# Individual points
yl <- c(-90, 90)
op <- par(xaxs="i", yaxs="i", mar=c(5,4.5,2,1))
#  par(mar=c(5,4.5,2,0))
plot(0,0, xlim=c(-180,180), ylim=yl, type="n", xlab="Longitude",
     ylab="Latitude", cex.lab=1.4)
plot(w, col="grey80", add=T, border=NA)
abline(h=0, lty=2, col="grey")
abline(h= -23.433, lty=2, col="grey"); abline(h= 23.433, lty=2, col="grey")
points(x2, y2, pch=19, col="darkred")
#  points(x2, y2, pch=19, col="red", cex=0.5)

box()
par(op)

write.csv(obis.corals, ("data/2019-11-05_obis_scleractinia.csv"), row.names = FALSE)

