library(here)
library(tidyverse)
library(sp)
library(rgeos)
library(rgdal)

#Get OBIS Data
#obis.corals<- occurrence("Scleractinia")
load(here("data", "original", "2019-11-13_obis_scleractinia.RData"))

# Check Maeandrina, Montastraea, Lophohelia, Tubastrea, Tubinaria
obis.corals$scientificName <- sub("Maeandrina", "Meandrina", obis.corals$scientificName)
obis.corals$scientificName <- sub("Montastraea", "Montastrea", obis.corals$scientificName)
obis.corals$scientificName <- sub("Lophohelia", "Lophelia", obis.corals$scientificName)
obis.corals$scientificName <- sub("Caulastraea", "Caulastrea", obis.corals$scientificName)
obis.corals$scientificName <- sub("Tubinaria", "Turbinaria", obis.corals$scientificName)

# Correct misplace spaces
obis.corals$scientificName [obis.corals$scientificName =="Platygyra daedalea "] <- "Platygyra daedalea"

#depth
obis.corals$depth[obis.corals$depth>10000] <- NA

# Solve Problem with bathymetry
# American data might be in feet rather than meters
#    St. Croix, USVI Benthic Composition and Monitoring Data (2002 - Present)
x <- which(obis.corals$datasetName==
             "St. Croix, USVI Benthic Composition and Monitoring Data (2002 - Present)")
obis.corals$depth[x] <- obis.corals$depth[x]*0.3048       

#    Biogeographic Characterization of Benthic Habitat Communities within the Flower Garden Banks National Marine Sanctuary (2006 - Present)
# Maximum depth 35 m according to website
x <- which(obis.corals$datasetName==
             "Biogeographic Characterization of Benthic Habitat Communities within the Flower Garden Banks National Marine Sanctuary (2006 - Present)")
obis.corals$depth[x] <- obis.corals$depth[x]*0.3048

#    St. John, USVI Benthic Composition and Monitoring Data (2002 - Present)
# Maximum depth 20 m according to Harborne et al. (2006, Ecology)
x <- which(obis.corals$datasetName==
             "St. John, USVI Benthic Composition and Monitoring Data (2002 - Present)")
obis.corals$depth[x] <- obis.corals$depth[x]*0.3048

# La Parguera, Puerto Rico Benthic Composition and Monitoring Data (2002 - Present)
# Maximum depth 100 ft according to website
x <- which(obis.corals$datasetName==
             "La Parguera, Puerto Rico Benthic Composition and Monitoring Data (2002 - Present)")
obis.corals$depth[x] <- obis.corals$depth[x]*0.3048

# exclude high latitude collection with strange data (lat, long mixup?
obis.corals <- subset(obis.corals, decimalLatitude!= -55.5333)

# Remove unknown coordinates # they all have coordinates
coral <- subset(obis.corals, is.na(decimalLatitude)==F)

# remove records not identified to species level
coral <- subset(coral, is.na(species)==F)

# delete extinct species
extinct <- c("Sphenotrochus intermedius")
coral <- subset(coral, scientificName%in%extinct==F)

# Get valid species names from OBIS to update coral species list
species <- subset(coral, select=c(family, scientificName, scientificNameAuthorship))
species <- unique(species)

val.sp <- readxl::read_excel(here("data", "coral_species_new.xlsx"))
obis.corals <- merge(coral, val.sp, by="scientificName")

# remove nomina dubia
x <- which(obis.corals$valid=="dubium")
obis.corals <- obis.corals[-x,]

# correct species names (and other information)
x <- which(!is.na(obis.corals$validName))

obis.corals$scientificName[x] <- obis.corals$validName[x]
obis.corals$family.x[x] <- obis.corals$family.y[x]
obis.corals$genus.x[x] <- obis.corals$genus.y[x]

# remove superfluent columns
obis.corals <- subset(obis.corals, select= c(-family.y, -validName, -valid, -scientificNameAuthorship.y, -genus.y))

# delete last remaining error
obis <- subset(obis.corals, decimalLatitude > -50)

# exclude strange record south of Australia
obis <- subset(obis, decimalLatitude!= -40.5292) 
# obis <- subset(obis, latitude!= -40.52917)

#  exclude occurrences on land
obis <- subset(obis, decimalLongitude!= -82.225 & decimalLongitude!= -86.000 &
                 decimalLongitude!= -79.9989 & decimalLongitude!= -49 & decimalLongitude!= 134)        

# Shift Madagaskar occurrences from land to sea
obis$decimalLongitude[obis$decimalLongitude==46.9] <- 44.15

# Shift Kenya occurrences from land to sea
# obis$longitude[obis$latitude== -0.0237] <- 40.97
# obis$latitude[obis$latitude== -0.0237] <- -2.27
obis$decimalLongitude[obis$decimalLatitude== -0.0236] <- 41.02
obis$decimalLatitude[obis$decimalLatitude== -0.0236] <- -2.1

# Delete A. cervicornis from Indian Ocean
del <- which(obis$scientificName=="Acropora cervicornis" & obis$decimalLongitude> 10)
obis <- obis[-del,]


# Omit Virgin Islands National Park Coral Transplant Study with just three species reported
obis <- subset(obis, datasetName!="Virgin Islands National Park Coral Transplant Study")

nrow(obis) # occurrences after cleaning
nrow(unique(subset(obis, select=c("decimalLatitude", "decimalLongitude"))))
length(levels(factor(obis$scientificName)))

# filter our rapid assessments
obis <- subset (obis, datasetName!="CRED Rapid Ecological Assessments of Coral Population in the Pacific Ocean" &
                  datasetName!="Atlantic and Gulf Rapid Reef Assessment - Benthic" &
                  datasetName!="CRED Rapid Ecological Assessments of Coral Health and Disease in the Pacific Ocean 2005-2008" &
                  datasetName!="CRED Rapid Ecological Assessment of Benthic Habitat Cover in the Pacific Ocean 2005-2010" &
                  datasetName!="CRED Rapid Ecological Assessments of Coral Population in the Pacific Ocean 2007-2010" &
                  datasetName!="CRED Rapid Ecological Assessment of Invertebrate in the Pacific Ocean")

x <- subset(obis, datasetName=="St. Croix, USVI Benthic Composition and Monitoring Data (2002 - Present)" | 
              datasetName=="St. John, USVI Benthic Composition and Monitoring Data (2002 - Present)" |
              datasetName=="La Parguera, Puerto Rico Benthic Composition and Monitoring Data (2002 - Present)")

# Filter out Carribean monitoring, because the depress rarefied diversity # revised 14.12.2011
obis <- subset (obis, datasetName!="St. Croix, USVI Benthic Composition and Monitoring Data (2002 - Present)" &
                  datasetName!="St. John, USVI Benthic Composition and Monitoring Data (2002 - Present)" &
                  datasetName!="La Parguera, Puerto Rico Benthic Composition and Monitoring Data (2002 - Present)")


obis <- rbind(obis, x)

nrow(obis) # occurrences after cleaning
nrow(unique(subset(obis, select=c("decimalLatitude", "decimalLongitude"))))
length(levels(factor(obis$scientificName)))

# check occurrences falling on land 
w <- chronosphere::fetch("NaturalEarth", "land")[1,] #continents only 
dat <- obis %>% mutate(row_n = 1:nrow(.))
coordinates(dat) <- ~ decimalLongitude + decimalLatitude
proj4string(dat) <- proj4string(w)

land_points <- over(w, dat)

oc <- obis[-as.numeric(na.omit(land_points$row_n)),]

oc$decimalLongitude <- round(oc$decimalLongitude, 4)
oc$decimalLatitude <- round(oc$decimalLatitude, 4)  
oc <- unique(oc)

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

# load file to remove strange coordinates
rem.oc <- read.csv("/home/nussaibah/Dropbox/nuss/iucn_traits/obis_clean/obis.cleaning.csv", header=T)
rem.oc <- subset(rem.oc, select=c(round.long, round.lat))
rem.oc <- unique(rem.oc)
rem.oc <- paste(rem.oc$round.long, rem.oc$round.lat)

x2 <- oc$decimalLongitude
y2 <- oc$decimalLatitude

coord <- paste(x2, y2)  
x <- which(coord %in% rem.oc)

if (length(x) > 0) oc <- oc[-x,]

# obis.corals <- oc %>% filter(depth <=150)

# download only those in iucn from obis
iucn.corals <- read.csv(here("data", "iucn_resolved.csv"), stringsAsFactors = FALSE)
obis.corals <- subset(oc, scientificName %in% iucn.corals$valid)

length(unique(obis.corals$scientificName))

write.csv(obis.corals %>% dplyr::select(scientificName, decimalLongitude, decimalLatitude), here("data", "2019-11-05_obis_scleractinia.csv"), row.names = FALSE)
