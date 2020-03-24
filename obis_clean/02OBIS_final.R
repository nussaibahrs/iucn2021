obis <- read.csv("/home/nussaibah/Dropbox/nuss/iucn_traits/obis_clean/obis.cleaning.csv")

# filter for z corals 
obis <- subset(obis, symbio!="az")

# Exclude ap corals if wanted
obis <- subset(obis, symbio!= "ap")

View(obis[grep("Mussismilia h", obis$tname),])

# round coordinates to 4 digits
obis$latitude <- round(obis$latitude, 4)
obis$longitude <- round(obis$longitude, 4)  

# delete last remaining error
obis <- subset(obis, latitude > -50)

# exclude strange record south of Australia
obis <- subset(obis, latitude!= -40.5292) 
# obis <- subset(obis, latitude!= -40.52917)

#  exclude occurrences on land
obis <- subset(obis, longitude!= -82.225 & longitude!= -86.000 &
                 longitude!= -79.9989 & longitude!= -49 & longitude!= 134)        

# Shift Madagaskar occurrences from land to sea
obis$longitude[obis$longitude==46.9] <- 44.15

# Shift Kenya occurrences from land to sea
# obis$longitude[obis$latitude== -0.0237] <- 40.97
# obis$latitude[obis$latitude== -0.0237] <- -2.27
obis$longitude[obis$latitude== -0.0236] <- 41.02
obis$latitude[obis$latitude== -0.0236] <- -2.1

# Delete A. cervicornis from Indian Ocean
del <- which(obis$tname=="Acropora cervicornis" & obis$longitude> 10)
obis <- obis[-del,]



# Omit Virgin Islands National Park Coral Transplant Study with just three species reported
obis <- subset(obis, resname!="Virgin Islands National Park Coral Transplant Study")

# Export without removing rapid assessments
#  write.table(obis, file="C:/Daten/GBIF/OBIS_analyse2.csv", row.names=F, sep=";")



nrow(obis) # occurrences after cleaning
nrow(unique(subset(obis, select=c("latitude", "longitude"))))
length(levels(factor(obis$tname)))




# filter our rapid assessments
obis <- subset (obis, resname!="CRED Rapid Ecological Assessments of Coral Population in the Pacific Ocean" &
                  resname!="Atlantic and Gulf Rapid Reef Assessment - Benthic" &
                  resname!="CRED Rapid Ecological Assessments of Coral Health and Disease in the Pacific Ocean 2005-2008" &
                  resname!="CRED Rapid Ecological Assessment of Benthic Habitat Cover in the Pacific Ocean 2005-2010" &
                  resname!="CRED Rapid Ecological Assessments of Coral Population in the Pacific Ocean 2007-2010" &
                  resname!="CRED Rapid Ecological Assessment of Invertebrate in the Pacific Ocean")

hist(obis$latitude)

x <- subset(obis, resname=="St. Croix, USVI Benthic Composition and Monitoring Data (2002 - Present)" | 
              resname=="St. John, USVI Benthic Composition and Monitoring Data (2002 - Present)" |
              resname=="La Parguera, Puerto Rico Benthic Composition and Monitoring Data (2002 - Present)")
x$latitude <- round(x$latitude, 1)
x$longitude <- round(x$longitude, 1)
nrow(unique(subset(x, select=c("latitude", "longitude", "tname"))))

x <- subset(x, duplicated(subset(x, select=c("latitude", "longitude", "tname")))==F)


# Filter out Carribean monitoring, because the depress rarefied diversity # revised 14.12.2011
obis <- subset (obis, resname!="St. Croix, USVI Benthic Composition and Monitoring Data (2002 - Present)" &
                  resname!="St. John, USVI Benthic Composition and Monitoring Data (2002 - Present)" &
                  resname!="La Parguera, Puerto Rico Benthic Composition and Monitoring Data (2002 - Present)")


obis <- rbind(obis, x)

nrow(obis) # occurrences after cleaning
nrow(unique(subset(obis, select=c("latitude", "longitude"))))
length(levels(factor(obis$tname)))

hist(obis$latitude)

# Plot a control map
x1 <- obis$longitude
y1 <- obis$latitude
source("maponly.R")


# Create separate genus and species vectors
x <- strsplit(as.character(obis$tname), " ", fixed=T) # output is list
gen <- character()
sp <- character()

# This takes time
for (i in 1:nrow(obis))   gen[i] <- x[[i]][1]
for (i in 1:nrow(obis))   sp[i] <- x[[i]][2]




# Export this clean file
write.table(obis, file="C:/Daten/GBIF/OBIS_analyse.csv", row.names=F, sep=";")
