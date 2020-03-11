# Use the OBIS R package

require(robis)
library(here)

robis
 records <- occurrence(scientificname = "Scleractinia") # retrieve records directly from OBIS
 save(records, file ="Scleractinia.RData")

 
# load("Scleractinia.RData")
# remove fields without information
load(here("data", "original", "2019-11-13_obis_scleractinia.RData"))
records <- obis.corals[,1:59] # WARNING: needs adujstment so not to lose longitude and taxonomic level (now in a different column)
save(records, file ="Scleractinia.RData")
# load again
load("Scleractinia.RData")


# Remove unknown coordinates # they all have coordinates
coral <- subset(records, is.na(decimalLatitude)==F)

# remove records not identified to species level
coral <- subset(coral, is.na(species)==F)

# delete extinct species
 extinct <- c("Sphenotrochus intermedius")
  coral <- subset(coral, scientificName%in%extinct==F)

# Get valid species names from OBIS to update coral species list
  species <- subset(coral, select=c(family, scientificName, scientificNameAuthorship))
  species <- unique(species)


   # Get information on symbiotic status and clean taxonomy
  #      cor.sp <- read.table(file="C:/Daten/PBDB/coral-species.csv", header=TRUE, sep=";", quote="\"")
  #      cor.sp$gensp <- paste(cor.sp$genus, cor.sp$species)
  #      cor.sp <- subset(cor.sp, select=c(symbio, genus, gensp))

 #  new_specs <- merge(species, cor.sp, by.x="scientificName", by.y="gensp", all.x=T)      
#  write.csv(new_specs, "coral_species_new.csv", row.names=F)        

   # Manual corrections are cross-checking with Worms
   # coral$scientificName[coral$scientificName==""]
 
 # Merge with cleaned taxonomy
 require(openxlsx)
 val.sp <- read.xlsx(here("data", "coral_species_new.xlsx"), 1)

  coral2 <- merge(coral, val.sp, by="scientificName")

  # remove nomina dubia
  x <- which(coral2$valid=="dubium")
  coral2 <- coral2[-x,]

  # correct species names (and other information)
  x <- which(!is.na(coral2$validName))
 
  coral2$scientificName[x] <- coral2$validName[x]
  coral2$family.x[x] <- coral2$family.y[x]
  coral2$genus.x[x] <- coral2$genus.y[x]

  # remove superfluent columns
   coral2 <- subset(coral2, select= c(-family.y, -validName, -valid, -scientificNameAuthorship.y, -genus.y))
   corals.OBIS <- coral2
   
   
  # How many valid species in OBIS?
   length(levels(factor(coral2$scientificName)))
   
  # Save the clean file
   save(corals.OBIS, file =here("data", "OBIS_cleaned.RData"))
   
   
#########################
   # Start here 
   here("data", "OBIS_cleaned.RData")

   # filter for z corals 
   obis <- subset(corals.OBIS, symbio=="z")
   # delete last remaining error
   obis <- subset(obis, decimalLatitude > -50)
   

    
    