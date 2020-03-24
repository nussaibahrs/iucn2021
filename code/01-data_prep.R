##### RESOLVING TAXONOMY WITH THE WORMS DATABASE ####
library(here)
library(tidyverse)
library(robis)
library(rworldmap)
library(rgeos)
library(sp)
library(rgdal)

source(here("code", "functions.R"))

#### IUCN
iucn <- read.csv(here("data", "2019-11-04_iucn.csv"), stringsAsFactors = FALSE) %>%
  filter(phylumName == "CNIDARIA")
head(iucn)

#Check names with WoRMS
namesToMatch <- iucn$scientificName

valid <- match_taxa(namesToMatch, verbose=FALSE)
iucn$valid <- valid

iucn[which(is.na(valid)),]$valid <- iucn[which(is.na(valid)),]$scientificName #NA generated one keep original sci. name
write.csv(iucn, here("data", "iucn_resolved.csv"), row.names = FALSE)

#Get OBIS Data
#obis.corals<- occurrence("Scleractinia")
load(here("data", "original", "2019-11-13_obis_scleractinia.RData"))

# depth < 150
# obis.corals <- obis.corals %>% filter(depth <=150)

# check occurrences falling on land 
w <- chronosphere::fetch("NaturalEarth", "land")[1,] #continents only 
dat <- obis.corals %>% mutate(row_n = 1:nrow(.))
coordinates(dat) <- ~ decimalLongitude + decimalLatitude
proj4string(dat) <- proj4string(w)

land_points <- over(w, dat)

obis.corals <- obis.corals[-as.numeric(na.omit(land_points$row_n)),]

#clean data according to Kiessling
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

# download only those in iucn from obis
iucn.corals <- read.csv(here("data", "iucn_resolved.csv"), stringsAsFactors = FALSE)
obis.corals <- subset(coral2, scientificName %in% iucn.corals$valid)


write.csv(obis.corals %>% dplyr::select(scientificName, decimalLongitude, decimalLatitude, depth, symbio), here("data", "2019-11-05_obis_scleractinia.csv"), row.names = FALSE)

#### CoralTraits Database
traits.corals <- read.csv(here("data", "original", "ctdb_1.1.1_data.csv"), stringsAsFactors = FALSE)

#Check names with WoRMS
namesToMatch <- traits.corals$specie_name %>% unique()
length(namesToMatch)

valid <- match_taxa(namesToMatch, chunk = 10, verbose = FALSE)

traits.corals$valid <- plyr::mapvalues(traits.corals$specie_name, from=namesToMatch, to=valid)

write.csv(traits.corals, here("data", "ctdb_resolved.csv"), row.names =FALSE)

#### START HERE FOR MERGING ####

###### geographic range
library(fields)
library(sp)
library(here)
library(tidyverse)
library(magrittr)

obis.corals <- read.csv(here("data", "2019-11-05_obis_scleractinia.csv"),
                        stringsAsFactors = FALSE) 


scle.sp <- unique(obis.corals$scientificName)
length(scle.sp)

great_circle <- c()

for (i in 1:length(scle.sp)){
  temp <- obis.corals[obis.corals$scientificName == scle.sp[i],] 
  
  #returns the maximum great circle distance for species | mean radius of earth is 6371 km 
  great_circle[i] <- max(fields::rdist.earth(unique(temp[,-1]), 
                                     miles=FALSE)) 
}

great_circle <- cbind(sp=scle.sp, range=great_circle)


#merge data
iucn.corals <- read.csv(here("data", "iucn_resolved.csv"), stringsAsFactors = FALSE)
head(iucn.corals)

selected_traits <- c("Depth lower",
                     "Growth form typical", "Corallite width maximum")

traits.corals <- read.csv(here("data", "ctdb_resolved.csv"), stringsAsFactors = FALSE)

sort(unique(traits.corals$trait_name))
unique(traits.corals[traits.corals$trait_name == "Corallite width maximum",]$standard_unit) #all in cm


traits.corals <- read.csv(here("data", "ctdb_resolved.csv"), stringsAsFactors = FALSE) %>%
  
  filter(
    #choose specific traits
    trait_name %in% c(selected_traits,"Zooxanthellate"),
    
    #choose only those in iucn
    valid %in% iucn.corals$valid)

#choose only symbiotic corals
corals.z <- traits.corals %>% filter(value %in% c("zooxanthellate")) %>% 
  pull(specie_name) %>% unique()

length(corals.z)

traits.corals <- traits.corals[traits.corals$specie_name %in% corals.z,]
head(traits.corals)

### growth forms
# traits.corals %>% filter(trait_name == "Growth form typical") %>%
#   group_by(value) %>% tally() %>%
#   arrange(n) #%T>% write.csv(here("data", "growth_forms.csv"), row.names = FALSE)

growth <- read.csv(here("data", "growth_forms.csv"), stringsAsFactors = FALSE)

#replace with new branching categories
traits.corals[traits.corals$trait_name == "Growth form typical",]$value <- plyr::mapvalues(traits.corals[traits.corals$trait_name == "Growth form typical",]$value, 
                 from=growth$value, to=growth$branching)

table(traits.corals[traits.corals$trait_name == "Growth form typical",]$value)

###### generate traits table for selected traits
traits.sp <- unique(traits.corals$valid)
traits.up <- as.data.frame(matrix(ncol=4, nrow=length(traits.sp)))
colnames(traits.up) <- c("sp", selected_traits)

for(i in 1:length(traits.sp)){
  temp <- traits.corals[traits.corals$valid == traits.sp[i],]
  
  traits.up$sp[i] <- traits.sp[i]
  
  for (j in selected_traits){
    pl <- which(temp$trait_name == j)
    
    if (length(pl) > 0){
      #maximum depth or corallite width
      if(j == "Depth lower"){ 
        traits.up[i,j] <- max(as.numeric(temp[pl,]$value), na.rm = TRUE) #maximum in case of multiple
      }
      
      if(j == "Corallite width maximum"){ 
        traits.up[i,j] <- max(as.numeric(temp[pl,]$value), na.rm=TRUE) #maximum in case of multiple
      }
      
      
      if(j == "Growth form typical"){ 
        if (length(unique(temp$value[pl])) > 1) { #check if they are all the same
          traits.up[i,j] <- temp[pl,] %>% sample_n(1) %>% pull(value)
        } else {
          traits.up[i,j] <- unique(temp[pl,]$value) #get single branching from multiple
        }
      }
    }
    
  }
  
}

# binary form of depth 
#traits.up$`Depth lower` <- ifelse(traits.up$`Depth lower` <=30, "shallow", "deep")
head(traits.up)

#checking for duplicates
single_entries <- iucn.corals %>% distinct(scientificName,valid, redlistCategory)

dupl <- single_entries %>%
  group_by(valid) %>% summarise(n=length(unique(redlistCategory))) %>%
  filter(n > 1) %>% pull(valid)

# resolving for synonynmisation duplication - 
# assign the lowest threat level and ditch all others
#conflicting categories
iucn_order <- unique(iucn.corals$redlistCategory)[c(2,1,4,3,5,6)]

for (i in 1:length(dupl)){
  temp_status <- single_entries %>% filter(valid == dupl[i] & redlistCategory != "Data Deficient") %>% 
    distinct(redlistCategory) %>% pull(redlistCategory)
  
  if (length(temp_status) > 1) {
    temp_sp <- single_entries %>% filter(valid == dupl[i])
    
    #check threat level
    temp_thr <- temp_sp %>% group_by(scientificName) %>% mutate(n = which(iucn_order == redlistCategory),
                                                                valid_n = ifelse(scientificName == valid, 1, 0))
    
    single_entries[single_entries$valid == dupl[i],]$redlistCategory <- iucn_order[min(temp_thr$n)]
  } else {
    single_entries[single_entries$valid == dupl[i],]$redlistCategory <- temp_status
  }
}

iucn.corals <- iucn.corals %>%
  select(valid, redlistCategory) %>%
#remove all duplicates
 filter(!valid %in% dupl) %>%
  #add updated status of duplicates
  bind_rows(single_entries %>% select(-scientificName)) %>%
  distinct()

#### merge
df.corals <- traits.up %>% left_join(iucn.corals, by=c("sp"="valid")) %>%
  filter(!is.na(redlistCategory)) %>%
  left_join(as.data.frame(great_circle, stringAsFactors = FALSE), by=c("sp"="sp"))

#check for NAs in traits
na.count <- apply(df.corals[-c(1,5)],1, function(x) {sum(is.na(x))})
df.corals <- df.corals[-which(na.count > 2),]

#save data
df.corals <- df.corals[,c("sp", "redlistCategory", selected_traits, "range")] %>%
  setNames(c("sp", "iucn", "max_depth", "branching", "corallite", "range")) %T>%
  write.csv(here("data", "traits_iucn.csv"), row.names = FALSE)

