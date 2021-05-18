##### RESOLVING TAXONOMY WITH THE WORMS DATABASE ####
library(tidyverse)
library(robis)
library(rworldmap)
library(rgeos)
library(sp)
library(rgdal)

source(file.path("code", "functions.R"))

#### IUCN
iucn <- read.csv(file.path("data", "traits_iucn.csv"), stringsAsFactors = FALSE) 

#### CoralTraits Database
traits.corals <- read.csv(file.path("data", "ctdb_1.1.1_data.csv"), stringsAsFactors = FALSE)

#Check names with WoRMS
namesToMatch <- traits.corals$specie_name %>% unique()
length(namesToMatch)

n <- 20 #doesn't work with a certain number of names e.g. > 30
chunks <- c(seq(0,length(namesToMatch) - length(namesToMatch)%%n, n), length(namesToMatch))

valid <- list()

for (i in 1:(length(chunks)-1)){
  temp <- namesToMatch[(chunks[i]+1):chunks[i+1]]
  tx <- worrms::wm_records_names(name=temp)
  
  for (t in 1:length(tx)){
    if(nrow(tx[[t]]) > 0) tx[[t]] <- tx[[t]][,c("scientificname", "valid_name")] else tx[[t]]<- cbind.data.frame(
      scientificname = temp[t],valid_name=NA)
    
    }
  valid <- append(valid, tx)
}

valid <- do.call(rbind.data.frame, valid)

#check for NAs in valid names
na.tx <- as.data.frame(valid[is.na(valid$valid_name),])

na.tx$manual <- match_taxa(na.tx$scientificname, verbose=FALSE, chunk=20)
write.csv(na.tx, file.path("data", "ctdb_na.csv"), row.names = FALSE)

# after manual check
na.tx <- read.csv(file.path("data", "ctdb_na.csv"), stringsAsFactors = FALSE) %>% 
  filter(!is.na(manual))
na.tx$manual <- sub("\\(.*\\) ", "", na.tx$manual)
na.tx$valid_name <- na.tx$manual


valid <- rbind(valid[!is.na(valid$valid_name),],
      na.tx[,c(1,2)])

traits.corals$valid <- plyr::mapvalues(traits.corals$specie_name, from=valid$scientificname, to=valid$valid_name)


write.csv(traits.corals %>% filter(!is.na(valid)), file.path("data", "ctdb_resolved.csv"), row.names =FALSE)

#### START file.path FOR MERGING ####

###### geographic range
library(fields)
library(sp)
library(file.path)
library(tidyverse)
library(magrittr)

obis.corals <- read.csv(file.path("data", "2019-11-05_obis_scleractinia.csv"),
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
iucn.corals <- read.csv(file.path("data", "traits_iucn.csv"), stringsAsFactors = FALSE)
head(iucn.corals)

selected_traits <- c("Depth lower",
                     "Growth form typical", "Corallite width maximum")

traits.corals <- read.csv(file.path("data", "ctdb_resolved.csv"), stringsAsFactors = FALSE)

sort(unique(traits.corals$trait_name))
unique(traits.corals[traits.corals$trait_name == "Corallite width maximum",]$standard_unit) #all in cm


traits.corals <- read.csv(file.path("data", "ctdb_resolved.csv"), stringsAsFactors = FALSE) %>%
  
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
#   arrange(n) #%T>% write.csv(file.path("data", "growth_forms.csv"), row.names = FALSE)

growth <- read.csv(file.path("data", "growth_forms.csv"), stringsAsFactors = FALSE)

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
colnames(traits.up)[1] <- "sp"
#checking for duplicates

#### merge
df.corals <- traits.up %>% left_join(iucn.corals, by=c("sp"="valid.name")) %>%
  filter(!is.na(redlistCategory)) %>%
  left_join(as.data.frame(great_circle, stringAsFactors = FALSE), by=c("sp"="sp"))

#save data
df.corals <- df.corals[,c("sp", "redlistCategory", x, selected_traits[1], "range")] %>%
  setNames(c("sp", "iucn", x, "max_depth", "branching", "range")) %T>%
  write.csv(file.path("data", "traits_iucn_mdmm.csv"), row.names = FALSE)

