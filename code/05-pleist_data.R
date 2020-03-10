library(here)
library(tidyverse)
library(magrittr)

source(here("code", "functions.R"))

#### Data prep: Pleistocene Data ####
pleist.loc <- readxl::read_excel(here("data", "original", "budd-pleistocene.xls"), sheet=1, skip = 1) %>%
  select(-total) %>% #remove last column
  slice(1:n()-1) %>% #remove last row
  mutate(SpeciesName = sub("^[^_]*_", "", SpeciesName),
         sp = as.character(paste(GenusName, SpeciesName)))

head(pleist.loc)
pleist.loc$sp

for (t in c("-ann", "-cav", "-mass", "-br", "-2", "aff\\.", ".cf")){
  pleist.loc$sp <- sub(t, "", pleist.loc$sp)
}

# check names in worms
sp_update <- match_taxa(pleist.loc$sp, chunk = 20, verbose = FALSE) 

sp <- cbind.data.frame(old=pleist.loc$sp, new=sp_update %>% replace_na(""))
sp <- sp %>% filter(new == "") #%T>% write.csv(here("data", "coral_indet.csv")) #those that were not matched in Worms

pleist.loc$sp2 <- sp_update

#check with pbdb
url <- "https://paleobiodb.org/data1.2/taxa/single.txt?name=%s"

sp_pbdb <- list()

for (i in 1:nrow(sp)){
  url2 <- sprintf(url, sp$old[i])
  url2 <- gsub(" ", "%20", url2)
  
  sp_pbdb[[i]] <-tryCatch({
    read.csv(url(url2))
  }, error = function(e) {
    cat('no')
  })
}

sp_pbdb <- sp_pbdb[-which(sapply(sp_pbdb, is.null))]
sp_pbdb <- do.call(rbind, sp_pbdb)

#add the accepted name
for(i in 1:nrow(sp_pbdb)){
  pleist.loc$sp2[pleist.loc$sp == sp_pbdb$taxon_name[i]] <- as.character(sp_pbdb$accepted_name[i])
}

##### still unknown, need to check
pleist.loc[is.na(pleist.loc$sp2),c("GenusName", "SpeciesName", "sp", "sp2")] <- c("Antillia dentata", "Caulastraea portoricensis", "Agaricia undata", "Montastraea brevis", "Montastraea trinitatis",
  "Undaria pusilla")

pleist.loc[,c("GenusName", "SpeciesName", "sp2")] %>%
  setNames(c("genus", "species", "name")) %T>%
  write.csv(here("data", "pleist_sp.csv"), row.names = FALSE)

# check traits from here: https://fossils.its.uiowa.edu/database/corals/combined/

##### traits: branching, corallite diameter, range, Z or AZ  
ctdb <- read.csv(here("data", "ctdb_resolved.csv"), stringsAsFactors = FALSE)

sp_ctdb <- pleist.loc$sp2[which(pleist.loc$sp2 %in% ctdb$valid & !is.na(pleist.loc$sp2)),]
length(sp_ctdb)

#not in the database
pleist.loc$sp2[which(!pleist.loc$sp2 %in% ctdb$valid & !is.na(pleist.loc$sp2)),] %>% length()
  
selected_traits <- c("Depth lower",
                     "Growth form typical", "Corallite width maximum")

traits.corals <- ctdb %>% filter(valid %in% sp_ctdb) %>%
  filter(trait_name %in% c(selected_traits, "Zooxanthellate"))

#choose only symbiotic corals
corals.z <-traits.corals %>% filter(value %in% c("zooxanthellate", "both")) %>% 
  pull(specie_name) %>% unique()

length(corals.z)

traits.corals <- traits.corals[traits.corals$specie_name %in% corals.z,]

### growth forms
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

##### adding environment #####
pleist.age <- readxl::read_excel(here("data", "original", "budd-pleistocene.xls"), sheet=2) %>%
  janitor::clean_names()

pleist.loc <- readxl::read_excel(here("data", "original", "budd-pleistocene.xls"), sheet=1, skip = 1) %>%
  select(-total) %>% #remove last column
  slice(1:n()-1) %>% #remove last row
  mutate(SpeciesName = sub("^[^_]*_", "", SpeciesName),
         sp = as.character(paste(GenusName, SpeciesName)))

#merge location and age
pleist.df <- pleist.loc %>% select(-GenusName, -SpeciesName) %>%
  reshape2::melt(id.vars="sp", variable.name = "loc", value.name="n") %>%
  filter(!is.na(n)) %>%
  left_join(pleist.age, by="loc") %>%
  mutate(environment = ifelse(environment==1, "shallow", "deep")) %>%
  group_by(sp, environment) %>%
  tally() %>%
  ungroup() %>%
  reshape2::dcast(sp~environment) 

pleist.df[is.na(pleist.df)] <- 0

pleist.df <- pleist.df %>% mutate(env = case_when(deep > shallow ~ "deep",
                                     shallow > deep ~ "shallow" ,
                                     shallow == deep ~ "shallow"
                     ))


pleist.sp <- read.csv(here("data", "pleist_sp2.csv"), stringsAsFactors = FALSE) %>%
  select(-comments)


#add to corrected species
pleist.df <- pleist.sp %>% mutate(sp = paste(genus, species)) %>%
  left_join(pleist.df, by="sp") %>%
  mutate(max_depth = case_when(env == "deep" ~ 50,
                               env == "shallow" ~ 30,
                               TRUE ~ 0))


#### branching
pleist.df$branching <- plyr::mapvalues(pleist.df$shape, from=unique(pleist.df$shape), c("HB", "HB", "LB", "LB", "LB", "LB", "MB", "MB", "MB"))

pleist.df <- pleist.df %>% select(name, corallite, branching, max_depth)

##### get range from pbdb
sp <- pleist.df$name

url <- "https://paleobiodb.org/data1.2/occs/list.txt?base_name=%s&show=paleoloc,coords"

sp_pbdb <- list()

for (i in 1:length(sp)){
  cat("\r", i, "out of", length(sp))
  url2 <- sprintf(url, sp[i])
  url2 <- gsub(" ", "%20", url2)
  
  sp_pbdb[[i]] <-tryCatch({
    read.csv(url(url2))
  }, error = function(e) {
    cat('no')
  })
}

sp_pbdb <- sp_pbdb[-which(sapply(sp_pbdb, is.null))]
sp_pbdb <- do.call(rbind, sp_pbdb) %>% filter(max_ma < 5.4 & #time constraints
                                                !is.na(paleolng)) #may have to recompute paleolng

#calculate geographic range
library(fields)

sp <- unique(sp_pbdb$accepted_name)
great.circle <- c()

for(i in 1:length(sp)){
  temp <- sp_pbdb %>% filter(accepted_name == sp[i]) %>% 
    distinct(paleolng, paleolat)
  
  great.circle[i] <- max(rdist.earth(temp))
}


#### add proportional range size calculations here
pbdb_all <- read.csv("https://paleobiodb.org/data1.2/occs/list.csv?max_ma=5.333&min_ma=0.117&show=paleoloc") %>%
  filter(!is.na(paleolng)) %>%
  distinct(paleolng, paleolat)

pbdb_range <- max(rdist.earth(pbdb_all))

pleist.df <- pleist.df %>%
  left_join(cbind.data.frame(sp, range=great.circle, prop_range = great.circle/pbdb_range), by=c("valid_name" = "sp")) %T>% 
  write.csv(here("data", "pleist_resolved.csv"), row.names = FALSE) 