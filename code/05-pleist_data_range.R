# Set up environment ------------------------------------------------------
# helper packages
library(here)
library(tidyverse)
library(magrittr)

# worms registry
library(worrms)

# load user functions
source(here("code", "functions.R"))


# Data --------------------------------------------------------------------
pleist.df <- read.table(here("data", "original", "range.txt")) %>%
  setNames(c("sp", "bottom", "top")) %>%
  mutate(genusName = sub("(.*)_.*", "\\1", sp),
         speciesName = sub(".*_(.*)", "\\1", sp)) 

pleist.df$genusName <- gsub("-.*","\\1",pleist.df$genusName)
pleist.df$speciesName <- gsub("-.*","\\1",pleist.df$speciesName)
pleist.df$sp <- paste(pleist.df$genusName, pleist.df$speciesName)


n <- c(seq(1,nrow(pleist.df),20), nrow(pleist.df))

wm_sp <- c()

for(i in 1:(length(n)-1)){
  
  ind <- seq(n[i],n[i+1]-1,1)
  
  #create sp df
  temp.sp <- as.data.frame(pleist.df$sp[ind], stringsAsFactors = FALSE)
  colnames(temp.sp) <- "sp"
  temp.sp$new <- NA
  temp.sp$extinct <- NA
  temp.sp$source <- NA
  temp.sp$notes <- NA
  
  temp <- wm_records_names(name = temp.sp$sp)
  
  #check for accepted taxa only
  for (l in 1:length(temp)){
    
    if(nrow(temp[[l]] > 0)){
      nn <- which(temp[[l]]$status == "accepted")
      
      if(length(nn) > 0){
        temp.sp$new[l] <- temp[[l]][nn,]$valid_name
        temp.sp$extinct[l] <- case_when(temp[[l]][nn,]$isExtinct==1~"extinct", 
                                        temp[[l]][nn,]$isExtinct==0~"extant")
        temp.sp$source[l] <- "worms"
      } else temp.sp$notes[l] <- paste(temp[[l]][1,]$status[1], "in worms")
    } 
  }
  
  wm_sp <- rbind.data.frame(wm_sp, temp.sp)
}

#check with pbdb
url <- "https://paleobiodb.org/data1.2/taxa/single.txt?name=%s"

for (i in 1:nrow(wm_sp)){
  if(is.na(wm_sp$new[i]) | is.na(wm_sp$extinct[i])){
    url2 <- sprintf(url, wm_sp$sp[i])
    url2 <- gsub(" ", "%20", url2)
    
    temp <- tryCatch({
      read.csv(url(url2))
    }, error = function(e) {
      cat("\r", i)
    })
    
    if (length(grep("Montastraea", wm_sp$sp[i])) > 0 & is.null(temp)) { #alternate spelling if not found with the other one
      url2 <- sprintf(url, gsub("Montastraea", "Montastrea", wm_sp$sp[i]))
      url2 <- gsub(" ", "%20", url2)
      temp <- tryCatch({
        read.csv(url(url2))
      }, error = function(e) {
        cat("\r", i)
      })
    }
    
    if(!is.null(temp)){
      wm_sp$new[i] <- as.character(temp$accepted_name)
      wm_sp$extinct[i] <- as.character(temp$is_extant)
      wm_sp$source[i] <- "pbdb"
    }
  }
}

wm_sp <- wm_sp[-grep("sp\\.", wm_sp$sp),]
wm_sp <- wm_sp[-grep("aff\\.", wm_sp$sp),]

write.csv(wm_sp, here("data","pleist_sp.csv"), row.names = FALSE)

#MANUAL CHECK PERFORMED HERE

# Resolved species --------------------------------------------------------
#read range data
pleist.df <- read.table(here("data", "original", "range.txt")) %>%
  setNames(c("sp", "bottom", "top")) %>%
  mutate(genusName = sub("(.*)_.*", "\\1", sp),
         speciesName = sub(".*_(.*)", "\\1", sp),
         bottom = ifelse(bottom == 0.1, 0, abs(bottom)),
         top = ifelse(top == 0.1, 0, abs(top))) 

pleist.df$genusName <- gsub("-.*","\\1",pleist.df$genusName)
pleist.df$speciesName <- gsub("-.*","\\1",pleist.df$speciesName)
pleist.df$sp <- paste(pleist.df$genusName, pleist.df$speciesName)

#get accepted species name and extinct status
pleist.sp <- read.csv(here("data", "pleist_sp.csv"), stringsAsFactors = FALSE) %>%
  mutate(new = case_when(is.na(new) & accepted == "yes" & accepted.name==""~sp,
                        accepted.name!=""~accepted.name,
         TRUE~new),
         extinct = extinct.final,
         new=trimws(new)) %>%
  select(sp, new, extinct) %>%
  distinct() %>%
  rename(valid_name = new)

pleist.df <- pleist.df %>% left_join(pleist.sp) %>%
  filter(!is.na(valid_name)) %>%
  group_by(valid_name) %>%
  summarise(bottom=max(bottom), top=min(top), globally.extinct = unique(extinct)) %>%
  mutate(regionally.extinct = ifelse(top > 0, "extinct", "extant"))

##### adding environment #####
pleist.age <- readxl::read_excel(here("data", "original", "budd-pleistocene.xls"), sheet=2) %>%
  janitor::clean_names()

pleist.loc <- readxl::read_excel(here("data", "original", "budd-pleistocene.xls"), sheet=1, skip = 1) %>%
  select(-total) %>% #remove last column
  slice(1:n()-1) %>% #remove last row
  mutate(SpeciesName = sub("^[^_]*_", "", SpeciesName))

pleist.loc$GenusName <- gsub("-.*","\\1",pleist.loc$GenusName)
pleist.loc$SpeciesName <- gsub("-.*","\\1",pleist.loc$SpeciesName)
pleist.loc$sp <- paste(pleist.loc$GenusName , pleist.loc$SpeciesName)

#merge valid_name, location and age
pleist.loc <- pleist.loc %>% left_join(pleist.sp) %>%
  filter(!is.na(valid_name))%>% 
  select(-GenusName, -SpeciesName) %>%
  reshape2::melt(id.vars="valid_name", variable.name = "loc", value.name="n") %>%
  filter(!is.na(n)) %>%
  left_join(pleist.age, by="loc") %>%
  mutate(environment = ifelse(environment==1, "shallow", "deep")) %>%
  group_by(valid_name, environment) %>%
  count(environment) %>%
  ungroup() %>%
  filter(!is.na(environment)) %>%
  reshape2::dcast(valid_name~environment) 

pleist.loc[is.na(pleist.loc)] <- 0

pleist.loc <- pleist.loc %>% mutate(env = case_when(deep > shallow ~ "deep",
                                                  shallow > deep ~ "shallow" ,
                                                  shallow == deep ~ "shallow"
))


pleist.sp2 <- read.csv(here("data", "original", "pleist_sp.csv"), stringsAsFactors = FALSE) %>%
  select(genus, species, name, corallite, shape, max_depth) 

pleist.sp2$genus <- gsub("-.*","\\1",pleist.sp2$genus)
pleist.sp2$species <- gsub("-.*","\\1",pleist.sp2$species)
pleist.sp2$sp <- paste(pleist.sp2$genus  , pleist.sp2$species)

#add to corrected species
pleist.sp2 <- pleist.sp2 %>% mutate(sp = paste(genus, species)) %>%
  left_join(pleist.sp) %>%
  left_join(pleist.loc %>% select(valid_name, env)) %>%
  mutate(max_depth = case_when(env == "deep" ~ 50,
                               env == "shallow" ~ 30,
                               TRUE ~ 0)) %>%
  filter(!is.na(valid_name))


#### branching
pleist.sp2$branching <- plyr::mapvalues(pleist.sp2$shape, from=unique(pleist.sp2$shape), c("HB", "HB", "LB", NA, "LB", "LB", NA, "MB", "MB", "MB", NA))

pleist.df <- pleist.sp2 %>% left_join(pleist.df) %>% select(valid_name, max_depth, branching, corallite, globally.extinct, regionally.extinct)

  ##### get range from pbdb
sp <- pleist.df$valid_name

url <- "https://paleobiodb.org/data1.2/occs/list.txt?base_name=%s&show=paleoloc,coords"

sp_pbdb <- list()

for (i in 1:length(sp)){
  cat("\r", i, "out of", length(sp))
  url2 <- sprintf(url, sp[i])
  url2 <- gsub(" ", "%20", url2)
  
  sp_pbdb[[i]] <-tryCatch({
    read.csv(url(url2)) %>% mutate(valid_name = sp[i])
  }, error = function(e) {
    cat('\r no')
  })
}

sp_pbdb <- sp_pbdb[-which(sapply(sp_pbdb, is.null))]
sp_pbdb <- do.call(rbind, sp_pbdb) %>% filter(max_ma < 5.333 & min_ma > 0.118 & #time constraints
                                                !is.na(paleolng)) #may have to recompute paleolng

#calculate geographic range
library(fields)

sp <- unique(sp_pbdb$valid_name)
great.circle <- c()

for(i in 1:length(sp)){
  temp <- sp_pbdb %>% filter(accepted_name == sp[i]) %>% 
    distinct(paleolng, paleolat)
  
  great.circle[i] <- max(rdist.earth(temp))
}

great.circle[is.infinite(great.circle)] <- NA

#### add proportional range size calculations here
pbdb_all <- chronosphere::fetch("pbdb") %>%
  filter(max_ma < 5.333 & min_ma > 0.118 & !is.na(paleolng)) %>%
  distinct(paleolng, paleolat)

pbdb_range <- max(rdist.earth(pbdb_all))

pleist.df <- pleist.df %>%   
  left_join(cbind.data.frame(sp, range=great.circle, prop_range = great.circle/pbdb_range), by=c("valid_name" = "sp")) %T>% 
  write.csv(here("data", "pleist_resolved.csv"), row.names = FALSE) 
