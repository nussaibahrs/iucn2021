library(here)
library(tidyverse)
library(caret)

##### load model 
load(here("output", "ML_binary.RData"))
final <- mods[[6]] #1:neural network, 6:gbm

#### PLEISTOCENE ####
pleist.df <- read.csv(here("data", "pleist_resolved.csv"), stringsAsFactors = FALSE)

#### process data
pleist.df <- predict(preProcValues, pleist.df) %>% 
  select (-range) %>%
  rename(range=prop_range) #choose proportional range


pred <- predict(final, newdata=na.omit(pleist.df))

##### add extinction
pleist.loc <- readxl::read_excel(here("data", "original", "budd-pleistocene.xls"), sheet=1, skip = 1) %>%
  select(-total) %>% #remove last column
  slice(1:n()-1) %>% #remove last row
  mutate(SpeciesName = sub("^[^_]*_", "", SpeciesName),
         sp = as.character(paste(GenusName, SpeciesName)))

pleist.age <- readxl::read_excel(here("data", "original", "budd-pleistocene.xls"), sheet=2) %>%
  janitor::clean_names()

pleist.env <- pleist.loc %>% select(-GenusName, -SpeciesName) %>%
  reshape2::melt(id.vars="sp", variable.name = "loc", value.name="n") %>%
  filter(!is.na(n)) %>%
  left_join(pleist.age, by="loc") 
  
pleist.sp <- read.csv(here("data", "pleist_sp.csv")) %>%
  select(genus, species, name) %>%
  mutate(sp=paste(genus, species),
         name = trimws(name))

pleist.env %>%
  left_join(pleist.sp %>% mutate(name = trimws(name)), by=c("sp"="sp")) %>%
  group_by(name) %>%
  summarise(bottom = max(age_start),
         top = min(age_end), 
         mid = (bottom + top)/2) %>%
  left_join(cbind(na.omit(pleist.df), pred) %>%
              select(name, pred)) %>%
  filter(!is.na(pred)) %>%
  View()


# Data-deficient ----------------------------------------------------------
### load data
df.corals <- read.csv(here("data", "traits_iucn.csv"), stringsAsFactors = FALSE)  %>%
  #remove all Data Deficient corals
  filter(iucn == "Data Deficient")

dd.pred <- na.omit(df.corals)

#preprocessing data
dd.pred <- predict(preProcValues, dd.pred)

#predict status
dd.pred$status <- as.character(predict(final, dd.pred))

dd.pred <- df.corals %>% select(sp) %>%
  left_join(dd.pred %>% select(sp, status)) %>%
  mutate(status = ifelse(is.na(status), "DD", status)) %>%
  distinct(sp, status)

barplot(table(dd.pred$status))



