library(here)
library(tidyverse)

coral_traits <- read.csv(here("data", "ctdb_resolved.csv"), stringsAsFactors = FALSE) %>%
  filter(trait_name %in% c("Growth rate", "Colony maximum diameter"))

unique(coral_traits$standard_unit)

#### growth 

coral_traits <- coral_traits %>% 
  mutate(value <- as.numeric(value),
         
         #convert to mm/yr & mm
         value2 = case_when(standard_unit=="mm yr^-1" ~ value,
                            standard_unit=="mm d^-1"~ value/10,
                            standard_unit=="mm month^-1" ~ value*12,
                            standard_unit=="mm per 7 months^-1"~ (value/7)*12,
                            standard_unit== "mm per 6 months^-1"~ value*2,
                            standard_unit=="mm per 4 months^-1"~ value*3, 
                            standard_unit=="mm per 3 months^-1"~ value*4, 
                            standard_unit=="mm per 8 months^-1"~ (value/8)*12,
                            standard_unit=="mm" ~ value,
                            standard_unit=="cm" ~ value*10
                            )
         )


coral_age <- coral_traits %>%
  group_by(specie_name, trait_name) %>%
  
  #get mean colony size and growth rate
  summarise(mean=mean(value2)) %>%
  pivot_wider(names_from = trait_name,
              values_from=mean) %>%
  ungroup() %>%
  setNames(c("sp", "size", "growth_rate")) %>%
  mutate(age=size/growth_rate) %>%
  filter(!is.na(age)) 

median(coral_age$age)
mean(coral_age$age)

range(coral_age$age)
n=200
hist(coral_age[coral_age$age < n,]$age, breaks = seq(0,n, 10))

########

coral_age %>% filter(age < 100) %>%
  summarise(g=mean(age))

coral_age2[,1]/coral_age2[,2]
