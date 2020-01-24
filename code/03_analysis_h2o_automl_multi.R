# Set up environment ------------------------------------------------------
# helper packages
library(here)
library(tidyverse)

# setting up machine learning models
library(h2o)

# Data  -------------------------------------------------------------------

# classification data
df.corals <- read.csv(here("data", "traits_iucn.csv"), stringsAsFactors = FALSE) %>% mutate(range = ifelse(range > 0, range, NA))
df.corals <- df.corals %>% mutate_if(is.numeric, function (x) as.numeric(scale(x))) #normalise numerical predictors

df <- df.corals %>%
  na.omit() %>%
  
  #omit species name
  select(-sp) %>%
  
  #remove all Data Deficient corals
  filter(iucn != "Data Deficient") %>%
  
  #treat as multi-class
  mutate(
    #reclassify iucn status 
    iucn = case_when(iucn == "Least Concern" ~ "LC",         
                     iucn == "Vulnerable" ~ "V",           
                     iucn == "Near Threatened"   ~ "NT",    
                     iucn == "Endangered"   ~ "E",         
                     iucn == "Critically Endangered"~ "E"),
    branching = factor(branching, levels=c("NB", "LB", "MB", "HB")),
    iucn = factor(iucn, levels=c("LC", "NT", "V", "E")))  

# Set up h2o session ------------------------------------------------------

# initialize h2o session
h2o.init(max_mem_size = "16g", nthreads=6)

# convert to h2o object
df.h2o <- as.h2o(df)

seed.no <- 42 #for reproducibility

# Split h2o data into training, validation, and test frames
data.split <- h2o.splitFrame(df.h2o,ratios = .75, seed=seed.no)

train <- data.split[[1]]  # For training
test <- data.split[[2]] # For final evaluation of model performance

# variable names for response & features
y <- "iucn"
x <- setdiff(names(df), y)


# Models ------------------------------------------------------------------

run_time <- 30 #maximum time to run automl

# * Automatic Machine Learning ----------------------------------------------
system.time(aml <- h2o.automl(y=y, x=x, 
                              training_frame = train,
                  # keep_cross_validation_fold_assignment = TRUE,
                  # keep_cross_validation_predictions = TRUE,
                  seed = seed.no,
                  # balance_classes = TRUE,
                  # nfolds=10,
                  max_runtime_secs = run_time))


# * Save models -----------------------------------------------------------

folder_name <-"model2"
unlink(here("output", folder_name), recursive = TRUE)

leaderboard <- as.data.frame(aml@leaderboard)

for(i in 1:nrow(leaderboard)) {
  aml1 <- h2o.getModel(aml@leaderboard[i, 1]) # get model object in environment
  h2o.saveModel(object = aml1, here("output", folder_name)) # pass that model object to h2o.saveModel as an argument
  
}

write.csv(leaderboard, here("output", folder_name, "leaderboard.csv"), row.names = FALSE)


# Terminate h2o session ---------------------------------------------------
h2o.shutdown()

