setwd("/home/nussaibah/Documents/Projects/TERSANE/coralTrace/traits_iucn")

library(here)
library(tidyverse)
library(caret)
#library(doMC)
#library(parallel)

# configure multicore
#n_cores <- detectCores()
#registerDoMC(cores = 4)

### load data
df.corals <- read.csv(here("data", "traits_iucn.csv"), stringsAsFactors = FALSE) # %>%
# mutate(corallite = ifelse(corallite > 30, NA, corallite))

#number of NAs in each column
100- apply(df.corals, 2, function(x) sum(is.na(x))) / nrow(df.corals)  *100

# Get rid of data deficient corals
df.corals <- df.corals %>%
  na.omit() %>%
  
  #omit species name
  select(-sp) %>%
  
  #remove all Data Deficient corals
  filter(iucn != "Data Deficient") %>%
  
  #treat as binary
  mutate(
    #reclassify iucn status 
    iucn = case_when(iucn == "Least Concern" ~ "LC",         
                     iucn == "Vulnerable" ~ "V",           
                     iucn == "Near Threatened"   ~ "NT",    
                     iucn == "Endangered"   ~ "E",         
                     iucn == "Critically Endangered"~ "E"),
    branching = factor(branching, levels=c("NB", "LB", "MB", "HB"))) 

#set seed for reproducibility

#divide between training and test
set.seed(42)
smp <- sample(1:nrow(df.corals), round(nrow(df.corals) * 0.75))

train <- df.corals[smp,]
test <- df.corals[-smp,]

table(df.corals$branching)
table(df.corals$iucn)

#preprocessing
preProcValues <- preProcess(df.corals,
                            method= c("center", "scale", "nzv")
)

preProcValues

trainTransformed <- predict(preProcValues, train)
testTransformed <- predict(preProcValues, test)

#cross validation
# set control parameters for training 
reps <- 5
cv <- 5
metric <- "AUC"

trCtrl <- trainControl(method = "repeatedcv", 
                       number = cv, 
                       repeats = reps, 
                       classProbs = TRUE,
                       summaryFunction = multiClassSummary,
                       sampling = "up",
                       search = "grid",
                       savePredictions = T)

#####  Boosted Trees ####
gbmGrid <-  expand.grid(interaction.depth = 1, #maximum depth of trees, i.e. highest level of variable interactions allowed
                        n.trees = (10:20)*50, #number of trees
                        shrinkage = seq(0.1, 1, 0.1), #learning rate
                        n.minobsinnode = 10:30 #minimum number of observations in the terminal nodes of the trees
)

set.seed(42)
gbmFit <- train(iucn ~ ., data = trainTransformed, 
                method = "gbm", 
                metric = metric,
                trControl = trCtrl, 
                tuneGrid = gbmGrid)

save(gmbFit, file = here("output", "ML_multi_gbm.RData"))
