# Set up environment ------------------------------------------------------
# helper packages
library(tidyverse)
library(magrittr)

# setting up machine learning models
library(h2o)


# Data  and set up h20 session -------------------------------------------------------------------
source(file.path("code", "01-model_fmm", "00-load_data.R"))

# * Automatic Machine Learning ----------------------------------------------
system.time(aml <- h2o.automl(y=y, x=x, 
                              training_frame = train,
                              keep_cross_validation_fold_assignment = TRUE,
                              keep_cross_validation_predictions = TRUE,
                              seed = seed.no,
                              stopping_metric = "AUC",
                              balance_classes = TRUE,
                              nfolds=5,
                              max_runtime_secs = run_time))


# * Save models -----------------------------------------------------------
folder_name <-folder 
unlink(file.path("output", folder_name), recursive = TRUE)

leaderboard <- as.data.frame(aml@leaderboard)

for(i in 1:nrow(leaderboard)) {
  aml1 <- h2o.getModel(aml@leaderboard[i, 1]) # get model object in environment
  h2o.saveModel(object = aml1, file.path("output", folder_name)) # pass that model object to h2o.saveModel as an argument
  
}

write.csv(leaderboard, file.path("output", folder_name, "leaderboard.csv"), row.names = FALSE)

# Find best model ---------------------------------------------------------
folder_name <- folder
leaderboard <- read.csv(file.path("output", folder_name, "leaderboard.csv"), stringsAsFactors = FALSE)

mods <- list() 

for(i in 1:nrow(leaderboard)){
  mods[[i]] <- h2o.loadModel(file.path("output", folder_name, leaderboard$model_id[i]))
}

# * Find probability threshold  -----------------
cutoff=seq(0, 1, 0.01) #threshold values
n=nrow(leaderboard) #number of models

opt_par <- list()

for (i in 1:n){
  temp <- as.data.frame(matrix(0, nrow=(length(cutoff)-1), ncol=10))
  colnames(temp) <- c("fname", "model", "cutoff.train", "cutoff.test", "youden.train", "youden.test", "AUC.train",
                      "AUC.test", "misclass_rate", "h")
  temp$model <-mods[[i]]@algorithm
  temp$fname <- mods[[i]]@model_id
  
  for (c in 2:length(cutoff)){
    
    res.train <- h2o.predict(mods[[i]], train)
    res.train <- as.data.frame(res.train)$p1
    
    res.test <- h2o.predict(mods[[i]], test)
    res.test <- as.data.frame(res.test)$p1
    
    #all values lower than cutoff value will be classified as 0 (NT in this case)
    train.labels <- ifelse(res.train < cutoff[c], 0, 1)
    test.labels <- ifelse(res.test < cutoff[c], 0, 1)
    true.train.labels <- as.data.frame(train)$extinct
    true.test.labels <- as.data.frame(test)$extinct
    
    met.train <- hmeasure::HMeasure(true.train.labels,train.labels) 
    met.test <- hmeasure::HMeasure(true.test.labels,test.labels)   
    
    
    temp[c,c("cutoff.train", "cutoff.test")] <- cutoff[c]
    
    #error metrics
    mod.roc.train <- pROC::roc(true.train.labels, train.labels) 
    mod.roc.test <- pROC::roc(true.test.labels, test.labels) 
    #performance
    temp[c, c("AUC.train", "AUC.test", 
              "misclass_rate", "youden.train", "youden.test", "h")] <- c(pROC::auc(mod.roc.train), #AUC
                                                                         pROC::auc(mod.roc.test), #AUC
                                                                         met.test$metrics$ER, #misclassification rate = 1 - Accuracy
                                                                         met.train$metrics$Youden,
                                                                         met.test$metrics$Youden, # Youden Index
                                                                         met.test$metrics$H
              )
  }
  opt_par[[i]] <- temp
  
}

opt_par <- do.call(rbind, opt_par)

#maximisation based on youden index
temp <- opt_par %>% group_by(fname) %>% summarise(youden.train=max(youden.train), youden.test=max(youden.test)) %>%
  na.omit()

opt_train <- temp %>% dplyr::select(fname, youden.train) %>% 
  left_join(opt_par %>% dplyr::select(fname, AUC.train, youden.train, cutoff.train), by=c("fname", "youden.train"))

opt_test <- temp %>% dplyr::select(fname, youden.test) %>% 
  left_join(opt_par %>% dplyr::select(fname, AUC.test, youden.test, cutoff.test), by=c("fname", "youden.test"))

opt_par_final <- full_join(opt_train, opt_test)%>%
  group_by(fname) %>%
  mutate(cutoff = (cutoff.test+ cutoff.train)/2) %>%
  arrange(desc(AUC.test), desc(AUC.train))  %T>%
  write.csv(file.path("output", paste0("Table_S_maximised_threshold_model_fmm_", hours, "h.csv")), row.names = FALSE) 


# Terminate h2o session ---------------------------------------------------
h2o.shutdown(prompt=FALSE)
