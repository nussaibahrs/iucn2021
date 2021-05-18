# Set up environment ------------------------------------------------------
# helper packages
library(tidyverse)
library(magrittr)

# setting up machine learning models
library(h2o) # version 3.28.0.2

# packages for explaining our ML models
library(pROC)
library(hmeasure)
library(caret)

source(file.path("code", "functions.R"))

#palette
u_col <- ggsci::pal_uchicago("default")(9)[c(2,5,4,3,1)]

# Data  and set up h20 session -------------------------------------------------------------------
source(file.path("code", "02-modern_corals", "00-load_data.R"))
table(df$extinct)

# Read leaderboard data ------------------------------------------------------
folder_name <- paste0("model_mmm_", hours, "h")
leaderboard <- read.csv(file.path("output", paste0("Table_S_maximised_threshold_model_mmm_", hours, "h.csv")), stringsAsFactors = FALSE) %>% arrange(desc(AUC.test))

# ** Save performance results ----------------------------------------------
#results details
res_summ <- table(sub("\\_.*", "", unique(leaderboard$fname))) %>%
  as.data.frame() %>%
  setNames(c("algorithm", "n"))

leaderboard <- leaderboard[leaderboard$cutoff > 0.2 & leaderboard$AUC.train > 0.7,]

#number of models to choose from
n=nrow(leaderboard)
auc.df <- data.frame(n = 1:n,
                     algorithm =NA,
                     fname=NA,
                     cutoff = NA,
                     AUC.train=NA,
                     AUC.test=NA)

for (i in 1:n){
  auc.df$fname[i] <- leaderboard$fname[i]
  auc.df$cutoff[i] <- leaderboard$cutoff[i]
  mod <- h2o.loadModel(file.path("output", folder_name, leaderboard$fname[i]))
  
  auc.df$algorithm[i] <- mod@algorithm
  
  res.train <- h2o.predict(mod , train)
  res.train <- as.data.frame(res.train)$p1
  
  res.test <- h2o.predict(mod , test)
  res.test <- as.data.frame(res.test)$p1
  
  #all values lower than cutoff value will be classified as 0 (NT in this case)
  train.labels <- ifelse(res.train < leaderboard$cutoff[i], 0, 1)
  test.labels <- ifelse(res.test < leaderboard$cutoff[i], 0, 1)
  true.train.labels <- as.data.frame(train)$extinct
  true.test.labels <- as.data.frame(test)$extinct
  
  auc.df$AUC.train[i] <- as.numeric(auc(true.train.labels, train.labels))
  auc.df$AUC.test[i] <- as.numeric(auc(true.test.labels, test.labels))
}

auc.df <- auc.df %>% 
  group_by(fname) %>%
  filter(AUC.test == max(AUC.test)) %>%
  arrange(desc(AUC.test))

auc.df <- auc.df[duplicated(auc.df$fname)==FALSE,] #remove duplicates

#checking for precision and recall
train.precision <- c()
test.precision <- c()

train.recall <- c()
test.recall <- c()

for (w in 1:10){ #top ten models
  win <- auc.df$n[w]
  aml_leader <- h2o.loadModel(file.path("output", folder_name, auc.df$fname[w]))
  
  res.train <- h2o.predict(aml_leader , train)
  res.train <- as.data.frame(res.train)$p1
  
  res.test <- h2o.predict(aml_leader , test)
  res.test <- as.data.frame(res.test)$p1xxxnhj
}

train.Fscore = 2 * train.precision * train.recall / (train.precision + train.recall)
test.Fscore = 2 * test.precision * test.recall / (test.precision + test.recall)

auc.df <- cbind(auc.df[1:10,], train.Fscore=train.Fscore, test.Fscore=test.Fscore) %>% 
  arrange(AUC.test)

write.csv(auc.df, file.path("output", "Table_S_model_performance_mmm.csv"), row.names=FALSE)
