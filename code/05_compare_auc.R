# Set up environment ------------------------------------------------------
# helper packages
library(here)
library(tidyverse)
library(magrittr)

# setting up machine learning models
library(h2o) # version 3.28.0.2

# packages for explaining our ML models
library(pROC)
library(verification)

#rescaling
scaled.new <- function (newdata, olddata) {(newdata-min(olddata, na.rm = TRUE))/(max(olddata, na.rm = TRUE)-min(olddata, na.rm = TRUE))}
unscale <- function(newdata, olddata){min(olddata, na.rm = TRUE) + newdata * (max(olddata, na.rm = TRUE) - min(olddata, na.rm=TRUE))}

# Data  -------------------------------------------------------------------

# classification data
df.corals <- read.csv(here("data", "traits_iucn.csv"), stringsAsFactors = FALSE) %>% mutate(
  range = ifelse(range > 0, range, NA),
  corallite_ori = corallite,
  depth_ori = max_depth,
  #max_depth = ifelse(max_depth <=30, "shallow", "deep") %>% as.factor(),
  range_ori = range)

#df.corals <- df.corals %>% mutate_at(c("max_depth", "corallite", "range"), function (x) as.numeric(scale(x))) #normalise numerical predictors
df.corals <- df.corals %>% mutate_at(c("max_depth", "corallite", "range"), function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))})

df <- df.corals %>%
  na.omit() %>%
  
  #remove all Data Deficient corals
  filter(iucn != "Data Deficient") %>%
  
  #treat as binary
  mutate(
    status =iucn,
    #reclassify iucn status
    iucn = case_when(iucn == "Least Concern" ~ "NT",
                     iucn == "Vulnerable" ~ "T",
                     iucn == "Near Threatened"   ~ "NT",
                     iucn == "Endangered"   ~ "T",
                     iucn == "Critically Endangered"~ "T"),
    branching = factor(branching, levels=c("NB", "LB", "MB", "HB")),
    iucn = ifelse(iucn == "T", 1, 0) %>% as.factor())

# Set up h2o session ------------------------------------------------------

# initialize h2o session
h2o.init(max_mem_size = "8g", nthreads=6)

# convert to h2o object
df.h2o <- as.h2o(df)

seed.no <- 42 #for reproducibility

# Split h2o data into training, validation, and test frames
data.split <- h2o.splitFrame(df.h2o,ratios = 0.75, seed=seed.no)

train <- data.split[[1]]  # For training
test <- data.split[[2]] # For final evaluation of model performance


# ROC --------------------------------------------------------------

auc_val <- function(type, w, dataset="train"){
  
  hours <- 8
  folder <- paste0("model_",type,"_norm_", hours, "h")
  folder_name <- folder
  
  leaderboard <- read.csv(here("output", paste0("Table_S_maximised_threshold_model_", type, "_norm_", hours, "h.csv")), stringsAsFactors = FALSE) %>% arrange(desc(AUC.test))
  auc.df <- read.csv(here("output", paste0("Table_S_model_performance_",type,".csv"))) %>% filter(cutoff > 0.32)
  win <- auc.df$n[w]
  
  aml_leader <- h2o.loadModel(here("output", folder_name, leaderboard$fname[win]))
  
  res.train <- h2o.predict(aml_leader , train)
  res.train <- as.data.frame(res.train)$p1
  
  res.test <- h2o.predict(aml_leader , test)
  res.test <- as.data.frame(res.test)$p1
  
  #all values lower than cutoff value will be classified as 0 (NT in this case)
  train.labels <- ifelse(res.train < leaderboard$cutoff[win], 0, 1)
  test.labels <- ifelse(res.test < leaderboard$cutoff[win], 0, 1)
  true.train.labels <- as.data.frame(train)$iucn
  true.test.labels <- as.data.frame(test)$iucn
  
  if (dataset=="train"){
    auc <- roc(true.train.labels, train.labels)  
    roc.p <- wilcox.test(train.labels[true.train.labels == 1], 
                         train.labels[true.train.labels == 0], 
                         alternative = "great")
  } else {
    auc <- roc(true.test.labels, test.labels)
    roc.p <- wilcox.test(test.labels[true.test.labels == 1], 
                         test.labels[true.test.labels == 0], 
                         alternative = "great")
  }
  
  return(list(auc, roc.p))
}

type <- combn(c("binary","morph", "dist"),2)
win <- c(binary=1, morph=2, dist=3)

aucCompare <- list()

for(d in c("train", "test")){
  res <- matrix(ncol=6, nrow=ncol(type))
  
  for (t in 1:2){
    temp <- type[,t]
    
    roc1 <- auc_val(temp[1],win[temp[1]],d)
    roc2 <- auc_val(temp[2],win[temp[2]],d)
    
    rocTest <- roc.test(roc1[[1]], roc2[[1]])
    
    res[t,] <- c(d, temp, rocTest$method, rocTest$statistic, rocTest$p.value)
  }
  
  colnames(res) <- c("dataset", "roc1", "roc2", "method", "statistic", "pvalue")
  aucCompare[[d]] <- as.data.frame(res)
}

write.csv(do.call(rbind, aucCompare), 
          here("output", "Table_S_auc_comparison.csv"), row.names=FALSE)

