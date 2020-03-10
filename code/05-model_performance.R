# Set up environment ------------------------------------------------------
# helper packages
library(here)
library(tidyverse)

# setting up machine learning models
library(h2o)

# packages for explaining our ML models
library(pROC)
library(hmeasure)
library(pdp)
library(vip)
library(DALEX)

# Data  -------------------------------------------------------------------

#palette
u_col <- ggsci::pal_uchicago("default")(9)[c(2,5,4,3,1)]

# classification data
df.corals <- read.csv(here("data", "traits_iucn.csv"), stringsAsFactors = FALSE)

df <- df.corals %>%
  na.omit() %>%
  
  #omit species name
  select(-sp) %>%
  
  #remove all Data Deficient corals
  filter(iucn != "Data Deficient") %>%
  
  #treat as binary
  mutate(
    #reclassify iucn status 
    iucn = case_when(iucn == "Least Concern" ~ "NT",         
                     iucn == "Vulnerable" ~ "T",           
                     iucn == "Near Threatened"   ~ "NT",    
                     iucn == "Endangered"   ~ "T",         
                     iucn == "Critically Endangered"~ "T"),
    branching = factor(branching, levels=c("NB", "LB", "MB", "HB")),
    iucn = ifelse(iucn == "T", 1, 0) %>% as.factor()) 

df[, c("max_depth", "corallite", "range")] <- apply(df[, c("max_depth", "corallite", "range")], 2, scale)
# Set up h2o session ------------------------------------------------------

# initialize h2o session
h2o.init(max_mem_size = "8g", nthreads=6)

# convert to h2o object
df.h2o <- as.h2o(df)set.seed(0)   # For reproducibility of train/test split

# Split h2o data into training, validation, and test frames
data.split <- h2o.splitFrame(df.h2o,ratios = .75)

train <- data.split[[1]]  # For training
test <- data.split[[2]] # For final evaluation of model performance

# variable names for response & features
y <- "iucn"
x <- setdiff(names(df), y)

# Set up h2o session ------------------------------------------------------

# initialize h2o session
h2o.init(max_mem_size = "8g", nthreads=6)

# load models
mods <- dir(here("output", "model"))
h2o_mods <- list()

for (i in 1:length(mods)){
  h2o_mods[[i]] <- h2o.loadModel(here("output", "model", mods[i]))
}


# Model performance -------------------------------------------------------

# * Model metrics after threshold optimisation ------------------------------

cutoff=seq(0, 1, 0.01) #threshold values
opt_par <- as.data.frame(matrix(0, nrow=length(mods), ncol=8))
colnames(opt_par) <- c("fname", "model", "cutoff", "youden", "bestTune", 
                       "AUC", "misclass_rate", "h")

for (i in 1:length(mods)){
  opt_par$model[i] <- h2o_mods[[i]]@algorithm
  
  for (c in 2:length(cutoff)){
    
    res <- h2o.predict(h2o_mods[[i]], test)
    res <- as.data.frame(res)$p1
    
    #all values lower than cutoff value will be classified as 0 (NT in this case)
    test.labels <- ifelse(res < cutoff[c], 0, 1)
    true.labels <- as.data.frame(test)$iucn
    
    met <- HMeasure(true.labels,test.labels)   
    
    if (opt_par[i,]$youden < met$metrics$Youden) {
      opt_par[i,]$fname <- mods[[i]]
      opt_par[i,]$cutoff <- cutoff[c]
      
      #error metrics
      mod.roc <- roc(true.labels, test.labels) 
      
      #performance
      opt_par[i, c("cutoff", "AUC", 
                            "misclass_rate", "youden", "h")] <- c(cutoff[c],
                                auc(mod.roc), #AUC
                                met$metrics$ER, #misclassification rate = 1 - Accuracy
                                met$metrics$Youden, # Youden Index
                                met$metrics$H
      )
      
      temp <- h2o_mods[[i]]@parameters
      temp <- temp[!names(temp) %in% c("model_id","nfolds","keep_cross_validation_predictions" ,"keep_cross_validation_fold_assignment",
                                       "fold_assignment","balance_classes", "overwrite_with_best_model", "stopping_rounds", "x", "y", "seed")]
      temp <- paste(names(temp), temp, sep="=")
      opt_par$bestTune[i] <- paste(temp, collapse = ", ")
    }
  }
  

}

opt_par %>% 
  mutate(model_name = paste(model, 1:15, sep="_")) %>% #change when correct names available
  arrange(desc(AUC)) %T>%
  write.csv(here("output", "opt_par.csv"), row.names = FALSE) %>%
  mutate(row=row_number(),
         misclass_rate = 1-misclass_rate) %>%
  select(row, model_name, AUC, misclass_rate) %>%
  reshape2::melt(id.vars =c("row", "model_name")) %>%
  ggplot(aes(x=reorder(model_name,row), y=value, fill=variable)) +
  scale_y_continuous(breaks=seq(0,1, 0.1)) +
  geom_bar(stat="identity", position = "dodge") +
  theme(axis.text.x = element_text(angle=45, hjust=1))

# * Model agnostic metrics
# create a data frame with just the features
features <- df %>% select(-iucn)
#  Create a numeric vector with the actual responses
response <- as.numeric(df$iucn)

# Create custom predict function that returns the predicted values as a
#    vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results[[3L]])
}

i=2 #replace by loop

# create model agnostic object
dalex_explainer <- DALEX::explain(
  model = h2o_mods[[i]],
  data = features,
  y = response,
  predict_function = pred,
  label = mods[[i]]
)

#global feature importance
dalex_vip <- ingredients::feature_importance(dalex_explainer, 
                                 loss_function = loss_root_mean_square)
plot(dalex_vip)

#plot(vip1, vip2, vip3,...)


# partial dependence
library(pdp)
library(DALEX)

pdp_rf <- ingredients::partial_dependency(dalex_explainer, variables = "corallite")

plot(pdp_rf) +
  ggtitle("Partial Dependency profile for corallite") #do for multiple ones

