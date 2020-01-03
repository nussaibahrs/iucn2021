library(here)
library(tidyverse)
library(caret)
library(doMC)
library(parallel)

# configure multicore
n_cores <- detectCores()
registerDoMC(cores = 4)

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
reps <- 10
cv <- 10
metric <- "AUC"

trCtrl <- trainControl(method = "repeatedcv", 
                       number = cv, 
                       repeats = reps, 
                       classProbs = TRUE,
                       summaryFunction = multiClassSummary,
                       sampling = "up",
                       search = "grid",
                       savePredictions = T)

##### RANDOM FOREST ####
rfGrid <- expand.grid(mtry = c(1: 20))

mx_nd <- 5:15
store_maxnode <- list()

for (maxnodes in mx_nd) {
  set.seed(42)
  rf_maxnode <- train(iucn~., data=trainTransformed,
                 method = "rf",
                 metric = metric,
                 tuneGrid = rfGrid,
                 trControl = trCtrl,
                 importance = TRUE,
                 nodesize = 14,
                 maxnodes = maxnodes,
                 ntree = 300)
current_iteration <- toString(maxnodes)
store_maxnode[[current_iteration]] <- rf_maxnode
}

results_mtry <- resamples(store_maxnode)
res <- as.data.frame(summary(results_mtry)$statistics$AUC)
maxnode <- mx_nd[which(res$Mean == max(res$Mean))]

store_maxtrees <- list()
ntr <- c(250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 2000)
for (ntree in ntr) {
  set.seed(42)
  rf_maxtrees <- train(iucn~., data=trainTransformed,
                       method = "rf",
                       metric = metric,
                       tuneGrid = rfGrid,
                       trControl = trCtrl,
                       importance = TRUE,
                       nodesize = 14,
                       maxnodes = 15,
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
}

results_tree <- resamples(store_maxtrees)
res <- as.data.frame(summary(results_tree)$statistics$AUC)
ntree <- ntr[which(res$Mean == max(res$Mean))]

set.seed(42)
rfFit <- train(iucn~., data=trainTransformed,
               method = "rf",
               metric = metric,
               tuneGrid = rfGrid,
               trControl = trCtrl,
               importance = TRUE,
               nodesize = 14, #might need to think about choose nodesize as well
               maxnodes = maxnode,
               ntree = ntree)


#### SUPPORT VECTOR MACHINE ####
svmGrid <- expand.grid(sigma = seq(0,1, 0.05), 
                       C = seq(0,1, 0.05))
set.seed(42)
svmFit <- train(iucn~., data=trainTransformed,
      method = "svmRadial",
      metric = metric,
      tuneGrid = svmGrid,
      trControl = trCtrl)

#### KNN ####
knnGrid <- expand.grid(k = seq(1,50, 1))
set.seed(42)
knnFit <- train(iucn ~ ., data = trainTransformed, 
                method = "knn", 
                metric = metric,
                trControl = trCtrl, 
                tuneGrid = knnGrid)



#### Linear: Naive Bayes ####
nbGrid <- expand.grid(
  usekernel = c(TRUE, FALSE),
  fL = 0:10, #Factor for Laplace correction, default factor is 0, i.e. no correction.
  adjust = seq(0, 10, by = 1)
)
set.seed(42)
nbFit <- train(iucn ~ ., data = trainTransformed, 
               method = "nb", 
               metric = metric,
               trControl = trCtrl, 
               tuneGrid = nbGrid)


#####  Boosted Trees ####
gbmGrid <-  expand.grid(interaction.depth = 1:5, #maximum depth of trees, i.e. highest level of variable interactions allowed
                        n.trees = (5:30)*50, #number of trees
                        shrinkage = seq(0.01, 1, 0.1), #learning rate
                        n.minobsinnode = 10:50 #minimum number of observations in the terminal nodes of the trees
                        )
set.seed(42)
gbmFit <- train(iucn ~ ., data = trainTransformed, 
                method = "gbm", 
                metric = metric,
                trControl = trCtrl, 
                tuneGrid = gbmGrid)

##### NEURAL NETWORK ####
nnetGrid <-  expand.grid(size = seq(from = 1, to = 15, by = 1),
                         decay = seq(from = 0.005, to = 1, by = 0.005))

set.seed(42)
nnetFit <- train(iucn~., data=trainTransformed, method = 'nnet', 
                 metric=metric,
                 trControl = trCtrl,
                 prox=TRUE,allowParallel=TRUE, verbose=TRUE,
                 maxit = 1000,
                 tuneGrid=nnetGrid)

#### save
mods <- list(nnetFit, rfFit, svmFit, knnFit, nbFit, gbmFit)

save(preProcValues, train, test, trainTransformed, testTransformed, mods, file = here("output", "ML_multi.RData"))

  #### MODEL PERFORMANCE VS TEST DATASET ####
library(here)
library(pROC)
library(tidyverse)
library(hmeasure)
library(caret)

#load models 
load(here("output", "ML_binary.RData"))

#compute performance
mods.performance <- as.data.frame(matrix(ncol=5, nrow=length(mods)))
colnames(mods.performance) <- c("model", "AUC", "misclass_rate", "youden", "h")

for (i in 1:length(mods)){
  res <- predict(mods[[i]], newdata=testTransformed, type="prob")
  mod.roc <- roc(ifelse(testTransformed[,"iucn"] == "T", 1, 0), res[[2]])
  
  
  #error metrics
  true.labels <- relabel(as.factor(testTransformed$iucn))
  test.labels <- relabel(predict(mods[[i]], newdata=testTransformed))
  met <- HMeasure(true.labels,test.labels)               
  # plot(mod.roc)
  
  #performance
  mods.performance[i,] <- c(mods[[i]]$method, 
                            mod.roc$auc, #AUC
                           met$metrics$ER, #misclassification rate = 1 - Accuracy
                          met$metrics$Youden, # Youden Index
                          met$metrics$H
                          )
}


mods.performance %>% mutate(n=1:6) %>%
  arrange(desc(AUC))

########
# Comparing Multiple Models
# Having set the same seed before running gbm.tune and xgb.tune
# we have generated paired samples and are in a position to compare models 
# using a resampling technique.
# (See Hothorn at al, "The design and analysis of benchmark experiments
# -Journal of Computational and Graphical Statistics (2005) vol 14 (3) 
# pp 675-699) 

#needs the seeds to be the same
rValues <- resamples(mods, modelNames = mods.performance$model)
rValues$values

summary(rValues)

bwplot(rValues,metric="AUC")	# boxplot
dotplot(rValues,metric="AUC")	# dotplot
#splom(rValues,metric="ROC")

###

