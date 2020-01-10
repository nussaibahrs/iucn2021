# Set up environment ------------------------------------------------------
# helper packages
library(here)
library(tidyverse)

# setting up machine learning models
library(h2o)

# packages for explaining our ML models
library(caret)
library(pdp)
library(vip)
library(iml)
library(DALEX)
library(lime)

# initialize h2o session
h2o.init()

# Data  -------------------------------------------------------------------

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

# convert to h2o object
df.h2o <- as.h2o(df)
df.h2o[, c("max_depth", "corallite", "range")] <- scale(df.h2o[, c("max_depth", "corallite", "range")])

set.seed(0)   # For reproducibility of train/test split

# Split h2o data into training, validation, and test frames
data.split <- h2o.splitFrame(df.h2o,ratios = c(.7,.2))

train <- data.split[[1]]  # For training
valid <- data.split[[2]]  # For validating trained models and comparing different hyperparameter vectors
test <- data.split[[3]] # For final evaluation of model performance

# variable names for response & features
y <- "iucn"
x <- setdiff(names(df), y)

# define parameter search criteria grid
search.criteria <- list(
  max_models = 100,
  strategy = "RandomDiscrete",
  max_runtime_secs = 900,
  stopping_rounds = 10,
  stopping_tolerance = 0.005,
  stopping_metric = "AUC",
  seed=seed.no
)


# Models ------------------------------------------------------------------
seed.no <- 42

#### glm
# specify alpha values
alphas <- seq(0, 1, 0.1)

# * GLM -------------------------------------------------------------------

# Parameter tuning:
glm.tune.grid <- h2o.grid(
  algorithm = "glm",
  hyper_params = list(alpha = alphas),
  y = y,
  x = x,
  grid_id = "glm.tune.grid",
  stopping_metric = "AUC", 
  training_frame = train,
  validation_frame = valid,
  nfolds = 10,
  fold_assignment = "Modulo",
  lambda_search = TRUE,
  family = "binomial",
  seed=seed.no,
  keep_cross_validation_predictions = TRUE
)

glm_perf <- h2o.getGrid(grid_id = "glm.tune.grid", 
            sort_by = "auc", 
            decreasing = FALSE)

# Best model chosen by validation error: 
glm <- h2o.getModel(glm_perf@model_ids[[1]])

# Best model chosen by validation error: 
caret::confusionMatrix(h2o.predict(best_model, test) %>% as.data.frame() %>% pull(predict), 
                as.data.frame(test)$iucn, positive = "1")



# * GBM ---------------------------------------------------------
# Tuned parameters:
# define hyperparameter grid
gbm.grid <- list(
  learn_rate = seq(.03,.05,.01),
  max_depth = c(8:10),
  sample_rate = seq(0.7, 1.0, 0.1),
  col_sample_rate = seq(0.4, 1.0, 0.1)
)

# search parameter grid
gbm.tune.grid <- h2o.grid(
  algorithm = "gbm",
  grid_id = "gbm.tune.grid",
  y = y,
  x = x,
  hyper_params = gbm.grid,
  training_frame = train,
  validation_frame = valid,
  score_each_iteration = TRUE,
  nfolds = 10,
  fold_assignment = "Modulo",
  keep_cross_validation_fold_assignment = TRUE,
  keep_cross_validation_predictions = TRUE,
  search_criteria = search.criteria
)

gbm_perf <- h2o.getGrid(grid_id = "gbm.tune.grid", 
                        sort_by = "auc", 
                        decreasing = FALSE)

# Best model chosen by validation error: 
gbm <- h2o.getModel(gbm_perf@model_ids[[1]])

# Best model chosen by validation error: 
caret::confusionMatrix(h2o.predict(best_model, test) %>% as.data.frame() %>% pull(predict), 
                       as.data.frame(test)$iucn, positive = "1")


# * Random Forest ---------------------------------------------------------

rf.grid <- list(ntrees = seq(50, 500, by = 50),
     #mtries = seq(3, 5, by = 1),
     # max_depth = seq(10, 30, by = 10),
     # min_rows = seq(1, 3, by = 1),
     # nbins = seq(20, 30, by = 10),
     sample_rate = c(0.55, 0.632, 0.75))

# search parameter grid
rf.tune.grid <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf.tune.grid",
  y = y,
  x = x,
  hyper_params = rf.grid,
  training_frame = train,
  validation_frame = valid,
  score_each_iteration = TRUE,
  nfolds = 10,
  fold_assignment = "Modulo",
  keep_cross_validation_fold_assignment = TRUE,
  keep_cross_validation_predictions = TRUE,
  search_criteria = search.criteria,
  binomial_double_trees = TRUE,
  balance_classes = TRUE
)



rf_perf <- h2o.getGrid(grid_id = "rf.tune.grid", 
                        sort_by = "auc", 
                        decreasing = FALSE)

# Best model chosen by validation error: 
rf <- h2o.getModel(rf_perf@model_ids[[1]])

h2o.saveModel(rf, path = "model")
# Deep learning neural network --------------------------------------------
hyper_params <- list(
  activation = c("Rectifier", "Maxout", "Tanh", "RectifierWithDropout", "MaxoutWithDropout", "TanhWithDropout"), 
  hidden = list(c(5, 5, 5, 5, 5), c(10, 10, 10, 10), c(50, 50, 50), c(100, 100, 100)),
  epochs = c(50, 100, 200),
  l1 = c(0, 0.00001, 0.0001), 
  l2 = c(0, 0.00001, 0.0001),
  rate = c(0, 01, 0.005, 0.001),
  rate_annealing = c(1e-8, 1e-7, 1e-6),
  rho = c(0.9, 0.95, 0.99, 0.999),
  epsilon = c(1e-10, 1e-8, 1e-6, 1e-4),
  momentum_start = c(0, 0.5),
  momentum_stable = c(0.99, 0.5, 0),
  input_dropout_ratio = c(0, 0.1, 0.2),
  max_w2 = c(10, 100, 1000, 3.4028235e+38)
)

dl_grid <- h2o.grid(algorithm = "deeplearning", 
                    x = x,
                    y = y,
                    grid_id = "dl_grid",
                    training_frame = train,
                    validation_frame = valid,
                    nfolds = 10,                           
                    fold_assignment = "Modulo",
                    hyper_params = hyper_params,
                    search_criteria = search_criteria
)


# * SVM - binary only -----------------------------------------------------

h2o.psvm()


# * k-means  ------------------------------------------------------------------

h2o.kmeans(
# * Ensemble --------------------------------------------------------------

# ensemble
ensemble <- h2o.stackedEnsemble(
  x = x,
  y = y,
  training_frame = train,
  validation_frame = valid,
  metalearner_nfolds = 5,
  model_id = "ensemble",
  base_models = list(glm, rf, gbm),
  metalearner_algorithm = "glm"
)


# Model performance -------------------------------------------------------

# model performance
h2o.auc(glm, xval = TRUE)
## [1] 0.8363927
h2o.auc(rf, xval = TRUE)
## [1] 0.7882236
h2o.auc(gbm, xval = TRUE)
## [1] 0.810503
h2o.auc(ensemble, xval = TRUE)
## [1] 0.8285921


# Function for collecting cross-validation results: 

results_cross_validation <- function(h2o_model) {
  h2o_model@model$cross_validation_metrics_summary %>% 
    as.data.frame() %>% 
    select(-mean, -sd) %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate_all(as.character) %>% 
    mutate_all(as.numeric) %>% 
    select(Accuracy = accuracy, 
           AUC = auc, 
           Precision = precision, 
           Specificity = specificity, 
           Recall = recall, 
           Logloss = logloss) %>% 
    return()
}

# Model performance based on test data: 
pred_class <- h2o.predict(default_rf, test) %>% as.data.frame() %>% pull(predict)

# Function for calculating confusion matrix: 

my_cm_com_rf <- function(thre) {
  du_bao_prob <- h2o.predict(default_rf, test) %>% as.data.frame() %>% pull(p1)
  du_bao <- case_when(du_bao_prob >= thre ~ "1", 
                      du_bao_prob < thre ~ "0") %>% as.factor()
  cm <- caret::confusionMatrix(du_bao, as.factor(as.data.frame(test)$iucn), positive = "1")
  return(cm)
  
}

# Set a range of threshold for classification: 
my_threshold <- c(0.10, 0.15, 0.35, 0.5)
results_list_rf <- lapply(my_threshold, my_cm_com_rf)


# Model Agnostic Procedures -----------------------------------------------
# 1. create a data frame with just the features
features <- as.data.frame(df) %>% select(-iucn)
# 2. Create a numeric vector with the actual responses
response <- df$iucn

# 3. Create custom predict function that returns the predicted values as a
#    vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results[[3L]])
}

# example of prediction output
pred(gbm, features) %>% head()

# GBM predictor object
iml_predictor_gbm <- Predictor$new(
  model = gbm, 
  data = features,
  y = response,
  predict.fun = pred,
  class = "classification"
)
# GBM explainer
dalex_explainer_gbm <- DALEX::explain(
  model = gbm,
  data = features,
  y = response,
  predict_function = pred,
  label = "gbm"
)

# Global feature importance
dalex_vip_glm <- variable_importance(dalex_explainer_glm, n_sample = -1) 
dalex_vip_rf  <- variable_importance(dalex_explainer_rf, n_sample = -1)
dalex_vip_gbm <- variable_importance(dalex_explainer_gbm, n_sample = -1)
dalex_vip_ensemble <- variable_importance(dalex_explainer_ensemble, n_sample = -1)
plot(dalex_vip_glm, dalex_vip_rf, dalex_vip_gbm, dalex_vip_ensemble, max_vars = 10)

# partial dependence

pdp_fun <- function(object, newdata) {
  # compute partial dependence 
  pd <- mean(predict(object, as.h2o(newdata))[[3L]])
  # return data frame with average predicted value
  return(as.data.frame(pd))
}
# partial dependence values
pd_df <- partial(
  ensemble, 
  pred.var = "OverTime", 
  train = features,
  pred.fun = pdp_fun
)

# partial dependence
pd_df
##   OverTime      yhat
## 1       No 0.1137813
## 2      Yes 0.3071652
# partial dependence plot
autoplot(pd_df)

# partial dependence values
partial(
  ensemble, 
  pred.var = "Age", 
  train = features,
  pred.fun = pdp_fun,
  grid.resolution = 20
) %>% 
  autoplot(rug = TRUE, train = features) + 
  ggtitle("Age")


