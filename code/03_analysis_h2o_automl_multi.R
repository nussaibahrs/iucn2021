# Set up environment ------------------------------------------------------
# setting up machine learning models
library(h2o)

# Data  -------------------------------------------------------------------
setwd("/media/adam/data/shared/Dropbox/WorkSpace/Help/2020-01-24_NussIUCN/")

# classification data
df.corals <- read.csv(file.path("data", "traits_iucn.csv"), stringsAsFactors = FALSE) 
df.corals$range[df.corals$range == 0] <- NA

df.corals$range <- as.numeric(df.corals$range)
df.corals$corallite <- as.numeric(df.corals$corallite)

df <- subset(na.omit(df.corals[,-1]), iucn != "Data Deficient")
df$branching <-  factor(df$branching, levels=c("NB", "LB", "MB", "HB"))
df$iucn <- as.factor(df$iucn)

levels(df$iucn) <- c("EN", "EN", "LC", "NT", "VU")


# Set up h2o session ------------------------------------------------------

# initialize h2o session
h2o.init(max_mem_size = "8g", nthreads=6)

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

run_time <- 60*60*8 #maximum time to run automl

# * Automatic Machine Learning ----------------------------------------------
system.time(aml <- h2o.automl(y=y, x=x, 
                              training_frame = train,
                  keep_cross_validation_fold_assignment = TRUE,
                  keep_cross_validation_predictions = TRUE,
                  stopping_metric = "mean_per_class_error",
                  balance_classes = TRUE,
                  nfolds=10,
                  max_runtime_secs = run_time,
                  max_models = 500))


# * Save models -----------------------------------------------------------

folder_name <-"model_multi"
system(paste("mkdir ", "export/", folder_name, sep=""))

leaderboard <- as.data.frame(aml@leaderboard)

for(i in 1:nrow(leaderboard)) {
  aml1 <- h2o.getModel(aml@leaderboard[i, 1]) # get model object in environment
  h2o.saveModel(object = aml1, file.path("/media/adam/data/shared/Dropbox/WorkSpace/Help/2020-01-24_NussIUCN/", "export", folder_name)) # pass that model object to h2o.saveModel as an argument
  
}

write.csv(leaderboard, file.path("export", folder_name, "leaderboard.csv"), row.names = FALSE)


# Terminate h2o session ---------------------------------------------------
h2o.shutdown()



