library(dplyr)

# Data --------------------------------------------------------------------

# classification data
df.corals <- read.csv(file.path("data", "pbdb_pleist.csv"), stringsAsFactors = FALSE) 
df.corals[df.corals == ""] <- NA


df <- df.corals %>%
  #na.omit() %>%
  #treat as binary
  mutate(
    #reclassify iucn status 
    extinct = case_when(extinct == "extant" ~ 0,         
                        TRUE~1),
    extinct=factor(extinct, levels = c(0,1)),
    branching = factor(branching, levels=c("NB", "LB", "MB", "HB"))) %>% 
  #omit species name
  dplyr::select(corallite, branching, budding, family, integration, extinct) %>% 
  ungroup()

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
y <- "extinct"
x <- c("corallite","branching", "budding", "family", "integration")
labsx <- c("Corallite\nDiameter", "Degree of\nbranching", "Budding Type", "Family", "Corallite\nIntegration")
# Models ------------------------------------------------------------------
hours <- 8
folder <- folder_name <- paste0("model_fmm_", hours, "h")
run_time <-60 * 60 * hours #maximum time to run automl