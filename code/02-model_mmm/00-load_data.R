library(dplyr)
library(h2o)

# Data --------------------------------------------------------------------

df.corals <- read.csv(file.path("data", "traits_iucn.csv")) 
colnames(df.corals)[2] <- "iucn"

df <- df.corals %>% 
  mutate(
    #reclassify iucn status 
    iucn = case_when(iucn == "Least Concern" ~ "LC",         
                     iucn == "Vulnerable" ~ "VU",           
                     iucn == "Near Threatened"   ~ "NT",    
                     iucn == "Endangered"   ~ "EN",         
                     iucn == "Critically Endangered"~ "CR",
                     iucn == "Data Deficient"~ "DD"),
    iucn = factor(iucn, levels = c("DD", "LC", "NT", "VU", "EN", "CR")),
    branching = factor(branching, levels=c("HB", "MB", "LB", "NB"))
  )

df$extinct <- 1
df$extinct[df$iucn == "DD"] <- -1
df$extinct[df$iucn %in% c("LC", "NT")] <- 0

# Set up h2o session ------------------------------------------------------

# initialize h2o session
h2o.init(max_mem_size = "8g", nthreads=6)

# convert to h2o object
temp <- df[df$extinct > -1,]
temp$extinct <- factor(temp$extinct, levels=c(0,1))

df.h2o <- as.h2o(temp)
df.h2o$extinct <- as.factor(df.h2o$extinct)
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
folder <- paste0("model_mmm_", hours, "h")
run_time <-60 * 60 * hours #maximum time to run automl

