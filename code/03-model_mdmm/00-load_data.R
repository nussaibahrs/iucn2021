library(h2o)

# classification data
df.corals <- read.csv(file.path("data", "traits_iucn_mdmm.csv"), 
                      stringsAsFactors = FALSE)

df.corals <- df.corals %>% 
  #standardise numeric variables between 0 and 1
  mutate_if(is.numeric, function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))})

df <- df.corals %>%
  na.omit() %>%
  
  #omit species name
  dplyr::select(-sp) %>%
  
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

# Set up h2o session ------------------------------------------------------

# initialize h2o session
h2o.init(max_mem_size = "16g", nthreads=6)

# convert to h2o object
df.h2o <- as.h2o(df)

seed.no <- 42 #for reproducibility

# Split h2o data into training, validation, and test frames
data.split <- h2o.splitFrame(df.h2o,ratios = .75, seed=seed.no)

train <- data.split[[1]]  # For training
test <- data.split[[2]] # For final evaluation of model performance

# variable names for response & features
y <- "iucn"
x <- c("corallite","branching", "budding", "family", "integration", "max_depth", "range")
labsx <- c("Corallite\nDiameter", "Degree of\nbranching", "Budding Type", "Family", "Corallite\nIntegration", "Maximum Water Depth", "Geographic Range")
# Models ------------------------------------------------------------------
hours <- 8
folder <- paste0("model_mdmm_", hours, "h")
run_time <-60 * 60 * hours #maximum time to run automl

