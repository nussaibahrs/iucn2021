# Set up environment ------------------------------------------------------
# helper packages
library(here)
library(tidyverse)
library(ggthemes)
library(patchwork)
library(magrittr)
library(grid)

# setting up machine learning models
library(h2o)

# packages for explaining our ML models
library(pROC)
library(hmeasure)
library(pdp)
library(DALEX)
library(ingredients)

source(here("code", "functions.R"))
# Data  -------------------------------------------------------------------
#palette
u_col <- ggsci::pal_uchicago("default")(9)[c(2,5,4,3,1)]

# classification data
df.corals <- read.csv(here("data", "traits_iucn.csv"), stringsAsFactors = FALSE) %>% mutate(range = ifelse(range > 0, range, NA))
df.corals <- df.corals %>% mutate_if(is.numeric, function (x) as.numeric(scale(x))) #normalise numerical predictors

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

# Set up h2o session ------------------------------------------------------

# initialize h2o session
h2o.init(max_mem_size = "8g", nthreads=6)

# convert to h2o object
df.h2o <- as.h2o(df)

seed.no <- 42 #for reproducibility

# Split h2o data into training, validation, and test frames
data.split <- h2o.splitFrame(df.h2o,ratios = c(0.6, 0.2), seed=seed.no)

train <- data.split[[1]]  # For training
test <- data.split[[2]] # For final evaluation of model performance

# variable names for response & features
y <- "iucn"
x <- setdiff(names(df), y)

# Load Model --------------------------------------------------------------

# * Load 10 best models in h2o session ------------------------------------

folder_name <- "model"
leaderboard <- read.csv(here("output", folder_name, "leaderboard.csv"), stringsAsFactors = FALSE)

mods <- list() 

for(i in 1:nrow(leaderboard)){
  mods[[i]] <- h2o.loadModel(here("output", folder_name, leaderboard$model_id[i]))
}

# * Find probability threshold by maximising youden index -----------------
cutoff=seq(0, 1, 0.05) #threshold values
n=10 #number of models

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
    true.train.labels <- as.data.frame(train)$iucn
    true.test.labels <- as.data.frame(test)$iucn
    
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

temp <- opt_par %>% group_by(fname) %>% summarise(AUC.train=max(AUC.train), AUC.test=max(AUC.test)) %>%
  na.omit()

opt_train <- temp %>% select(fname, AUC.train) %>% 
  left_join(opt_par %>% select(fname, AUC.train, youden.train, cutoff.train), by=c("fname", "AUC.train"))

opt_test <- temp %>% select(fname, AUC.test) %>% 
  left_join(opt_par %>% select(fname, AUC.test, youden.test, cutoff.test), by=c("fname", "AUC.test"))

opt_par <- full_join(opt_train, opt_test) %>%
  group_by(fname) %>%
  mutate(cutoff = mean(cutoff.test, cutoff.train)) %>%
  filter(cutoff > 0.3)   %>%
  arrange(desc(AUC.train), desc(AUC.test))


aml_leader <-  h2o.getModel(opt_par$fname[1])
rm(list=("mods")) #remove models to save memory

# * Model agnostic metrics

# create a data frame with just the features
features <- as.data.frame(df) %>% select(-iucn)
#  Create a numeric vector with the actual responses
response <- as.numeric(df$iucn)

# Create custom predict function that returns the predicted values as a vector
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results[[3L]])
}


# * Global feature importance ---------------------------------------------
aml_varimp <- h2o.varimp(aml_leader) #calculate variable importance

# ** Fig 2: Variable Importance -------------------------------------------
aml_varimp <- aml_varimp %>% 
  as.data.frame() %>%
  filter(relative_importance > 0.11) %>%
  mutate(rank = nrow(.):1) %>%
  mutate(variable = gsub("\\.", ": ", variable),
         variable = gsub("_", " ", variable))

p_vi <- aml_varimp %>% 
  as.data.frame() %>%
  filter(relative_importance > 0.1) %>%
  mutate(rank = nrow(.):1) %>%
  ggplot(aes(x=reorder(variable, rank), y=scaled_importance)) +
  geom_bar(stat="identity", width = 0.6, fill = u_col[5]) +
  scale_y_continuous(expand=expand_scale(mult = c(0, .1))) +
  labs(x="Variables", y="Scaled Importance") +
  coord_flip() +
  theme_light(base_size = 15) +
  theme(axis.title = element_text(size = 12, face="bold"),
        axis.text.x = element_text(family = "Roboto Mono", size = 10),
        panel.grid = element_blank()) 

ggsave (here("figs", "Fig_02_var_imp.svg"), p_vi, h=4, w=8)

# *  partial dependence profile ---------------------------------------------------

# ** generate model agnostic object to explain response of features -------

dalex_explainer <- DALEX::explain(
  model = aml_leader,
  data = features,
  y = as.factor(response),
  predict_function = pred
)

# ** generate the partial dependence profiles for each variable ------------
pdp_corallite <- ingredients::partial_dependency(dalex_explainer, variables = "corallite")
pdp_branching <- ingredients::partial_dependency(dalex_explainer, variables = "branching", variable_type="categorical")
pdp_range <- ingredients::partial_dependency(dalex_explainer, variables = "range")
pdp_max_depth <- ingredients::partial_dependency(dalex_explainer, variables = "max_depth")


# ** Fig. 3: Partial dependency Plots -------------------------------------

# function to plot pdp
pdp_plot <- function (pdp, col="darkred"){
  lab <- levels(pdp$`_vname_`)
  pdp <- as.data.frame(pdp)
  
  if(is.factor(pdp$`_x_`)){
    p <- ggplot(pdp, aes(x=`_x_`, y=`_yhat_`)) + geom_bar(stat="identity", aes(fill=`_x_`), width = 0.6) + 
      labs(x=lab, y="Predicted Value") +
      scale_y_continuous(expand=expand_scale(mult = c(0, .1))) +
      theme_light(base_size = 15) +
      theme(axis.title = element_text(size = 12, face="bold"),
            axis.text.x = element_text(family = "Roboto Mono", size = 10),
            panel.grid = element_blank(), 
            legend.position = "none") 
  } else {
    p <- ggplot(pdp, aes(x=`_x_`, y=`_yhat_`)) + geom_line(col=col) + 
      labs(x=lab, y="Predicted Value") +
      theme_light(base_size = 15) +
      theme(axis.title = element_text(size = 12, face="bold"),
            axis.text = element_text(family = "Roboto Mono", size = 10),
            panel.grid = element_blank(), 
            legend.position = "none") 
  }
  
  return(p)
}

# plot pdp
p1 <- pdp_plot(pdp_corallite, col = u_col[2]) + labs(x="Corallite Diameter (normalised)")
p2 <- pdp_plot(pdp_branching) + scale_fill_manual(values=u_col[2:5]) + labs(x="Degree of branching")
p3 <- pdp_plot(pdp_range, col = u_col[3]) + labs(x="Range (normalised)")
p4 <- pdp_plot(pdp_max_depth, col = u_col[4]) + labs(x="Maximum depth (normalised)")

svg(here("figs", "Fig_03_model_metrics.svg"), w=8, h=8)
p1 + p2 + p3 + p4  + plot_annotation(tag_levels = "A")
dev.off()

# DD species ------------------------------------------------
dd <- df.corals %>% filter(iucn =="Data Deficient") %>% select (-iucn)

# * Prediction of threat status -------------------------------------------

res <- h2o.predict(aml_leader, as.h2o(dd))
res <- as.data.frame(res)$p1

# calculate threat status based on threshold cutoff
dd$status <- ifelse(res < opt_par$cutoff[1], "NT", "T")

#number of missing values for each species
dd$na_count <- apply(dd, 1, function(x) sum(is.na(x)))
dd[dd$na_count >2,]$status <- "DD"

#get genus for each species
dd$genus <- gsub( " .*$", "", dd$sp)

# ** Fig. 4: Threat status by genera (no. of sp > 2) -----------------------
dd_plot <- dd %>% filter(status != "DD") %>%
  mutate(na_count = ifelse(na_count == 0, "complete", "incomplete")) %>%
  group_by(genus, status, na_count) %>%
  tally() %>%
  group_by(genus) %>%
  mutate(n2 = sum(n),
         n3 = n/n2) %>%
  arrange(desc(n2)) %>%
  filter(n2 > 1) %>% 
  mutate(n3 =ifelse(status == "NT", -n3, n3),
         status2 = paste0(status, substring(na_count, 1, 1))) %>%
  mutate(status2 = ifelse(status2 %in% c("NTi", "Ti"), "DD", status2), 
         status2 = factor(status2, levels=c("NTc", "Tc","DD"))) %>%
  arrange(genus, status2)

dd_sum <- dd_plot %>% ungroup() %>% group_by(genus, status) %>% mutate(tot=sum(n), n=sum(n3)) %>%
  mutate(status2 = ifelse(status == "NT", "NTc", "Tc")) %>% distinct(genus, status2, n, tot)

p <- ggplot(dd_plot, aes(fill=status2, y=n3, x=reorder(genus, n2))) +
  geom_bar(stat="identity", alpha=0.5, width=0.6) +
  geom_point(data=dd_sum , 
             aes(x=genus, y=n, col=status2), cex=6, inherit.aes = FALSE) +
  geom_text(data=dd_sum, aes(x=genus, y=n, label=tot), col="white", fontface="bold", inherit.aes = FALSE) +
  scale_y_continuous(breaks = seq(-1,1, 0.5), labels = abs(seq(-1,1, 0.5))) +
  scale_x_discrete(expand=expand_scale(mult=0.05))+
  scale_color_manual(values=u_col[c(2,5)], labels=c("Not threatened", "Threatened"), guide=FALSE)+
  scale_fill_manual(values=u_col[c(2,5,1)], labels=c("Not threatened", "Threatened", "Data-deficient"))+
  labs(x="Genus", y="Proportion", fill="Status") +
  coord_flip() +
  geom_hline(yintercept = 0, col="darkgrey") +
  theme_light(base_size = 15) +
  theme(axis.title = element_text(size = 12, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(family = "Roboto Mono", size = 10),
        axis.text = element_text(family = "Roboto Mono", size = 10),
        panel.grid = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(2, 1, 1, 1), "lines"))

svg(here("figs", "Fig_04_dd_species.svg"), w=6, h=6)
p 
grid::grid.text("Threatened", x = unit(0.45, "npc"), y = unit(.95, "npc"), 
                gp=gpar(family="Roboto Mono", fontface="bold"), hjust=-1.04)
grid::grid.text("Not Threatened", x = unit(0.45, "npc"), y = unit(.95, "npc"), 
                gp=gpar(family="Roboto Mono", fontface="bold"), hjust=0.3)
dev.off()

# * Map of DD distribution risk -------------------------------------------
dd_occ <- read.csv(here("data", "2019-11-05_obis_scleractinia.csv"), stringsAsFactors = FALSE) %>%
  filter(scientificName %in% dd$sp) %>% #select only occurrences for dd species
  distinct(scientificName, decimalLatitude, decimalLongitude) %>%
  left_join(dd %>% dplyr::select(sp, status), by=c("scientificName" = "sp"))

#need script for obis cleaning here

library(icosa)

gr <- hexagrid(12, sp=TRUE)
gr

dd_occ$cell <- locate(gr, dd_occ[, c(3,2)])

me <- tapply(INDEX=dd_occ$cell, X=dd_occ$status, function (x) length(x[x == "T"])/length(x))

fl <- facelayer(gr)

fl[] <- me

world <- rworldmap::getMap()

svg(here("figs", "Fig_05_map_dd_risk.svg"), w=12, h=6)
plot(fl, col=c("blue", "yellow", "red"), legend=FALSE)
plot(world, col="lightgrey", border="darkgrey", add=TRUE)
dev.off()



# Plio-Pleistocene --------------------------------------------------------
pleist.df <- read.csv(here("data", "pleist_resolved.csv"), stringsAsFactors = FALSE) %>%
  mutate(corallite = as.numeric(gsub(">", "", corallite)), 
         corallite = as.numeric(scale(corallite)),          
         range = as.numeric(scale(prop_range)))
pleist.df[pleist.df$globally.extinct == "extinct" & pleist.df$regionally.extinct == "extant",]$regionally.extinct <- "extinct"
res <- h2o.predict(aml_leader, as.h2o(pleist.df))
res <- as.data.frame(res)$p1

# calculate threat status based on threshold cutoff
pleist.df$status <- ifelse(res < opt_par$cutoff[1], "NT", "T")

pleist.summ <- pleist.df %>% 
  mutate(status = ifelse(status == "T", "Threatened", "Not Threatened"),
         globally.extinct = tools::toTitleCase(globally.extinct),
         regionally.extinct = tools::toTitleCase(regionally.extinct)) %>%
  group_by(globally.extinct, regionally.extinct, status) %>%
  tally()


library(ggalluvial)

ggplot(pleist.summ,
       aes(y = n,
           axis1 = regionally.extinct, axis2 = globally.extinct, axis3 = status)) +
  scale_fill_manual(values = u_col[c(2,5)]) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", aes(fill = status), width=0.1) +
  geom_stratum(width = 0.1) +
  geom_text(stat = "stratum", infer.label = TRUE, reverse = TRUE, angle=90) +
  scale_x_continuous(breaks = 1:3, labels = c("Regionally", "Globally", "Predicted \nStatus"), expand = expand_scale(mult = 0)) +
  scale_y_continuous(expand = expand_scale(mult = 0)) +
  theme_light(base_size = 15) +
  labs(x="", y="No. of species", fill = "Status")+
  theme(axis.title = element_text(size = 12, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(family = "Roboto Mono", size = 10),
        axis.text = element_text(family = "Roboto Mono", size = 10, angle=45, hjust=1),
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal")

ggsave(here("figs", "Fig_06_pleist_predictions.svg"), w=6, h=6)
