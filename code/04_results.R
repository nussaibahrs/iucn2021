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

# * Load 20 best models in h2o session ------------------------------------

folder_name <- "model"
leaderboard <- read.csv(here("output", folder_name, "leaderboard.csv"), stringsAsFactors = FALSE)

mods <- list() 

for(i in 1:nrow(leaderboard)){
  mods[[i]] <- h2o.loadModel(here("output", folder_name, leaderboard$model_id[i]))
}

# * Find probability threshold by maximising youden index -----------------
cutoff=seq(0, 1, 0.05) #threshold values
n=20 #number of models

opt_par <- as.data.frame(matrix(0, nrow=n, ncol=7))
colnames(opt_par) <- c("fname", "model", "cutoff", "youden", 
                       "AUC", "misclass_rate", "h")
for (i in 1:n){
  opt_par$model[i] <-mods[[i]]@algorithm
  
  for (c in 2:length(cutoff)){
    
    res <- h2o.predict(mods[[i]], test)
    res <- as.data.frame(res)$p1
    
    #all values lower than cutoff value will be classified as 0 (NT in this case)
    test.labels <- ifelse(res < cutoff[c], 0, 1)
    true.labels <- as.data.frame(test)$iucn
    
    met <- hmeasure::HMeasure(true.labels,test.labels)   
    
    if (opt_par[i,]$youden < met$metrics$Youden) {
      opt_par[i,]$fname <- mods[[i]]@model_id
      opt_par[i,]$cutoff <- cutoff[c]
      
      #error metrics
      mod.roc <- pROC::roc(true.labels, test.labels) 
      
      #performance
      opt_par[i, c("cutoff", "AUC", 
                   "misclass_rate", "youden", "h")] <- c(cutoff[c],
                                                         pROC::auc(mod.roc), #AUC
                                                         met$metrics$ER, #misclassification rate = 1 - Accuracy
                                                         met$metrics$Youden, # Youden Index
                                                         met$metrics$H
                   )
    }
  }
}

#get final model
opt_par <- opt_par %>% arrange(misclass_rate, desc(AUC))
opt_par %>% select(model, cutoff, AUC, misclass_rate, youden) %>%
  mutate_if(is.numeric, round, digits=3) %>%
  mutate(misclass_rate = 1-misclass_rate) %>%
  set_names(c("model type", "threshold", "AUC", "Accuracy", "Youden Index")) %T>%
  write.csv(here("output", "Table_01_AutoML_results.csv"), row.names = FALSE)

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

p_vi <- aml_varimp %>% 
  as.data.frame() %>%
  filter(relative_importance > 0) %>%
  mutate(rank = nrow(.):1) %>%
  ggplot(aes(x=reorder(variable, rank), y=scaled_importance)) +
  geom_bar(stat="identity", width = 0.6, fill = u_col[5]) +
  scale_x_discrete(labels=rev(c("corallite diameter", "branching: LB", "branching: NB", "branching: HB", "range", "branching: MB", "maximum water depth"))) +
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
dd$status <- ifelse(res < opt_par$cutoff[1], "T", "NT")

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

mround <- function(x,base){
  base*round(x/base)
}

binwidth <- 5
  
dd_occ[,c(2,3)] <- apply(dd_occ[,c(2,3)], 2, mround, base=binwidth/2)

hexes <- hex_coord_df(dd_occ$decimalLongitude, dd_occ$decimalLatitude, width=5, height=4)

bin = hex_bin(dd_occ$decimalLongitude, dd_occ$decimalLatitude, var4=dd_occ$status, width=5, height=4)
hexes = hex_coord_df(x=bin$x, y=bin$y, width=attr(bin,"width"), height=attr(bin,"height"), size=rep(1, nrow(bin)))
hexes$prop = rep(bin$var, each=6)*100
world <- map_data("world")

bin <- hex_bin(dd_occ$decimalLongitude, dd_occ$decimalLatitude, var4=dd_occ$status, width=5, height=4, 
               func=function(x){t <- table(x); t["DD"]/sum(t)}) %>%
  mutate(dd = ifelse(var > 0.5, "DD", "DS"))
hexes$dd = rep(bin$dd, each=6)

hexes %>% left_join(hexes_dd, by="id") %>% head()

p <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  geom_polygon(data=hexes, aes(x=x, y=y, fill=prop, group=id, col=dd)) + coord_equal(expand=FALSE) +
  labs(x="Longitude", y="Latitude", fill="Percentage of \nDD species threatened") +
  scale_fill_gradient2(low=u_col[2],  mid=u_col[4], high=u_col[5], midpoint=50) +
  scale_color_manual(values=c("black", "white"), guide=FALSE)+
  theme_light(base_size = 15) +
  theme(axis.title = element_text(size = 12, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(family = "Roboto Mono", size = 10),
        axis.text = element_text(family = "Roboto Mono", size = 10),
        panel.grid = element_blank(), 
        legend.position = c(.15,0.07),
        legend.direction = "horizontal")
  

ggsave(here("figs", "Fig_05_map_dd_risk.svg"), p, width = 12, h=6)
