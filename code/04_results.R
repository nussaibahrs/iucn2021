# Set up environment ------------------------------------------------------
# helper packages
library(here)
library(tidyverse)
library(magrittr)

# setting up machine learning models
library(h2o)

# packages for explaining our ML models
library(pROC)
library(hmeasure)
library(pdp)
library(DALEX)
library(ingredients)

# plotting packages
library(icosa)
library(ggthemes)
library(patchwork)
library(grid)
library(ggalluvial)

#load functions
source(here("code", "functions.R"))

#rescaling
scaled.new <- function (newdata, olddata) {scale(newdata, attr(olddata, "scaled:center"), attr(olddata, "scaled:scale"))}
unscale <- function(newdata, olddata){
  newdata * attr(olddata, 'scaled:scale') + attr(olddata, 'scaled:center')
}

# Data  -------------------------------------------------------------------
#palette
u_col <- ggsci::pal_uchicago("default")(9)[c(2,5,4,3,1)]

# classification data
df.corals <- read.csv(here("data", "traits_iucn.csv"), stringsAsFactors = FALSE) %>% mutate(
  range = ifelse(range > 0, range, NA),
  corallite_ori = corallite,
  depth_ori = max_depth,
  range_ori = range)

df.corals <- df.corals %>% mutate_at(c("max_depth", "corallite", "range"), function (x) as.numeric(scale(x))) #normalise numerical predictors

df <- df.corals %>%
  na.omit() %>%
  
  #omit species name
  select(-sp) %>%
  
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
data.split <- h2o.splitFrame(df.h2o,ratios = c(0.6, 0.2), seed=seed.no)

train <- data.split[[1]]  # For training
test <- data.split[[2]] # For final evaluation of model performance

# variable names for response & features
y <- "iucn"
x <- c("max_depth", "branching", "corallite", "range" )

# Load Model --------------------------------------------------------------

# * Load 10 best models in h2o session ------------------------------------

folder_name <- "model"
leaderboard <- read.csv(here("output", folder_name, "leaderboard.csv"), stringsAsFactors = FALSE)

mods <- list() 

for(i in 1:nrow(leaderboard)){
  mods[[i]] <- h2o.loadModel(here("output", folder_name, leaderboard$model_id[i]))
}

# * Find probability threshold  -----------------
cutoff=seq(0, 1, 0.01) #threshold values
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

#maximisation based on youden index
temp <- opt_par %>% group_by(fname) %>% summarise(youden.train=max(youden.train), youden.test=max(youden.test)) %>%
  na.omit()

opt_train <- temp %>% select(fname, youden.train) %>% 
  left_join(opt_par %>% select(fname, AUC.train, youden.train, cutoff.train), by=c("fname", "youden.train"))

opt_test <- temp %>% select(fname, youden.test) %>% 
  left_join(opt_par %>% select(fname, AUC.test, youden.test, cutoff.test), by=c("fname", "youden.test"))

opt_par_final <- full_join(opt_train, opt_test) %T>%
  write.csv(here("output", "Table_S_maximised_threshold_models.csv"), row.names = FALSE) %>%
  group_by(fname) %>%
  mutate(cutoff = (cutoff.test+ cutoff.train)/2) %>%
  arrange(desc(AUC.train), desc(AUC.test)) %>%
  filter(cutoff > 0.3)

opt_par_final
win <- 1
aml_leader <-  h2o.getModel(opt_par_final$fname[win])


# * Confusion Matrix ------------------------------------------------------
res.train <- h2o.predict(aml_leader , train)
res.train <- as.data.frame(res.train)$p1

res.test <- h2o.predict(aml_leader , test)
res.test <- as.data.frame(res.test)$p1

#all values lower than cutoff value will be classified as 0 (NT in this case)
train.labels <- ifelse(res.train < opt_par_final$cutoff[win], 0, 1)
test.labels <- ifelse(res.test < opt_par_final$cutoff[win], 0, 1)
true.train.labels <- as.data.frame(train)$iucn
true.test.labels <- as.data.frame(test)$iucn

cf.train <- caret::confusionMatrix(as.factor(train.labels), true.train.labels)
cf.test <- caret::confusionMatrix(as.factor(test.labels), true.test.labels)

prop.table(cf.train$table, 2)
prop.table(cf.test$table, 2)


# * Misclassification -----------------------------------------------------
train.df <- train %>% as.data.frame %>% 
  mutate(pred = ifelse(train.labels==0, "NT", "T"))

test.df <- test %>% as.data.frame %>% 
  mutate(pred = ifelse(test.labels==0, "NT", "T"))


df <- bind_rows(train.df, test.df) %>%
  mutate(
    status = factor(status, levels=c("Least Concern", "Near Threatened", "Vulnerable", "Endangered", "Critically Endangered")),
         iucn = case_when(status == "Least Concern" ~ "NT",         
                          status == "Vulnerable" ~ "T",           
                          status == "Near Threatened"   ~ "NT",    
                          status == "Endangered"   ~ "T",
                          status == "Critically Endangered"   ~ "T")) 
# per iucn status
df %>% group_by(status, iucn, pred) %>%
  tally() %>%
  group_by(status) %>%
  mutate(tot=sum(n), 
         prop = n/sum(n)) %>%
  ungroup %>%
  reshape2::dcast(status+iucn~pred, value.var = "prop")

# per branching
df %>% 
  group_by(branching, iucn, pred) %>%
  tally() %>%
  group_by(branching) %>%
  mutate(tot=sum(n),
         prop = n/tot) %>%
  filter(iucn!=pred) %>%
  summarise(prop =sum(prop))

# per corallite diameter
df %>% 
  mutate( ints = cut(corallite_ori ,breaks = c(0, 5, 250))) %>%
  group_by(iucn, pred, ints) %>%
  tally() %>%
  ungroup() %>%
  group_by(ints) %>%
  mutate(prop = n/sum(n)) %>%
  filter(iucn != pred) %>%
  group_by(ints) %>% summarise(sum(prop))

# per range
df %>% 
  mutate( ints = cut(range_ori ,breaks=c(0, 1000, 20000))) %>%
  group_by(iucn, pred, ints) %>%
  tally() %>%
  ungroup() %>%
  group_by(ints) %>%
  mutate(prop = n/sum(n)) %>%
  filter(iucn != pred) %>%
  group_by(ints) %>% summarise(sum(prop))

# per range
df %>% 
  mutate( ints = cut(depth_ori ,breaks=c(0, 30, 165))) %>%
  group_by(iucn, pred, ints) %>%
  tally() %>%
  ungroup() %>%
  group_by(ints) %>%
  mutate(prop = n/sum(n)) %>%
  filter(iucn != pred) %>%
  group_by(ints) %>% summarise(sum(prop))

rm(list=("mods")) #remove models to save memory

# * Model agnostic metrics

# create a data frame with just the features
features <- as.data.frame(df) %>% select(-iucn)
#  Create a numeric vector with the actual responses
response <- ifelse(df$iucn == "NT", 0, 1)

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

# *  Partial dependence profile ---------------------------------------------------

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
pdp_plot <- function (pdp, olddata=NULL, levels=NULL, unscale=TRUE, col="darkred"){
  lab <- levels(pdp$`_vname_`)
  pdp <- as.data.frame(pdp)
  
  if (unscale == TRUE) {pdp$`_x_` <- unscale(pdp$`_x_`, olddata)}
  
  if(is.factor(pdp$`_x_`)){
    pdp$`_x_` <- factor(as.character(pdp$`_x_`), levels=levels)
    
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
p1 <- pdp_plot(pdp_corallite, olddata=scale(df$corallite_ori), col = u_col[2]) + labs(x="Corallite Diameter (mm)") + 
  coord_cartesian(xlim=c(0, 25)) + geom_hline(yintercept = opt_par_final$cutoff[win], linetype="dashed", col="darkgrey")

p2 <- pdp_plot(pdp_branching, unscale = FALSE, levels=c("NB", "LB", "MB", "HB")) + scale_fill_manual(values=u_col[2:5]) + 
  labs(x="Degree of branching")  + geom_hline(yintercept = opt_par_final$cutoff[win], linetype="dashed", col="darkgrey")

p3 <- pdp_plot(pdp_range, olddata=scale(df$range_ori), col = u_col[3]) + labs(x="Range (normalised)") + 
  geom_hline(yintercept = opt_par_final$cutoff[win], linetype="dashed", col="darkgrey")

p4 <- pdp_plot(pdp_max_depth, olddata=scale(df$depth_ori), col = u_col[4]) + labs(x="Maximum depth (m)") +
  coord_cartesian(xlim=c(0, 100)) + geom_hline(yintercept = opt_par_final$cutoff[win], linetype="dashed", col="darkgrey")

svg(here("figs", "Fig_03_model_metrics.svg"), w=8, h=8)
p1 + p2 + p3 + p4  + plot_annotation(tag_levels = "A")
dev.off()

# DD species ------------------------------------------------
dd <- df.corals %>% filter(iucn =="Data Deficient") %>% select (-iucn)

# * Prediction of threat status -------------------------------------------

res <- h2o.predict(aml_leader, as.h2o(dd))
res <- as.data.frame(res)$p1

# calculate threat status based on threshold cutoff
dd$status <- ifelse(res < opt_par_final$cutoff[win], "NT", "T")

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
             aes(x=genus, y=n, col=status2), cex=7.5, inherit.aes = FALSE) +
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
  filter(!between(decimalLatitude, 39, 45) & !between(decimalLongitude, -38, 30)) %>% #remove mediterranean ones
  filter(scientificName %in% dd$sp) %>% #select only occurrences for dd species
  distinct(scientificName, decimalLatitude, decimalLongitude) %>%
  left_join(dd %>% dplyr::select(sp, status), by=c("scientificName" = "sp"))

gr <- hexagrid(16, sp=TRUE) #2.5 degrees
dd_occ$cell <- locate(gr, dd_occ[, c(3,2)])

me <- tapply(INDEX=dd_occ$cell, X=dd_occ$status, function (x) length(x[x == "T"])/length(x[x=="NT"| x == "T"]))
mean(me, na.rm = TRUE)

world <- rworldmap::getMap()
cuts <- seq(0.0, 1, 0.2)
map_col <- colorRampPalette(u_col[c(2,4,5)])(length(cuts)-1)

svg(here("figs", "Fig_05_map_dd_risk.svg"), w=8, h=5.4)
plot(gr, border="NA")

for (i in 1:(length(cuts)-1)){
  n <- which(me >= cuts[i] & me < cuts[i+1])  
  temp <- names(n)
  
  plot(gr[temp], col=map_col[i], border="white", add=TRUE)
}

plot(world, col=scales::alpha("grey95", 0.8), border="grey80", add=TRUE)
cuts_labs <- cuts*100 #turn to percentage
legend("bottomleft", paste0(cuts_labs[1:(length(cuts)-1)], "-", cuts_labs[2:length(cuts)], "%"), fill=map_col, horiz=TRUE, border=NA, bty="n", 
       title = expression(bold("Percentage of DD species \nunder threat")))
dev.off()


# Plio-Pleistocene --------------------------------------------------------

pleist.df <- read.csv(here("data", "pleist_resolved.csv"), stringsAsFactors = FALSE) %>%
  mutate(corallite = as.numeric(gsub(">", "", corallite)))

pleist.df <- pleist.df %>% group_by(valid_name) %>%
  summarise(corallite = max(corallite), max_depth=max(max_depth), range = max(range), #prop_range = max(prop_range), 
            branching = unique(branching)[1], globally.extinct = unique(globally.extinct), 
            regionally.extinct = unique(regionally.extinct)) 

#scale based on original data
pleist.df <- pleist.df %>% 
  mutate(corallite = as.numeric(scaled.new(corallite, scale(df.corals$corallite_ori))),
         max_depth = as.numeric(scaled.new(max_depth, scale(df.corals$depth_ori))),
         range = as.numeric(scaled.new(max_depth, scale(df.corals$range_ori)))
         )

pleist.df[pleist.df$branching == "",]$branching <- NA 

pleist.df[pleist.df$globally.extinct == "extinct" & pleist.df$regionally.extinct == "extant",]$regionally.extinct <- "extinct"
res <- h2o.predict(aml_leader, as.h2o(pleist.df))
res <- as.data.frame(res)$p1

# calculate threat status based on threshold cutoff
pleist.df$status <- ifelse(res < opt_par_final$cutoff[win], "NT", "T")
#number of missing values for each species
pleist.df$na_count <- apply(pleist.df, 1, function(x) sum(is.na(x)))
pleist.df[pleist.df$na_count >2,]$status <- "DD"

prop.table(table(pleist.df$status[pleist.df$status != "DD"]))

pleist.summ <- pleist.df %>% 
  filter(status != "DD") %>%
  mutate(status = case_when(status == "T" ~ "Threatened", 
                            status == "NT" ~"Not Threatened"),
         globally.extinct = tools::toTitleCase(globally.extinct),
         regionally.extinct = tools::toTitleCase(regionally.extinct)) %>%
  group_by(globally.extinct, regionally.extinct, status) %>%
  tally() %>%
  ungroup()

# % of threaned status going extinct
pleist.summ %>% group_by(status, globally.extinct) %>% summarise(n=sum(n)) %>%
  ungroup() %>%
  left_join(pleist.summ %>% group_by(globally.extinct) %>% summarise(tot=sum(n))) %>%
  mutate(prop = n/tot)

pleist.summ %>% group_by(status, regionally.extinct) %>% summarise(n=sum(n)) %>%
  group_by(regionally.extinct) %>%
  mutate(prop=n/sum(n))

 
ggplot(pleist.summ,
       aes(y = n,axis1 = globally.extinct,
           axis2 = regionally.extinct,  axis3 = status)) +
  scale_fill_manual(values = u_col[c(2,5)]) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", aes(fill = status), width=0.1) +
  geom_stratum(width = 0.1) +
  geom_text(stat = "stratum", infer.label = TRUE, reverse = TRUE, angle=90) +
  scale_x_continuous(breaks = 1:3, labels = c("Globally", "Regionally", "Predicted \nStatus"), expand = expand_scale(mult = 0)) +
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

ggsave(here("figs", "Fig_06_pleist_predictions.svg"), w=6, h=8)
