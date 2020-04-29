# Set up environment ------------------------------------------------------
#setwd("/mnt/c/Users/nussa/Dropbox/nuss/iucn_traits")

# helper packages
library(here)
library(tidyverse)
library(magrittr)

# setting up machine learning models
library(h2o)

# packages for explaining our ML models
library(caret)
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
library(raster)
library(ggpubr)
library(ggrepel)

#rescaling
scaled.new <- function (newdata, olddata) {(newdata-min(olddata, na.rm = TRUE))/(max(olddata, na.rm = TRUE)-min(olddata, na.rm = TRUE))}
unscale <- function(newdata, olddata){min(olddata, na.rm = TRUE) + newdata * (max(olddata, na.rm = TRUE) - min(olddata, na.rm=TRUE))}

# Data  -------------------------------------------------------------------
#palette
u_col <- ggsci::pal_uchicago("default")(9)[c(2,5,4,3,1)]

# classification data
df.corals <- read.csv(here("data", "traits_iucn.csv"), stringsAsFactors = FALSE) %>% mutate(
  range = ifelse(range > 0, range, NA),
  corallite_ori = corallite,
  depth_ori = max_depth,
  range_ori = range)

#df.corals <- df.corals %>% mutate_at(c("max_depth", "corallite", "range"), function (x) as.numeric(scale(x))) #normalise numerical predictors
df.corals <- df.corals %>% mutate_at(c("max_depth", "corallite", "range"), function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))})

df <- df.corals %>%
  na.omit() %>%
  
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
data.split <- h2o.splitFrame(df.h2o,ratios = 0.75, seed=seed.no)

train <- data.split[[1]]  # For training
test <- data.split[[2]] # For final evaluation of model performance

# variable names for response & features
y <- "iucn"
x <- c("range", "max_depth")

# Load Model --------------------------------------------------------------

# * Load best model in h2o session ------------------------------------
hours <- 8
leaderboard <- read.csv(here("output", paste0("Table_S_maximised_threshold_model_dist_norm_", hours, "h.csv")), stringsAsFactors = FALSE) %>% arrange(desc(AUC.test))
folder <- paste0("model_dist_norm_", hours, "h")
folder_name <- folder

# ** Save performance results ----------------------------------------------
#results details
res_summ <- table(sub("\\_.*", "", unique(leaderboard$fname))) %>%
  as.data.frame() %>%
  setNames(c("algorithm", "n"))

write.csv(res_summ, here("output", "Table_S_results_summary_dist.csv"), row.names = FALSE)

#number of models to choose from
n=nrow(leaderboard)
auc.df <- data.frame(n = 1:n, 
                     algorithm =NA,
                     fname=NA,
                     cutoff = NA,
                     AUC.train=NA,
                     AUC.test=NA)

for (i in 1:n){
  auc.df$fname[i] <- leaderboard$fname[i]
  auc.df$cutoff[i] <- leaderboard$cutoff[i]
  mod <- h2o.loadModel(here("output", folder_name, leaderboard$fname[i]))
  
  auc.df$algorithm[i] <- mod@algorithm
  
  res.train <- h2o.predict(mod , train)
  res.train <- as.data.frame(res.train)$p1
  
  res.test <- h2o.predict(mod , test)
  res.test <- as.data.frame(res.test)$p1
  
  #all values lower than cutoff value will be classified as 0 (NT in this case)
  train.labels <- ifelse(res.train < leaderboard$cutoff[i], 0, 1)
  test.labels <- ifelse(res.test < leaderboard$cutoff[i], 0, 1)
  true.train.labels <- as.data.frame(train)$iucn
  true.test.labels <- as.data.frame(test)$iucn
  
  auc.df$AUC.train[i] <- as.numeric(auc(true.train.labels, train.labels))
  auc.df$AUC.test[i] <- as.numeric(auc(true.test.labels, test.labels))
}

auc.df <- auc.df %>% #filter(cutoff > 0.32) %>%
  group_by(fname) %>% 
  filter(AUC.test == max(AUC.test)) %>%
  arrange(desc(AUC.test)) 

auc.df <- auc.df[duplicated(auc.df$fname)==FALSE,] #remove duplicates

#checking for precision and recall
train.precision <- c()
test.precision <- c()

train.recall <- c()
test.recall <- c()

for (w in 1:10){
  win <- auc.df$n[w]
  aml_leader <- h2o.loadModel(here("output", folder_name, leaderboard$fname[win]))
  
  res.train <- h2o.predict(aml_leader , train)
  res.train <- as.data.frame(res.train)$p1
  
  res.test <- h2o.predict(aml_leader , test)
  res.test <- as.data.frame(res.test)$p1
  
  #all values lower than cutoff value will be classified as 0 (NT in this case)
  train.labels <- ifelse(res.train < leaderboard$cutoff[win], 0, 1)
  test.labels <- ifelse(res.test < leaderboard$cutoff[win], 0, 1)
  true.train.labels <- as.data.frame(train)$iucn
  true.test.labels <- as.data.frame(test)$iucn
  
  train.precision[w] <- precision(as.factor(train.labels), true.train.labels, levels=c(1,0))
  test.precision[w] <- precision(as.factor(test.labels), true.test.labels, levels=c(1,0))
  
  train.recall[w] <- recall(as.factor(train.labels), true.train.labels, levels=c(1,0))
  test.recall[w] <- recall(as.factor(test.labels), true.test.labels, levels=c(1,0))
  
  
}

train.Fscore = 2 * train.precision * train.recall / (train.precision + train.recall)
test.Fscore = 2 * test.precision * test.recall / (test.precision + test.recall)

auc.df <- cbind(auc.df[1:10,], train.Fscore=train.Fscore, test.Fscore=test.Fscore)
write.csv(auc.df, here("output", "Table_S_model_performance_dist.csv"), row.names=FALSE)

# ** Load performance results ---------------------------------------------
# reading output
auc.df <- read.csv(here("output", "Table_S_model_performance_dist.csv")) %>% filter(cutoff > 0.32)
w=2 #winning model

win <- auc.df$n[w]
aml_leader <- h2o.loadModel(here("output", folder_name, leaderboard$fname[win]))

#### TEXT FOR PAPER
paste0("Out of the ", length(unique(leaderboard$fname)), " models trained by the AutoML, the optimal model was a ", aml_leader@algorithm, 
       " with AUC of ", round(auc.df$AUC.train[w],2),  " on the training dataset and ", round(auc.df$AUC.test[w],2), " on the test dataset and", 
       " an average classification probability threshold (decision threshold) of ", auc.df$cutoff[w])

# * Confusion Matrix ------------------------------------------------------
res.train <- h2o.predict(aml_leader , train)
res.train <- as.data.frame(res.train)$p1

res.test <- h2o.predict(aml_leader , test)
res.test <- as.data.frame(res.test)$p1

#all values lower than cutoff value will be classified as 0 (NT in this case)
train.labels <- ifelse(res.train < leaderboard$cutoff[win], 0, 1)
test.labels <- ifelse(res.test < leaderboard$cutoff[win], 0, 1)
true.train.labels <- as.data.frame(train)$iucn
true.test.labels <- as.data.frame(test)$iucn

cf.train <- caret::confusionMatrix(as.factor(train.labels), true.train.labels)
cf.test <- caret::confusionMatrix(as.factor(test.labels), true.test.labels)

acc.train <- prop.table(cf.train$table, 2)
acc.test <- prop.table(cf.test$table, 2)

### TEXT FOR PAPER
paste0("The model correctly classified ", round(acc.train[2,2], 2)*100,"% and ",round(acc.test[2,2], 2)*100,
       "% of the threatened species for the training and test dataset respectively.")
paste0("It had a higher classification accuracy for non-threatened species, ",round(acc.train[1,1], 2)*100,
       "% and ",round(acc.test[1,1], 2)*100,"%, and overall accuracy of ", round(cf.train$overall[1], 2)*100,"% and ", round(cf.test$overall[1], 2)*100,
       "%, respectively on the training and test datasets.")

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
  mutate(ints = cut(corallite_ori ,breaks = c(0, 5, max(corallite_ori)))) %>%
  group_by(iucn, pred, ints) %>%
  tally() %>%
  ungroup() %>%
  group_by(ints) %>%
  mutate(prop = n/sum(n)) %>%
  filter(iucn != pred) %>%
  group_by(ints) %>% summarise(sum(prop))

# * Model agnostic metrics

# create a data frame with just the features
features <- as.data.frame(df) %>% dplyr::select(x)

#  Create a numeric vector with the actual responses
response <- ifelse(df$iucn == "NT", 0, 1)

# Create custom predict function that returns the predicted values as a vector
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results[[3L]])
}


# * Global feature importance ---------------------------------------------
# ** Fig 2: Variable Importance -------------------------------------------
if(aml_leader@algorithm != "stackedensemble") {
  aml_varimp <- h2o.varimp(aml_leader) 
} else {
  #get base model contributing to the most to the stackedensemble
  meta <- h2o.getModel(aml_leader@model$metalearner$name)
  meta <- sort(meta@model$coefficients, decreasing = TRUE)[1]
  meta <- h2o.getModel(names(meta))
  aml_varimp <- h2o.varimp(meta) 
}

aml_varimp <- aml_varimp %>%
  as.data.frame() %>%
  filter(variable != "branching: missing(NA)") %>%
  mutate(rank = nrow(.):1) %>%
  mutate(variable = gsub("\\.", ": ", variable),
         variable = gsub("_", " ", variable))

p_vi <- aml_varimp %>% 
  as.data.frame() %>%
  filter(variable != "branching: missing(NA)") %>%
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
            axis.text = element_text(size = 10),
            panel.grid = element_blank(), 
            legend.position = "none") 
  }
  
  return(p)
}

# plot pdp
p1 <- pdp_plot(pdp_corallite, olddata=df$corallite_ori, col = u_col[2]) + labs(x="Corallite Diameter (mm)") + 
  coord_cartesian(xlim=c(0, 25), ylim=c(0.2,0.6)) 

p2 <- pdp_plot(pdp_branching, unscale = FALSE, levels=c("NB", "LB", "MB", "HB")) + scale_fill_manual(values=u_col[2:5]) + 
  labs(x="Degree of branching")  

svg(here("figs", "Fig_S_model_metrics_morph.svg"), w=10, h=5)
p1 + p2  + plot_annotation(tag_levels = "A") + plot_layout(ncol=2)
dev.off()

# DD species ------------------------------------------------
dd <- df.corals %>% filter(iucn =="Data Deficient") %>% dplyr::select (-iucn)

# * Prediction of threat status -------------------------------------------

res <- h2o.predict(aml_leader, as.h2o(dd))
res <- as.data.frame(res)$p1

# calculate threat status based on threshold cutoff
dd$status <- ifelse(res < leaderboard$cutoff[win], "NT", "T")

#number of missing values for each species
dd$na_count <- apply(dd[,x], 1, function(x) sum(is.na(x)))
dd[dd$na_count >1,]$status <- "DD"

table(df$iucn) %>% prop.table()
table(dd$status)[-1] %>% prop.table()
table(c(df$iucn, dd$status))[-1] %>% prop.table()

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

p <- ggplot(dd_plot, aes(fill=status, y=n3, x=reorder(genus, n2))) +
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
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        panel.grid = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(2, 1, 1, 1), "lines"))

svg(here("figs", "Fig_S_morph_dd_species.svg"), w=6, h=6)
p 
grid::grid.text("Threatened", x = unit(0.45, "npc"), y = unit(.95, "npc"), 
                gp=gpar(family="Roboto Mono", fontface="bold"), hjust=-1.04)
grid::grid.text("Not Threatened", x = unit(0.45, "npc"), y = unit(.95, "npc"), 
                gp=gpar(family="Roboto Mono", fontface="bold"), hjust=0.3)
dev.off()

# * Map of DD distribution risk -------------------------------------------
gr <- hexagrid(16, sp=TRUE) #2.5 degrees
world <- rworldmap::getMap()
cuts <- seq(0.0, 1, 0.2)
map_col <- colorRampPalette(u_col[c(2,4,5)])(length(cuts)-1)


#### select regions
sp.pol <- gr@sp
sp.df <- data.frame(faces = rownames(gr@faces))

rownames(sp.df) <- sp.df$faces
sp.pol <- SpatialPolygonsDataFrame(sp.pol, sp.df)

#crop according to regions
carribean_ext <- extent(-100.546875,-40.649414,6.271618,37.230328) 
indo_pacific_ext <- extent(28.125000,240.820313,-45.336702,22.593726)
carribean <- crop(sp.pol, carribean_ext)
indo_pacific <- crop(sp.pol, indo_pacific_ext)


#### DD only
dd_occ <- read.csv(here("data", "2019-11-05_obis_scleractinia.csv"), stringsAsFactors = FALSE) %>%
  filter(!between(decimalLatitude, 39, 45) & !between(decimalLongitude, -38, 30)) %>% #remove mediterranean ones
  filter(scientificName %in% dd$sp) %>% #select only occurrences for dd species
  distinct(scientificName, decimalLatitude, decimalLongitude) %>%
  left_join(dd %>% dplyr::select(sp, status), by=c("scientificName" = "sp"))

dd_occ$cell <- locate(gr, dd_occ[, c(3,2)])
dd_occ$region <- ifelse(dd_occ$cell %in% carribean$faces, "carribean", "indo-pacific")

dd_occ <- dd_occ %>% distinct(region, scientificName, status) %>%
  filter(status != "DD") 

##### DS only
ds_occ <- read.csv(here("data", "2019-11-05_obis_scleractinia.csv"), stringsAsFactors = FALSE) %>%
  filter(!between(decimalLatitude, 39, 45) & !between(decimalLongitude, -38, 30)) %>% #remove mediterranean ones
  filter(scientificName %in% df$sp) %>% #select data sufficient species only
  distinct(scientificName, decimalLatitude, decimalLongitude) %>%
  left_join(df %>% dplyr::select(sp, iucn) %>% mutate(status = iucn) %>% dplyr::select(-iucn), by=c("scientificName" = "sp"))

ds_occ$cell <- locate(gr, ds_occ[, c(3,2)])

ds_occ$region <- ifelse(ds_occ$cell %in% carribean$faces, "carribean", "indo-pacific")

ds_occ <- ds_occ %>% distinct(region, scientificName, status) %>%
  filter(status != "DD") 

#### DS + DD
all_occ <- rbind(ds_occ, dd_occ)

#results
table(dd_occ$region, dd_occ$status) %>% prop.table(1)
table(ds_occ$region, ds_occ$status) %>% prop.table(1)
table(all_occ$region, all_occ$status) %>% prop.table(1)

table(all_occ$status) %>% prop.table()
table(ds_occ$status) %>% prop.table()


# Plio-Pleistocene --------------------------------------------------------

pleist.df <- read.csv(here("data", "pleist_resolved.csv"), stringsAsFactors = FALSE) %>%
  mutate(corallite = as.numeric(gsub(">", "", corallite)))

pleist.df <- pleist.df %>% group_by(valid_name) %>%
  summarise(corallite = max(corallite), max_depth=max(max_depth), range = max(prop_range), #prop_range = max(prop_range), 
            branching = unique(branching)[1], globally.extinct = unique(globally.extinct), 
            regionally.extinct = unique(regionally.extinct)) 

#scale based on original data
pleist.df <- pleist.df %>% 
  mutate(corallite = as.numeric(scaled.new(corallite, df.corals$corallite_ori)),
         max_depth = as.numeric(scaled.new(max_depth, df.corals$depth_ori))
         #range = as.numeric(scaled.new(max_depth, scale(df.corals$range_ori)))
  )

pleist.df[pleist.df$branching == "",]$branching <- NA 

pleist.df[pleist.df$globally.extinct == "extinct" & pleist.df$regionally.extinct == "extant",]$regionally.extinct <- "extinct"
res <- h2o.predict(aml_leader, as.h2o(pleist.df))
res <- as.data.frame(res)$p1

# calculate threat status based on threshold cutoff
pleist.df$status <- ifelse(res < leaderboard$cutoff[win], "NT", "T")

#number of missing values for each species
pleist.df$na_count <- apply(pleist.df[,x], 1, function(x) sum(is.na(x)))
pleist.df[pleist.df$na_count >1,]$status <- "DD"

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

#ggsave(here("figs", "Fig_06_pleist_predictions.svg"), w=6, h=6)

### per genus
df.ds <- df.corals %>%
  na.omit() %>%
  dplyr::select(sp, iucn) %>%
  #remove all Data Deficient corals
  filter(iucn != "Data Deficient") %>%
  
  #treat as binary
  mutate(
    #reclassify iucn status 
    status = case_when(iucn == "Least Concern" ~ "NT",         
                       iucn == "Vulnerable" ~ "T",           
                       iucn == "Near Threatened"   ~ "NT",    
                       iucn == "Endangered"   ~ "T",         
                       iucn == "Critically Endangered"~ "T")) 

df.dd <- dd %>% 
  filter(status != "DD") %>%
  dplyr::select(sp, status) %>%
  mutate(iucn = "DD")

df <- rbind(df.ds, df.dd)

df.comp <- pleist.df %>% 
  filter(status != "DD") %>%
  separate(valid_name, c("genus", "species")) %>%
  group_by(genus, status) %>%
  tally() %>%
  reshape2::dcast(genus~status) %>%
  group_by(genus) %>%
  replace_na(list(T=0)) %>%
  mutate(tot = sum(T, NT, na.rm=TRUE),
         T=T/tot) %>% 
  dplyr::select(genus, T, tot) %>%
  rename(fossil = T) %>%
  left_join(df %>% separate(sp, c("genus", "species")) %>% 
              group_by(genus, status) %>%
              tally() %>%
              reshape2::dcast(genus~status) %>%
              group_by(genus) %>%
              replace_na(list(T=0)) %>%
              mutate(totm = sum(T, NT, na.rm=TRUE),
                     T=T/totm)%>% 
              dplyr::select(genus, T, totm) %>%
              rename(modern = T), by="genus") %>%
  na.omit()

ggplot(df.comp, aes(x=fossil, y=modern, label=paste0(genus, "(", tot, "/", totm, ")"))) +
  geom_point() +
  geom_text_repel(cex=3.2, vjust=-0.1, hjust=-0.1) +
  labs(x="Proportion of Plio-Plestocene species threatened",
       y="Proportion of modern species threatened") +
  # geom_hline(yintercept = 0.305, linetype="dashed") +
  # geom_vline(xintercept = 0.305, linetype="dashed") +
  geom_abline(slope=1, intercept=0) +
  theme_light(base_size = 15) +
  theme(axis.title = element_text(size = 12, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(family = "Roboto Mono", size = 10),
        axis.text = element_text(family = "Roboto Mono", size = 10),
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal")

#ggsave(here("figs", "Fig_S_plio_modern_comp.svg"), w=7, h=7)
