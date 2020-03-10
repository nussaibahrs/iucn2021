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

df <- subset(na.omit(df.corals[,-1]), iucn != "Data Deficient")
df$branching <-  factor(df$branching, levels=c("NB", "LB", "MB", "HB"))
df$iucn <- as.factor(df$iucn)

levels(df$iucn) <- c("EN", "EN", "LC", "NT", "VU")
levels(df$iucn) <- c("LC", "NT", "VU", "EN")

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
x <- c("max_depth", "branching", "corallite", "range" )

# Load Model --------------------------------------------------------------

# * Load best model in h2o session ------------------------------------

folder_name <- "model_multi"
leaderboard <- read.csv(here("output", folder_name, "leaderboard.csv"), stringsAsFactors = FALSE)
aml_leader <- h2o.loadModel(here("output", folder_name, leaderboard$model_id[1])) #load top one
aml_leader@algorithm

# * Confusion Matrix ------------------------------------------------------
res.train <- h2o.predict(aml_leader, train)
res.test <- h2o.predict(aml_leader , test)

#all values lower than cutoff value will be classified as 0 (NT in this case)
train.labels <- as.data.frame(res.train)$predict
test.labels <- as.data.frame(res.test)$predict
true.train.labels <- as.data.frame(train)$iucn
true.test.labels <- as.data.frame(test)$iucn

cf.train <- caret::confusionMatrix(as.factor(train.labels), true.train.labels)
cf.test <- caret::confusionMatrix(as.factor(test.labels), true.test.labels)

prop.table(cf.train$table, 2)
prop.table(cf.test$table, 2)


# * Misclassification -----------------------------------------------------
train.df <- train %>% as.data.frame %>% 
  mutate(pred = train.labels)

test.df <- test %>% as.data.frame %>% 
  mutate(pred= test.labels)


df <- bind_rows(train.df, test.df)

# per iucn status
df %>% group_by(iucn, pred) %>%
  tally() %>%
  group_by(iucn) %>%
  mutate(tot=sum(n), 
         prop = n/sum(n)) %>%
  ungroup %>%
  reshape2::dcast(iucn~pred, value.var = "prop")

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

p_vi

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
  coord_cartesian(xlim=c(0, 25)) 

p2 <- pdp_plot(pdp_branching, unscale = FALSE, levels=c("NB", "LB", "MB", "HB")) + scale_fill_manual(values=u_col[2:5]) + 
  labs(x="Degree of branching")  

p3 <- pdp_plot(pdp_range, olddata=scale(df$range_ori), col = u_col[3]) + labs(x="Range (normalised)") + 
  geom_hline(yintercept = opt_par_final$cutoff[win], linetype="dashed", col="darkgrey")

p4 <- pdp_plot(pdp_max_depth, olddata=scale(df$depth_ori), col = u_col[4]) + labs(x="Maximum depth (m)") +
  coord_cartesian(xlim=c(0, 100)) 


p1 + p2 + p3 + p4  + plot_annotation(tag_levels = "A")

# DD species ------------------------------------------------
dd <- df.corals %>% filter(iucn =="Data Deficient") %>% select (-iucn) 

# * Prediction of threat status -------------------------------------------

res <- h2o.predict(aml_leader, as.h2o(dd))

# calculate threat status based on threshold cutoff
dd$status <- as.data.frame(res)$predict %>% as.character()

#number of missing values for each species
dd$na_count <- apply(dd[,x], 1, function(x) sum(is.na(x)))
dd[dd$na_count >2,]$status <- "DD"

#get genus for each species
dd$genus <- gsub( " .*$", "", dd$sp)

# ** Fig. 4: Threat status by genera (no. of sp > 2) -----------------------
dd %>% filter(status != "DD") %>%
  group_by(genus, status) %>%
  tally() %>%
  group_by(genus) %>%
  mutate(n2 = sum(n),
         n3 = n/n2) %>%
  arrange(desc(n2)) %>%
  filter(n2 > 1) 


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
res <- as.data.frame(res)$predict

# calculate threat status based on threshold cutoff
pleist.df$status <- res %>% as.character()

#number of missing values for each species
pleist.df$na_count <- apply(pleist.df[,x], 1, function(x) sum(is.na(x)))
pleist.df[pleist.df$na_count >2,]$status <- "DD"

prop.table(table(pleist.df$status[pleist.df$status != "DD"]))

pleist.summ <- pleist.df %>% 
  filter(status != "DD") %>%
  mutate(globally.extinct = tools::toTitleCase(globally.extinct),
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


### per genus
df.ds <- df.corals %>%
  na.omit() %>%
  select(sp, iucn) %>%
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
  select(sp, status) %>%
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
  select(genus, T, tot) %>%
  rename(fossil = T) %>%
  left_join(df %>% separate(sp, c("genus", "species")) %>% 
  group_by(genus, status) %>%
  tally() %>%
  reshape2::dcast(genus~status) %>%
    group_by(genus) %>%
    replace_na(list(T=0)) %>%
    mutate(totm = sum(T, NT, na.rm=TRUE),
           T=T/totm)%>% 
    select(genus, T, totm) %>%
    rename(modern = T), by="genus") %>%
  na.omit()
  
library(ggrepel)

ggplot(df.comp, aes(x=fossil, y=modern, label=paste0(genus, "(", tot, "/", totm, ")"))) +
  geom_point() +
  geom_text_repel(cex=3.2, vjust=-0.1, hjust=-0.1) +
  labs(x="Proportion of Plio-Plestocene species threatened",
       y="Proportion of modern species threatened") +
  # geom_hline(yintercept = 0.305, linetype="dashed") +
  # geom_vline(xintercept = 0.305, linetype="dashed") +
  geom_abline(slope=1, intercept=0) +
  geom_abline(slope=c(2, 0.5), intercept = 0, linetype="dashed") +
  theme_light(base_size = 15) +
  theme(axis.title = element_text(size = 12, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(family = "Roboto Mono", size = 10),
        axis.text = element_text(family = "Roboto Mono", size = 10),
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal")

