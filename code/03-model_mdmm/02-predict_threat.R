# Set up environment ------------------------------------------------------
# helper packages
library(tidyverse)
library(magrittr)

# setting up machine learning models
library(h2o) # version 3.28.0.2

# packages for explaining our ML models
library(pROC)
library(hmeasure)
library(pdp)
library(DALEX)
library(ingredients)
library(caret)

# plotting packages
library(icosa)
library(ggthemes)
library(patchwork)
library(grid)
library(raster)
library(ggpubr)

# Data  -------------------------------------------------------------------
#palette
u_col <- ggsci::pal_uchicago("default")(9)[c(2,5,4,3,1)]


# Load data ---------------------------------------------------------------
source(file.path("code", "03-modern_corals_mdmm", "00-load_data.R"))

folder_name <- folder

leaderboard <- read.csv(file.path("output", paste0("Table_S_maximised_threshold_mdmm_", hours, "h.csv")), stringsAsFactors = FALSE) %>% arrange(desc(AUC.test))


# ** Save performance results ----------------------------------------------
#results details
res_summ <- table(sub("\\_.*", "", unique(leaderboard$fname))) %>%
  as.data.frame() %>%
  setNames(c("algorithm", "n"))

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
  mod <- h2o.loadModel(file.path("output", folder_name, leaderboard$fname[i]))

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

auc.df <- auc.df %>% filter(cutoff > 0.32) %>%
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
aml_leader <- h2o.loadModel(file.path("output", folder_name, leaderboard$fname[win]))

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
write.csv(auc.df, file.path("output", "Table_S_model_performance_mdmm.csv"), row.names=FALSE)


# ** Load performance results ---------------------------------------------
# reading output
auc.df <- read.csv(file.path("output", "Table_S_model_performance_binary.csv"))
w=1 #winning model
win <- auc.df$n[w]
aml_leader <- h2o.loadModel(file.path("output", folder_name, leaderboard$fname[win]))

#plot for AUC cutoff
cutoff <- seq(0,1,0.01)


mod <- data.frame(cutoff=cutoff,
                  train = NA,
                  test = NA)

mod_cutoff <- leaderboard[leaderboard$fname == auc.df$fname[w] & leaderboard$cutoff == auc.df$cutoff[w], c("cutoff.train", "cutoff.test")] %>%
  setNames(c("train", "test")) %>%
  reshape2::melt()


for (c in 1:length(cutoff)){
  res.train <- h2o.predict(aml_leader , train)
  res.train <- as.data.frame(res.train)$p1

  res.test <- h2o.predict(aml_leader , test)
  res.test <- as.data.frame(res.test)$p1

  #all values lower than cutoff value will be classified as 0 (NT in this case)
  train.labels <- ifelse(res.train < cutoff[c], 0, 1)
  test.labels <- ifelse(res.test < cutoff[c], 0, 1)
  true.train.labels <- as.data.frame(train)$iucn
  true.test.labels <- as.data.frame(test)$iucn

  mod$train[c] <- as.numeric(auc(true.train.labels, train.labels))
  mod$test[c] <- as.numeric(auc(true.test.labels, test.labels))
}

p1 <- reshape2::melt(mod, id.vars="cutoff") %>%
  ggplot(aes(x=cutoff, y=value, col=variable)) +
  geom_rect(aes(xmin=0, xmax=1, ymin=0, ymax=0.7),fill="lightgrey", alpha=0.02, col=NA)+
  geom_line() +
  geom_vline(data=mod_cutoff, aes(xintercept=value, col=variable, linetype=variable), show.legend = FALSE) +
  coord_cartesian(ylim=c(0,1), expand=FALSE)+
  scale_color_manual(values=u_col[c(2,5)], labels=c("Training", "Test")) +
  labs(x="Cutoff Threshold", y="AUC", col="Dataset", linetype="Dataset")+
  theme_light(base_size = 15) +
  theme(axis.title = element_text(size = 12, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        panel.grid = element_blank(),
        # panel.border = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal")

ggsave(file.path("figs", "Fig_S_final_model.svg"), p1, w=5, h=4)


# * Model agnostic metrics ----

# create a data frame with just the features
features <- as.data.frame(df) %>% dplyr::select(all_of(x))

#  Create a numeric vector with the actual responses
response <- ifelse(df$iucn == "NT", 0, 1)

# Create custom predict function that returns the predicted values as a vector
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results[[3L]])
}

dalex_explainer <- DALEX::explain(
  model = aml_leader,
  data = features,
  y = as.numeric(as.character(response)),
  predict_function = pred
)


# * Global feature importance ---------------------------------------------
# ** Fig 2: Variable Importance -------------------------------------------
set.seed(42)
aml_varimp  <- model_parts(explainer = dalex_explainer, 
                           loss_function = loss_root_mean_square,
                           B = 100, #number of permutations
                           type = "difference")

varimp.df <- cbind.data.frame(loss=aml_varimp$dropout_loss, 
                              p=aml_varimp$permutation, 
                               lab=aml_varimp$variable) %>% 
  filter(lab %in% x)

varimp.df$lab <- plyr::mapvalues(varimp.df$lab, x, labsx)

varimp.summary <- varimp.df %>% group_by(lab) %>%  
  summarise(res=mean(loss)) %>% 
  arrange(desc(res))

p_vi <- ggplot() + 
  geom_bar(data=varimp.summary, aes(x=reorder(lab, res), y=res), stat="identity", fill=u_col[1], width=.6, alpha=0.8) +
  geom_boxplot(data=varimp.df, aes(x=lab, y=loss), fill=u_col[5], col=u_col[5], width=0.4, outlier.colour = NA) +
  scale_y_continuous(expand=expand_scale(mult = c(0, .1))) +
  coord_flip() +
  labs(x="Variables", y="Root mean square error (RMSE) after 1000 permutations") +
  theme_light(base_size = 15) +
  theme(axis.title = element_text(size = 12, face="bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        panel.grid = element_blank()
  )

p_vi

ggsave (file.path("figs", "Fig_S2_var_imp_mdmm.svg"), p_vi, h=4, w=8)

# DD species ------------------------------------------------
dd <- df.corals %>% filter(iucn =="Data Deficient") %>% dplyr::select (-iucn)

# * Prediction of threat status -------------------------------------------

res <- h2o.predict(aml_leader, as.h2o(dd))
res <- as.data.frame(res)$p1

# calculate threat status based on threshold cutoff
dd$status <- ifelse(res < leaderboard$cutoff[win], "NT", "T")

#number of missing values for each species
dd$na_count <- apply(dd[,x], 1, function(x) sum(is.na(x)))
dd[dd$na_count >2,]$status <- "DD"

table(dd$status)

table(df$iucn) %>% prop.table()
table(dd$status) %>% prop.table()
df$iucn2 <- df$iucn
table(c(ifelse(df$iucn==0, "NT", "T"), 
        dd$status)) %>% prop.table()

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
        axis.text.y = element_text(face="italic"),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(2, 1, 1, 1), "lines"))

# svg(file.path("figs", "Fig_04_dd_species.svg"), w=6, h=6)
p
grid::grid.text("Threatened", x = unit(0.45, "npc"), y = unit(.95, "npc"),
                gp=gpar(family="Roboto Mono", fontface="bold"), hjust=-1.04)
grid::grid.text("Not Threatened", x = unit(0.45, "npc"), y = unit(.95, "npc"),
                gp=gpar(family="Roboto Mono", fontface="bold"), hjust=0.3)


# Plio-Pleistocene --------------------------------------------------------
pleist.tr <- read.csv(file.path("data", "pbdb_pleist.csv"))

pleist.df <- read.csv(file.path("data", "pleist_resolved.csv"), stringsAsFactors = FALSE) %>%
  dplyr::select(valid.name=valid_name, max_depth, range)

pleist.df <- merge(pleist.tr, pleist.df)

res <- h2o.predict(aml_leader, as.h2o(pleist.df))
res <- as.data.frame(res)$p1

# calculate threat status based on threshold cutoff
pleist.df$status <- ifelse(res < leaderboard$cutoff[win], "extant", "extinct")

length(which(pleist.df$extinct != pleist.df$status))/nrow(pleist.df) #mismatch

pleist.t <- pleist.df[pleist.df$extinct == "extinct",]
length(which(pleist.t$extinct != pleist.t$status))/nrow(pleist.t)
