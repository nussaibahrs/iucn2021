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

source(file.path("code", "functions.R"))
#palette
u_col <- ggsci::pal_uchicago("default")(9)[c(2,5,4,3,1)]

# Data  and set up h20 session -------------------------------------------------------------------
source(file.path("code", "01-model", "00-load_data.R"))

# ** Load performance results ---------------------------------------------
# reading output
auc.df <- read.csv(file.path("output", "Table_S_model_performance_fmm.csv")) 


w=1 #winning model
win <- auc.df$n[w]

leaderboard <- read.csv(file.path("output", folder, "leaderboard.csv"), stringsAsFactors = FALSE)
aml_leader <- h2o.loadModel(file.path("output", folder, auc.df$fname[w]))

#plot for AUC cutoff
cutoff <- seq(0,1,0.01)


mod <- data.frame(cutoff=cutoff,
                  train = NA,
                  test = NA)


for (c in 1:length(cutoff)){
  res.train <- h2o.predict(aml_leader , train)
  res.train <- as.data.frame(res.train)$p1

  res.test <- h2o.predict(aml_leader , test)
  res.test <- as.data.frame(res.test)$p1

  #all values lower than cutoff value will be classified as 0 (NT in this case)
  train.labels <- ifelse(res.train < cutoff[c], 0, 1)
  test.labels <- ifelse(res.test < cutoff[c], 0, 1)
  true.train.labels <- as.data.frame(train)$extinct
  true.test.labels <- as.data.frame(test)$extinct

  mod$train[c] <- as.numeric(auc(true.train.labels, train.labels))
  mod$test[c] <- as.numeric(auc(true.test.labels, test.labels))
}

max.train <- mod$cutoff[which.max(mod$train)]
max.test <- mod$cutoff[which.max(mod$test)]

mod_cutoff <- cbind.data.frame(value=c(max.train, max.test), variable=c("train", "test"))

p1 <- reshape2::melt(mod, id.vars="cutoff") %>%
  ggplot(aes(x=cutoff, y=value, col=variable)) +
  geom_rect(aes(xmin=0, xmax=1, ymin=0, ymax=0.7),fill="lightgrey", alpha=0.02, col=NA)+
  geom_line() +
  geom_vline(data=mod_cutoff, aes(xintercept=value, col=variable, 
                                  linetype=variable), show.legend = FALSE) +
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

ggsave(file.path("figs", "Fig_S4_final_model.svg"), p1, w=5, h=4)


# * Confusion Matrix ------------------------------------------------------
res.train <- h2o.predict(aml_leader , train)
res.train <- as.data.frame(res.train)$p1

res.test <- h2o.predict(aml_leader , test)
res.test <- as.data.frame(res.test)$p1

#all values lower than cutoff value will be classified as 0 (NT in this case)
train.labels <- ifelse(res.train < auc.df$cutoff[w], 0, 1)
test.labels <- ifelse(res.test < auc.df$cutoff[w], 0, 1)
true.train.labels <- as.data.frame(train)$extinct
true.test.labels <- as.data.frame(test)$extinct

cf.train <- caret::confusionMatrix(as.factor(train.labels), true.train.labels, positive="1")
cf.test <- caret::confusionMatrix(as.factor(test.labels), true.test.labels, positive = "1")

acc.train <- prop.table(cf.train$table, 2)
acc.test <- prop.table(cf.test$table, 2)

write.csv(cbind(true = c(as.character(true.train.labels), 
                         as.character(true.test.labels)), 
                pred=c(train.labels, 
                       test.labels)), 
          file.path("output", "fossil_results.csv"), row.names = FALSE)

# Model agnostic metrics ---------------------------------------------

# create a data frame with just the features
features <- as.data.frame(df) %>% dplyr::select(x)

#  Create a numeric vector with the actual responses
response <- df$extinct


# * Global feature importance ---------------------------------------------

# ** generate model agnostic object to explain response of features -------

# Create custom predict function that returns the predicted values as a vector
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(as.numeric(as.character(results[[3L]])))
}

dalex_explainer <- DALEX::explain(
  model = aml_leader,
  data = features,
  y = as.numeric(as.character(response)),
  predict_function = pred
)

# ** Fig 2: Variable Importance -------------------------------------------
set.seed(42)
aml_varimp  <- model_parts(explainer = dalex_explainer, 
                           loss_function = loss_root_mean_square,
                           B = 50, #number of permutations
                           type = "difference")

varimp.df <- cbind.data.frame(loss=aml_varimp$dropout_loss, p=aml_varimp$permutation, lab=aml_varimp$variable) %>% 
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

ggsave (file.path("figs", "Fig_02_var_imp.svg"), p_vi, h=4, w=8)

# *  Partial dependence profile ---------------------------------------------------
# ** generate the partial dependence profiles for each variable ------------
pdp_corallite <- ingredients::partial_dependency(dalex_explainer, variables = "corallite")
pdp_branching <- ingredients::partial_dependency(dalex_explainer, variables = "branching", variable_type="categorical")
pdp_integration <- ingredients::partial_dependency(dalex_explainer, variables = "integration")
pdp_budding <- ingredients::partial_dependency(dalex_explainer, variables = "budding", variable_type="categorical")
pdp_family <- ingredients::partial_dependency(dalex_explainer, variables = "family", variable_type="categorical")

# ** Fig. 3: Partial dependency Plots -------------------------------------

p1 <- pdp_plot(pdp_corallite, col = u_col[2], unscale = FALSE) +
  labs(x="Corallite Diameter (mm)") +
  coord_cartesian(xlim=c(0, 25))

p2 <- pdp_plot(pdp_branching, unscale = FALSE, levels=c("NB", "LB", "MB", "HB")) + scale_fill_manual(values=u_col[2:5]) +
  labs(x="Degree of branching")

p3 <- pdp_plot(pdp_integration, unscale = FALSE)+
  labs(x="Level of corallite integration")

p4 <- pdp_plot(pdp_budding, unscale = FALSE, levels=c("intra", "extra", "both"))+
  labs(x="Budding type") + scale_fill_manual(values=u_col[2:5]) + 
  scale_x_discrete(labels=c("intra", "extra", "both"))

svg(file.path("figs", "Fig_02_model_metrics.svg"), w=8, h=8)
p3 + p1 + p2 + p4 + plot_annotation(tag_levels = "A")
dev.off()

