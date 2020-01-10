# Initialise environment --------------------------------------------------
# * load libraries ----
library(here)
library(tidyverse)
library(ggridges)
library(patchwork)

# * palette ----
u_col <- ggsci::pal_uchicago("default")(9)[c(2,5,4,3,1)]
u_cat <- u_col[c(1,5,2)]
#pie(1:5, col=u_col)

# * load and set up data ----
df.corals <- read.csv(here("data", "traits_iucn.csv"), stringsAsFactors = FALSE) %>%
  na.omit() %>%
  mutate(
    #reclassify iucn status 
    iucn = case_when(iucn == "Least Concern" ~ "NT",         
                     iucn == "Vulnerable" ~ "NT",           
                     iucn == "Near Threatened"   ~ "NT",    
                     iucn == "Endangered"   ~ "T",         
                     iucn == "Critically Endangered"~ "T",
                     iucn == "Data Deficient"~ "DD"),
    iucn = factor(iucn, levels = c("DD", "T", "NT")),
    branching = factor(branching, levels=c("HB", "MB", "LB", "NB"))
  )


# Figure 01 - Individual traits and their link to extinction risk  --------
# set theme
theme_set(theme_light(base_size = 15))

# * range size ----

# calculate world average
world_avg <- df.corals %>%
  summarize(avg = mean(range, na.rm = T)) %>%
  pull(avg)

reg_avg <- df.corals %>% group_by(iucn) %>%
  summarise(n=mean(range, na.rm=TRUE))

p1 <- ggplot(df.corals, aes(x = iucn, y = range, color = iucn)) +
  coord_flip() +
  #scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0.005)) +
  scale_color_manual(values=u_cat, guide=FALSE) + 
  labs(x = NULL, y = "Range (km)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face="bold"),
        axis.text.x = element_text(family = "Roboto Mono", size = 10),
        panel.grid = element_blank()) + 
  geom_hline(aes(yintercept = world_avg), color = "gray70", size = 0.6) +
  geom_segment(data=reg_avg, aes(x = iucn, xend = iucn,
                                 y = world_avg, yend = n, col=iucn),
               size = 0.8, inherit.aes = FALSE) +
  stat_summary(fun.y = mean, geom = "point", size = 5) +
  geom_jitter(size = 2, alpha = 0.1, width = 0.2) 

# * corallite diameter ----
p2 <- ggplot(na.omit(df.corals), aes(x = iucn, y = corallite, color = iucn)) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 11), expand = c(0.005, 0.005)) +
  scale_color_manual(values=u_cat, guide=FALSE) + 
  labs(x = NULL, y = "Corallite Diameter (mm)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face="bold"),
        axis.text.x = element_text(family = "Roboto Mono", size = 10),
        panel.grid = element_blank()) +
  #should remove outliers
  geom_boxplot(color="gray20", outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.1, width = 0.2) 

# * branching ----
p3 <- df.corals  %>% 
  group_by(iucn, branching) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count)) %>%
  ggplot(aes(iucn, fill=branching, y=perc*100)) +
  geom_bar(stat="identity", width = 0.5) +
  coord_flip() +
  scale_y_continuous(limits=c(0, 100), expand=c(0,0)) +
  scale_x_discrete(expand=expand_scale(add=c(.8, 1.75))) +
  scale_fill_manual(values=u_col[5:2], labels=c("High", "Moderate", "Low",  "None")) + 
  guides(fill = guide_legend(nrow = 2, byrow = T)) +
  labs(x = NULL, y = "Percentage", fill="Degree of \nbranching") +
  theme(legend.position=c(.5,.85),
        legend.direction="horizontal",
        axis.title = element_text(size = 12, face="bold"),
        axis.text.x = element_text(family = "Roboto Mono", size = 10),
        panel.grid = element_blank(),
        legend.title = element_text(size = 10, face="bold"),
        legend.text = element_text(size=10),
        legend.background = element_blank()) 

# * maximum depth ----
p4 <- ggplot(df.corals, aes(x = max_depth, y = iucn, fill = iucn, col=iucn)) +
  geom_density_ridges(alpha=0.2, scale=1.2) +
  theme_ridges() + 
  scale_fill_manual(values=u_cat, guide=FALSE) + 
  scale_color_manual(values=u_cat, guide=FALSE) + 
  scale_y_discrete(expand=expand_scale(add = c(0, 1.5))) +
  labs(y = NULL, x = "Maximum depth (m)") +
  theme_light(base_size = 15) +
  theme(axis.title = element_text(size = 12, face="bold"),
        axis.text.x = element_text(family = "Roboto Mono", size = 10),
        panel.grid = element_blank()) 

# * merge
svg(here("figs", "Fig_01_variables.svg"), w=8, h=8)
(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = "A")
dev.off()


# Fig 02 - Model performance (AUC of six models) --------------------------


