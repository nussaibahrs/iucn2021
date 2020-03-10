library(here)
library(tidyverse)
library(ggridges)
library(patchwork)

#palette
u_col <- ggsci::pal_uchicago("default")(9)[c(2,5,4,3,1)]
#pie(1:5, col=u_col)

### load data ####
df.corals <- read.csv(here("data", "traits_iucn.csv"), stringsAsFactors = FALSE) %>%
  na.omit() %>%
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

theme_set(theme_light(base_size = 15))
# % of observations per status
p0 <- df.corals %>% group_by(iucn) %>% 
  tally() %>%
  ungroup() %>%
  mutate(n = n/sum(n) * 100) %>%
  ggplot() +
  geom_bar(aes(x=iucn, y=n, fill=iucn), stat="identity") +
  coord_cartesian(expand=FALSE, ylim=c(0,40)) +
  labs(x="IUCN Red List Category", y="Percentage (%)") +
  scale_fill_manual(values=c(u_col, "black"), guide=FALSE) +
  theme(axis.title = element_text(face="bold"),
        panel.grid = element_blank()) +
  geom_vline(xintercept=c(1.5, 3.5), linetype="dashed")

p0 <- p0 + annotate("text", x=c(1,2.5, 5), y=38, label=c("Data-deficient", "Non-threatened", "Threatened"), fontface=2)
ggsave(here("figs", "Fig_SI_species_per_category.svg"), p0, w=8, h=8)

#CR To EN
df.corals[df.corals$iucn == "CR",]$iucn <- "EN"
df.corals$iucn <- droplevels(df.corals$iucn)

table(df.corals$iucn)

##### status vs range
world_avg <- df.corals %>%
  summarize(avg = median(range, na.rm = T)) %>%
  pull(avg)

levels(df.corals$iucn) <- c("DD", "N", "N", "T", "T")

reg_avg <- df.corals %>% group_by(iucn) %>%
  summarise(n=median(range, na.rm=TRUE))

p1 <- ggplot(df.corals, aes(x = iucn, y = range, color = iucn)) +
  coord_flip() +
  geom_jitter(size = 2, alpha = 0.1, width = 0.2) +
  scale_color_manual(values=u_col[c(1,2,5)], guide=FALSE) + 
  labs(x = NULL, y = "Range (km)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face="bold"),
        axis.text.x = element_text(family = "Roboto Mono", size = 10),
        panel.grid = element_blank()) + 
  geom_hline(aes(yintercept = world_avg), color = "gray70", size = 0.6) +
  geom_segment(data=reg_avg, aes(x = iucn, xend = iucn,
                                 y = world_avg, yend = n, col=iucn),
               size = 0.8, inherit.aes = FALSE) +
  stat_summary(fun.y = median, geom = "point", size = 5) 

# status vs corallite diameter
world_avg2 <- df.corals %>%
  summarize(avg = median(corallite, na.rm = T)) %>%
  pull(avg)

reg_avg2 <- df.corals %>% group_by(iucn) %>%
  summarise(n=median(corallite, na.rm=TRUE))

p2 <- ggplot(na.omit(df.corals), aes(x = iucn, y = corallite, color = iucn)) +
  coord_flip() +
  geom_jitter(size = 2, alpha = 0.1, width = 0.2) +
  scale_color_manual(values=u_col[c(1,2,5)], guide=FALSE) + 
  labs(x = NULL, y = "Corallite Diameter (mm)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face="bold"),
        axis.text.x = element_text(family = "Roboto Mono", size = 10),
        panel.grid = element_blank()) + 
  geom_hline(aes(yintercept = world_avg2), color = "gray70", size = 0.6) +
  geom_segment(data=reg_avg2, aes(x = iucn, xend = iucn,
                                 y = world_avg2, yend = n, col=iucn),
               size = 0.8, inherit.aes = FALSE) +
  stat_summary(fun.y = median, geom = "point", size = 5) +
  scale_y_log10(breaks=c(0.1, 0.5,1,5, 25, 50,100, 200), expand = c(0.005, 0.005),
                labels=c(0,0.5, 1,5,25,50,100,200))


#branching
p3 <- df.corals  %>% 
  group_by(iucn, branching) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count)) %>%
  ggplot(aes(iucn, fill=branching, y=perc*100)) +
  geom_bar(stat="identity", width = 0.6) +
  coord_flip() +
  scale_y_continuous(limits=c(0, 100), expand=c(0,0)) +
  scale_x_discrete(expand=expand_scale(add=c(.8, 1.75))) +
  scale_fill_manual(values=u_col[5:2], labels=c("High", "Moderate", "Low",  "None")) + 
  guides(fill = guide_legend(nrow = 2, byrow = T)) +
  labs(x = NULL, y = "Percentage", fill="Degree of \nbranching") +
  theme(legend.position=c(.5,.9),
        legend.direction="horizontal",
        axis.title = element_text(size = 12, face="bold"),
        axis.text.x = element_text(family = "Roboto Mono", size = 10),
        panel.grid = element_blank(),
        legend.title = element_text(size = 10, face="bold"),
        legend.text = element_text(size=10),
        legend.background = element_blank()) 

#### maximum depth
p4 <- ggplot(df.corals, aes(x = max_depth, y = iucn, fill = iucn, col=iucn)) +
  geom_density_ridges(alpha=0.2, scale=1.2) +
  theme_ridges() + 
  scale_fill_manual(values=u_col[c(1,2,5)], guide=FALSE) + 
  scale_color_manual(values=u_col[c(1,2,5)], guide=FALSE) + 
  scale_y_discrete(expand=expand_scale(add = c(0, 1.5))) +
  labs(y = NULL, x = "Maximum depth (m)") +
  theme_light(base_size = 15) +
  theme(axis.title = element_text(size = 12, face="bold"),
        axis.text.x = element_text(family = "Roboto Mono", size = 10),
        panel.grid = element_blank()) 

# merge
svg(here("figs", "Fig_01_variables.svg"), w=8, h=8)
(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = "A")
dev.off()
