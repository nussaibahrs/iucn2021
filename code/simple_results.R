library(ggplot2)

modern <- read.csv(file.path("output", "modern_results.csv"))
pal <-  c("#1d3557", "#FF8B01", "#e63946")
df <- setNames(data.frame(prop.table(table(modern$pred, modern$iucn))),
         c("predicted", "red_list", "vals"))

df$cat <- "none"
df$cat[df$predicted=="T" & df$red_list=="T"] <- "risk-risk"
df$cat[df$predicted=="NT" & df$red_list=="T"] <- "risk-no risk"
df$cat[df$predicted=="T" & df$red_list=="NT"] <- "risk-no risk"
df$cat[df$predicted=="NT" & df$red_list=="NT"] <- "no risk-no risk"
df$cat <- factor(df$cat, levels=c("no risk-no risk", "risk-no risk", "risk-risk"))


library(extrafont)

# Lower left: Not at risk from fossil traits and risk in red list. 
# Lower right. Not at risk from fossil traits and not at risk in red list.
# Upper left: predicted at risk from fossil traits and risk in red list. 
# upper right: predicted at risk from fossil traits and no risk in red list. 

df$y <- NA
df$x <- NA

df[df$predicted=="NT" & df$red_list=="T", c("x", "y")]<- c(1,1)
df[df$predicted=="NT" & df$red_list=="NT", c("x", "y")]<- c(2,1)
df[df$predicted=="T" & df$red_list=="T", c("x", "y")]<- c(1,2)
df[df$predicted=="T" & df$red_list=="NT", c("x", "y")]<- c(2,2)

df$predicted <- plyr::mapvalues(df$predicted, c("NT", "T"), c("Not at risk", "At risk"))
df$red_list <- plyr::mapvalues(df$red_list, c("NT", "T"), c("Not at risk", "At risk"))


ggplot(df, aes(y=y, x=x)) +
  geom_point(aes(size=vals, col=cat)) +
  geom_vline(xintercept = 1.5) +
  geom_hline(yintercept = 1.5) +
  # Lower left: Not at risk from fossil traits and risk in red list. 
  # Lower right. Not at risk from fossil traits and not at risk in red list.
  # Upper left: predicted at risk from fossil traits and risk in red list. 
  # upper right: predicted at risk from fossil traits and no risk in red list. 
  scale_x_continuous(breaks=c(1,2), labels=c("At risk", "Not at risk")) +
  scale_y_continuous(breaks=c(1,2), labels=c("Not at risk", "At risk")) +
  coord_cartesian(xlim=c(0.5,2.5), ylim=c(0.5,2.5)) +
  geom_text(aes(label=paste0(round(vals,2)*100, "%")), col="white") +
  scale_size_continuous(range=c(12,30)) +
  scale_color_manual(values = pal) +
  labs(x="As per the IUCN", y="Predicted using fossils traits") +
  theme_minimal(base_family = "Roboto") +
  theme(legend.position = "none",
        axis.text.y = element_text(angle=90, hjust=0.5),
        panel.grid = element_blank(),
        axis.title = element_text(face=2),
        axis.title.y = element_text(angle=90))

ggsave("comparison.png", w=5, h=3, dpi=600)


# Waffle plot -------------------------------------------------------------
library(waffle)
library(patchwork)

pal <- c("darkgrey",  "#fcbf49", "#f77f00","#d62828")

modern <- na.omit(modern[,c("pred", "iucn")])
modern$comb <- paste(modern$iucn, modern$pred, sep="-")

comb <- table(modern$comb)


p1 <- waffle(comb[c(1,3,2,4)]/1, rows=12, size=0.5, 
             colors=pal, 
             title="How many corals are currently at risk of extinction?",
             xlab="1 square = 1 coral species") +
  theme(plot.title = element_text(size=12, face="bold"),
        legend.position = "none"
        ) +
  ylim(-10, 18) + xlim(0,75)

df$perc <- round(df$vals*100)
df$perc <- paste0(df$perc, "%")

xp <- c(10,25, 35, 45)
yp <- c(15, -5, -4)

fonts <- 4


p1 +annotate("segment", x=xp[1],xend =xp[1] , y=yp[1], yend=yp[1]-3, linetype="dashed",
             col=pal[1]) +
  annotate("text", x=xp[1], y=16, label="Not threatened", size=fonts, fontface=1, col=pal[1],
           hjust=1) +
  
  
  annotate("segment", x=xp[2], xend=xp[2], y=0, yend=yp[2], linetype="dashed", col=pal[2]) +
  annotate("text", x=xp[2], y=yp[2]-3, label=paste0('atop(bold("',df$perc[3],'")~" could be at"~bold("lower risk"))'), parse=TRUE, 
           size=fonts, hjust=1, col=pal[2]) +
  annotate("text", x=xp[2], y=yp[2]-fonts, label="than currently thought",
           size=fonts, hjust=1, col=pal[2]) +
  
  annotate("segment", x=xp[3], xend=xp[3], y=0, yend=yp[3], linetype="dashed", col=pal[3]) +
  annotate("text", x=xp[3], y=yp[3]-1.5, label=paste0('atop(bold("',df$perc[2],'")~"',"that  are at low risk as per the IUCN",'")'), 
           size=fonts, hjust=0, parse=TRUE, col=pal[3]) +
  annotate("text", x=xp[3], y=yp[3]-fonts, label='atop("were predicted as being"~bold("threatened"))', 
           size=fonts, hjust=0,col=pal[3], parse=TRUE) +
  annotate("segment", x=xp[4], xend=xp[4], y=yp[1], yend=yp[1]-3, linetype="dashed",
           col=pal[4]) +
  annotate("text", x=xp[4], y=yp[1]+1, label=paste0('atop("Only"~bold("',df$perc[4],'")~"of at risk corals")'), 
           size=fonts, hjust=0, col=pal[4], parse=2) +
  annotate("text", x=xp[4], y=yp[1]+1, label="matched the model predictions", 
           size=fonts, hjust=0, col=pal[4])

ggsave("waffle_coral.svg", w=10, h=5)


