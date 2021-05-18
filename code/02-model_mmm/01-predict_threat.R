# Set up environment ------------------------------------------------------
# helper packages
library(tidyverse)
library(magrittr)
library(grid)
library(icosa)
library(raster)
library(sp)
library(extrafont)

# setting up machine learning models
library(h2o) # version 3.28.0.2

#palette
u_col <- ggsci::pal_uchicago("default")(9)[c(2,5,4,3,1)]

# Data and h2o setup--------------------------------------------------------------------
source(file.path("code", "02-modern_corals", "00-load_data.R"))

# Load model --------------------------------------------------------------
auc.df <- read.csv(file.path("output", "Table_S_model_performance_fmm.csv")) 
leaderboard <- read.csv(file.path("output", "model_fmm_8h", "leaderboard.csv"), stringsAsFactors = FALSE)

w=1 #winning model
win <- auc.df$n[w]
aml_leader <- h2o.loadModel(file.path("output", "model_fmm_8h", auc.df$fname[w]))

df <- df.corals %>%
  #treat as binary
  mutate(
    status =iucn,
    #reclassify iucn status
    iucn = case_when(iucn == "Least Concern" ~ "NT",
                     iucn == "Vulnerable" ~ "T",
                     iucn == "Near Threatened"   ~ "NT",
                     iucn == "Endangered"   ~ "T",
                     iucn == "Critically Endangered"~ "T"),
    branching = factor(branching, levels=c("NB", "LB", "MB", "HB")))

#comparing iucn
res <- h2o.predict(aml_leader, as.h2o(df))
res <- as.data.frame(res)$p1

df$pred <- ifelse(res < auc.df$cutoff[w], "NT", "T")

write.csv(df, file.path("output", "modern_results.csv"), row.names=FALSE)
prop.table(table(df$pred))
prop.table(table(df$iucn, df$pred), 1)

# comparison between extinction risk and conservation status for
  #* Least Concern with high extinction risk
  #* Endangered with low extinction risk

LC <- df %>%  filter(status == "Least Concern" & pred=="T")
EN <-  df %>%  filter(status %in% c("Endangered", "Critically Endangered") & pred=="NT")

#non data deficient
ndd <- df %>% filter(status !="Data Deficient") 

nrow(ndd[ndd$pred == ndd$iucn,])/nrow(ndd)

mismatch <- ndd[ndd$pred != ndd$iucn,]
table(mismatch$status)
table(ndd$status)

mismatch <- data.frame(table(mismatch$status)/table(ndd$status))

colnames(mismatch)[1] <- "iucn"

mismatch <- mismatch %>%  mutate(iucn = case_when(iucn == "Least Concern" ~ "LC",         
                                                  iucn == "Vulnerable" ~ "VU",           
                                                  iucn == "Near Threatened"   ~ "NT",    
                                                  iucn == "Endangered"   ~ "EN",         
                                                  iucn == "Critically Endangered"~ "CR",
                                                  iucn == "Data Deficient"~ "DD"),
                                 iucn = factor(iucn, levels = c("DD", "LC", "NT", "VU", "EN", "CR"))
) 

# fossil data
df.corals <- read.csv(file.path("data", "pbdb_pleist.csv"), stringsAsFactors = FALSE)
df.corals[df.corals == ""] <- NA

df.fossil <- df.corals %>%
  #na.omit() %>%
  #treat as binary
  mutate(
    #reclassify iucn status 
    extinct = case_when(extinct == "extant" ~ "NT",         
                        TRUE~"T"),
    extinct=factor(extinct, levels = c("NT","T")),
    branching = factor(branching, levels=c("NB", "LB", "MB", "HB"))) %>% 
  #omit species name
  dplyr::select(corallite, branching, budding, family, integration, extinct) %>% 
  ungroup()


res <- h2o.predict(aml_leader, as.h2o(df.fossil[,x]))
res <- as.data.frame(res)$p1

df.fossil$pred <- ifelse(res < auc.df$cutoff[w], "NT", "T")

predicted_ext <- c(fossil=prop.table(table(df.fossil[df.fossil$extinct=="T",]$pred))[2],
prop.table(table(df$status, df$pred), 1)[,2])
names(predicted_ext)[1] <- "Extinct taxa"

predicted_ext <- data.frame(cat=factor(names(predicted_ext), 
                                       levels=c("Extinct taxa", "Least Concern", "Near Threatened", "Vulnerable", "Endangered", "Critically Endangered", "Data Deficient")), 
                                       prop=predicted_ext)

theme_set(theme_light(base_size = 15))

ggplot(predicted_ext, aes(x=cat, y=prop)) +
  geom_point(size=3, col=u_col[5])+
  coord_cartesian(ylim=c(0,1)) +
  labs(x="", y="Proportion of coral species\nat risk of extinction") +
  theme(axis.title = element_text(face="bold"),
        axis.text.x = element_text(angle=45, hjust=1, family="Roboto"),
        panel.grid = element_blank()) +
  geom_vline(xintercept = 1.5, col=u_col[1]) +
  geom_vline(xintercept = 6.5, col=u_col[1], linetype="dashed")

ggsave(file.path("figs", "Fig_06_mismatch.svg"), w=5, h=5)

prop.table(table(ndd$pred))

# DD species ------------------------------------------------
dd <- df %>% filter(status =="Data Deficient") 

prop.table(table(dd$pred))

# * Prediction of threat status -------------------------------------------

res <- h2o.predict(aml_leader, as.h2o(dd))
res <- as.data.frame(res)$p1

# calculate threat status based on threshold cutoff
dd$pred <- ifelse(res < leaderboard$cutoff[win], "NT", "T")

table(dd$pred)

table(df$iucn) %>% prop.table()
table(dd$pred) %>% prop.table()
df$iucn2 <- df$iucn
table(c(df$iucn, dd$pred)) %>% prop.table()

# * Map of DD distribution risk -------------------------------------------
gr <- hexagrid(16, sp=TRUE) #2.5 degrees
world <- rworldmap::getMap()
world <- cleangeo::clgeo_Clean(world)
world <- raster::aggregate(world, dissolve=TRUE)


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


dd_occ <- read.csv(file.path("data", "2019-11-05_obis_scleractinia.csv"), stringsAsFactors = FALSE) %>%
  filter(scientificName %in% df$valid.name) %>% #select only occurrences for dd species
  distinct(scientificName, decimalLatitude, decimalLongitude) %>%
  left_join(df %>% dplyr::select(valid.name, pred, status), by=c("scientificName" = "valid.name"))

dd_occ$cell <- locate(gr, dd_occ[, c(2,3)])
dd_occ$region <- ifelse(dd_occ$cell %in% carribean$faces, 
                        "carribean", "indo-pacific")

dd_occ$iucn <- plyr::mapvalues(dd_occ$status, unique(dd_occ$status), c("NT", "NT", "T", "T", "DD", "T"))

#results
cat("Modern Corals");table(dd_occ$region, dd_occ$pred) %>% prop.table(1)
cat("Modern Corals, iucn");table(dd_occ$region, dd_occ$iucn) %>% prop.table(1)

me_iucn <- tapply(INDEX=dd_occ$cell, X=dd_occ$iucn, function (x) length(x[x == "T"])/length(x[x=="NT"| x == "T"]))
me_dd <- tapply(INDEX=dd_occ$cell, X=dd_occ$pred, function (x) length(x[x == "T"])/length(x[x=="NT"| x == "T"]))


##### plot
txt <- 1.2
asp=1.2

cuts <- seq(0,1, 0.2) #c(seq(0.0, 0.4, 0.1), 1)
map_col <- colorRampPalette(u_col[c(2,4,5)])(length(cuts)-1)


svg(file.path("figs", "Fig_05_map_dd_risk.svg"), w=7, h=5)
layout(matrix(c(1,2, 3), nrow=3, byrow=TRUE),
       heights = c(1,1,0.3))
par(mar=c(1,3,2,1))

#IUCN
plot(NULL, xlim=c(-180, 180), ylim=c(-30,30), axes=FALSE, 
     xlab="", ylab="",
     asp=asp, xaxs="i")

for (i in 1:(length(cuts)-1)){
  n <- which(me_iucn >= cuts[i] & me_iucn < cuts[i+1])
  temp <- names(n)
  
  plot(gr[temp], col=map_col[i], border="white", add=TRUE)
}

plot(world, col=scales::alpha("grey95", 0.5), border="grey80", add=TRUE)
mtext ("A", side=3, line=0, adj=-0.04, cex=txt)
box(col="darkgrey")

#Predicted
plot(NULL, xlim=c(-180, 180), ylim=c(-30,30), axes=FALSE, 
     xlab="", ylab="",
     asp=asp, xaxs="i")

for (i in 1:(length(cuts)-1)){
  n <- which(me_dd >= cuts[i] & me_dd < cuts[i+1])
  temp <- names(n)
  
  plot(gr[temp], col=map_col[i], border="white", add=TRUE)
}

plot(world, col=scales::alpha("grey95", 0.5), border="grey80", add=TRUE)

mtext ("B", side=3, line=0, adj=-0.04, cex=txt)
box(col="darkgrey")


#add legend
par(xpd=TRUE)
plot(0,0, type="n", axes=FALSE)
cuts_labs <-cuts[-c(6,7)]*100 #turn to percentage
legend("center", c(paste0(cuts_labs[1:(length(cuts)-3)], "-", 
                          cuts_labs[2:(length(cuts)-2)], "%"), "> 60%"), fill=map_col, horiz=TRUE, border=NA, bty="n",
       title = expression(bold("Percentage of coral species under threat")),
       cex=1.2)


dev.off()
