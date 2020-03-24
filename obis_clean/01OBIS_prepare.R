load("/home/nussaibah/Dropbox/nuss/iucn_traits/data/original/2019-11-13_obis_scleractinia.RData")
scler <- obis.corals
nrow(scler)

colnames(scler)[grep("decimal", colnames(scler))] <- c("latitude", "longitude")
scler$tname <- scler$scientificName
coral <- subset(scler, is.na(latitude)==F)

# Isolate genus name
gen <- sub(" .*", "", coral$tname)
gen <- table(gen)

# Check Maeandrina, Montastraea, Lophohelia, Tubastrea, Tubinaria
coral$tname <- sub("Maeandrina", "Meandrina", coral$tname)
coral$tname <- sub("Montastraea", "Montastrea", coral$tname)
coral$tname <- sub("Lophohelia", "Lophelia", coral$tname)
coral$tname <- sub("Caulastraea", "Caulastrea", coral$tname)
coral$tname <- sub("Tubinaria", "Turbinaria", coral$tname)

# Correct misplace spaces
coral$tname[coral$tname=="Platygyra daedalea "] <- "Platygyra daedalea"

# View(coral[grep("Mussismilia h", coral$tname),])
# 
# # Revisions based on Wallace (1999) and Veron (2000) and WoRMS (2010)
# syno <- read.csv("c:/daten/PBDB/coral-synonyms.csv", sep=";", header=T) # carefully checked until 10 occurrences
# 
# 
# # Takes 1 minute
# for (i in 1:nrow(syno)) coral$tname <- sub(syno$OBIS.name[i], syno$valid.name[i], coral$tname, fixed=T)

sp1 <- factor(coral$tname)

# Get information on symbiotic status and clean taxonomy
cor.sp <- readxl::read_excel("/home/nussaibah/Dropbox/nuss/iucn_traits/data/coral_species_new.xlsx")
cor.sp <- subset(cor.sp, is.na(valid))
cor.sp <- unique(cor.sp)
Species <- cor.sp$scientificName
cor.sp <- cbind(cor.sp, Species)


# Problem with NMNH Invertebrate Zoology Collections: Fossils included!
# Delete NMNH Invertebrate Zoology Collections that do nat have physical data
coral <- subset(coral, resname!="NMNH Invertebrate Zoology Collections" | temperature > 0)

# which species in coral master-list are not in OBIS
# Species[Species %in% coral$tname==F & cor.sp$symbio!="az"]
# which species in OBIS are not in coral masterlist
sort(table(coral$tname[coral$tname %in% Species==F])) # by occurrences
table(coral$tname[coral$tname %in% Species==F]) # alphabetical

# Drop unknown taxa 
# Note many sp identifications in OBIS data
coral2 <- merge(coral, cor.sp, by.x="tname", by.y="Species")

sp2 <- factor(coral2$tname)

Species[Species %in% sp2 ==F]

x <- as.data.frame(sort(table(sp1[sp1 %in% sp2 ==F])))
write.table(x, file="c:/daten/R/Pleistocene/obis_not_included.csv", sep=";")


# Look for individual species names
View(coral[grep("harttii", coral$tname),])

View(coral2[grep("Mussismilia h", coral2$tname),])
View(coral[grep("Mussismilia h", coral$tname),])
nrow(coral2)

# limit to reef corals
# coral <- subset(coral, symbio=="") 
# limit to az corals
# coral <- subset(coral, symbio=="az")    


# Consider removing all rounded coordinates (currently not implemented)
# coral2 <- subset(coral2, latitude == round(latitude) & longitude== round(longitude)))

# make bathymetry numeric
coral2$depth <- as.character(coral2$depth)
coral2$depth <- as.numeric(coral2$depth)
# replace greater than 10000 m with NA
coral2$depth[coral2$depth>10000] <- NA

# Solve Problem with bathymetry
# American data might be in feet rather than meters
#    St. Croix, USVI Benthic Composition and Monitoring Data (2002 - Present)
x <- which(coral2$resname==
             "St. Croix, USVI Benthic Composition and Monitoring Data (2002 - Present)")
coral2$depth[x] <- coral2$depth[x]*0.3048       

#    Biogeographic Characterization of Benthic Habitat Communities within the Flower Garden Banks National Marine Sanctuary (2006 - Present)
# Maximum depth 35 m according to website
x <- which(coral2$resname==
             "Biogeographic Characterization of Benthic Habitat Communities within the Flower Garden Banks National Marine Sanctuary (2006 - Present)")
coral2$depth[x] <- coral2$depth[x]*0.3048

#    St. John, USVI Benthic Composition and Monitoring Data (2002 - Present)
# Maximum depth 20 m according to Harborne et al. (2006, Ecology)
x <- which(coral2$resname==
             "St. John, USVI Benthic Composition and Monitoring Data (2002 - Present)")
coral2$depth[x] <- coral2$depth[x]*0.3048

# La Parguera, Puerto Rico Benthic Composition and Monitoring Data (2002 - Present)
# Maximum depth 100 ft according to website
x <- which(coral2$resname==
             "La Parguera, Puerto Rico Benthic Composition and Monitoring Data (2002 - Present)")
coral2$depth[x] <- coral2$depth[x]*0.3048

# round coordinates to 4 digits
coral2$latitude <- round(coral2$latitude, 4)
coral2$longitude <- round(coral2$longitude, 4)  

# exclude high latitude collection with strange data (lat, long mixup?
coral2 <- subset(coral2, latitude!= -55.5333)




# Export this clean file
# write.table(coral, file="C:/Daten/GBIF/OBIS2_clean.csv", row.names=F, sep=";")
# write.table(coral2, file="C:/Daten/GBIF/OBIS3_clean.csv", row.names=F, sep=";")
# write.table(coral2, file="C:/Daten/GBIF/OBIS4_clean.csv", row.names=F, sep=";")
write.table(coral2, file="C:/Daten/GBIF/OBIS5_clean.csv", row.names=F, sep=";")

#############

# Calculation of mean physical parameters
# omit deep stuff
# x <- which(coral2$symbio!="az" & coral2$depth>100)
#   coral3 <- coral2[-x,] 

# Use all data
coral3 <- coral2

# N
N <- table(coral3$tname)
# Bathymetry
min.dep <- tapply(coral3$depth, coral3$tname, min, na.rm=T)
max.dep <- tapply(coral3$depth, coral3$tname, max, na.rm=T)
median.dep <- tapply(coral3$depth, coral3$tname, median, na.rm=T)
# Temperature
min.t <- tapply(coral3$temperature, coral3$tname, min, na.rm=T)
max.t <- tapply(coral3$temperature, coral3$tname, max, na.rm=T) 
median.t <- tapply(coral3$temperature, coral3$tname, median, na.rm=T) 
# Salinity
min.sal <- tapply(coral3$salinity, coral3$tname, min, na.rm=T)
max.sal <- tapply(coral3$salinity, coral3$tname, max, na.rm=T)
median.sal <- tapply(coral3$salinity, coral3$tname, median, na.rm=T)
# Nitrate
min.nit <- tapply(coral3$nitrate, coral3$tname, min, na.rm=T)
max.nit <- tapply(coral3$nitrate, coral3$tname, max, na.rm=T)
median.nit <- tapply(coral3$nitrate, coral3$tname, median, na.rm=T)
# Phosphate
min.pho <- tapply(coral3$phosphate, coral3$tname, min, na.rm=T)
max.pho <- tapply(coral3$phosphate, coral3$tname, max, na.rm=T)
median.pho <- tapply(coral3$phosphate, coral3$tname, median, na.rm=T)
# Latitude
max.lat.n <- tapply(coral3$latitude, coral3$tname, max, na.rm=T)
min.lat <- tapply(abs(coral3$latitude), coral3$tname, min, na.rm=T) # closest to equator
median.lat <- tapply(abs(coral3$latitude), coral3$tname, median, na.rm=T)
max.lat.s <- tapply(coral3$latitude, coral3$tname, min, na.rm=T)

res <- cbind(N, min.dep, max.dep, median.dep, min.t, max.t, median.t, 
             min.sal, max.sal, median.sal, min.nit, max.nit, median.nit, 
             min.pho, max.pho, median.pho,
             max.lat.n, min.lat, max.lat.s, median.lat)

res2 <- as.data.frame(cbind(names(N), res))
colnames(res2) <- c("spec", colnames(res))

sp <- subset(cor.sp, select=c(symbio, Species))
res3 <- merge(res2, sp, by.x="spec", by.y="Species", all.x=TRUE)

# edit(res2)

# Link to higher taxa
sepk <- read.csv("c:/daten/PBDB/Scleract.csv", sep=";", header=T)

# conect to higher taxa for taxonomy
gen <- unlist(strsplit(as.character(res3$spec), " "))
i.gen <- seq(1, 2*length(rownames(res3)), by=2)
gen <- gen[i.gen]
tax <- subset(sepk, GENUS %in% gen)
tax <- subset(tax, select=c("INTEGRAT", "GROWTH", "MORPH",
                            "CLADE", "CONF_CLADE", "FAMILY_NAM", "GENUS"))
sp1 <- as.data.frame(cbind(gen, as.character(res3$spec)))
sp.tax <- merge(sp1, tax, by.x="gen", by.y="GENUS")
# colnames(sp.tax) <- c("gen", "sp", "fam")

res4 <- merge(sp.tax, res3, by.x="V2", by.y="spec")

# for all data
# old download
# write.table(res2, file="c:/daten/r/output/Obis4_taxa.csv", sep=";", row.names=F)
# new donwload
write.table(res4, file="c:/daten/r/Pleistocene/Obis5_taxa.csv", sep=";", row.names=F)

# Limit to ap and z corals for Study with JP and Morana
res5 <- subset(res4, symbio!="az")

NAs <- res5 == "Inf" | res5 == "-Inf"
res5[NAs] <- NA

# add sex from Kerr et al. (2011)
sex <- read.table(file="c:/daten/PBDB/coral-sex.csv", sep=";", header=T)
res5 <- merge(res5, sex, by.x="V2", by.y="species", all.x=T)

# Add paleo information
PBDB <- read.table(file="C:/Daten/PBDB/Coral_occs_July12.csv", header=TRUE, sep=",", quote="\"")
PBDB$gensp <- paste(PBDB$occurrence.genus_name, PBDB$occurrence.species_name)
oldest <- as.data.frame(tapply(PBDB$ma_max, PBDB$gensp, max))
colnames(oldest) <- "ma.oldest.record"

res6 <- merge(res5, oldest, by.x="V2", by.y="row.names", all.x=T)


write.table(res6, file="c:/daten/r/Pleistocene/Coral_species.csv", sep=";", row.names=F)

# check high latitude occurrences of non-az corals
which(max.lat>40 & res3$symbio!="az")
which(max.dep>80 & res3$symbio!="az")

View(subset(coral3, symbio!="az" & depth>80))  

