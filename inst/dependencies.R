# Attachments ----
to_install <- c("dplyr","ggthemes","hmeasure","magrittr","patchwork","pROC","tidyverse","ggsci","caret","DALEX","ggpubr","grid","icosa","ingredients","pdp","plyr","raster","reshape2","cleangeo","extrafont","rworldmap","scales","sp","fields","rgdal","rgeos","robis","worrms","ggplot2")
  for (i in to_install) {
    message(paste("looking for ", i))
    if (!requireNamespace(i)) {
      message(paste("     installing", i))
      install.packages(i)
    }
  }

# The following two commands remove any previously installed H2O packages for R.
if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }
if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }

# Next, we download packages that H2O depends on.
pkgs <- c("RCurl","jsonlite")
for (pkg in pkgs) {
  if (! (pkg %in% rownames(installed.packages()))) { install.packages(pkg) }
}

# Now we download, install and initialize the H2O package for R.
install.packages("h2o", type="source", repos="https://h2o-release.s3.amazonaws.com/h2o/rel-yu/2/R")