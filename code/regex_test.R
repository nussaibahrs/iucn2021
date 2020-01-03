library(here)
library(tidyverse)
library(magrittr)

dat <- read.csv(here("data", "pleist_sp.csv"), stringsAsFactors = FALSE)

url <- "https://fossils.its.uiowa.edu/database/corals/combined/%s.htm"

corallite <- c()
shape <- c()

for (i in 1:nrow(dat)){
  cat("\r", i, "out of", nrow(dat))
  
  sp <- gsub("\\.", "_", paste(dat[i,c("genus", "species")], collapse = ""))
  
  url2 <- sprintf(url, sp)
  
  x <- tryCatch({readLines(url(url2))}, error = function(e){cat("no"); return(NULL)})
  
  if(!is.null(x)){
    n <- grep("Corallite Diameter", x)
    
    corallite[i] <- ifelse(length(n) > 0, sub(".*\\b(\\d+.\\d+).*", "\\1", x[n]), NA) #extract last number
    
    n <- grep("Colony Shape", x)
    
    shape[i] <- ifelse(length(n) > 0, sub("^.+Colony Shape:</B>\\s*([a-z, ]+).*<BR>", "\\1", x[n]), NA)
  }
}

dat <- cbind(dat[,1:3], corallite, shape) %T>% write.csv(here("data", "pleist_sp.csv"), row.names = FALSE)

#### check with ctbd
ctdb <- read.csv(here("data", "ctdb_resolved.csv"), stringsAsFactors = FALSE)


dat$max_depth <- NA

for (i in 1:nrow(dat)){
  if (is.na(dat$corallite[i]) & dat$name[i] %in% ctdb$valid){
    temp <- ctdb %>% filter(valid == dat$name[i] & trait_name == "Corallite width maximum") %>%
      pull(value)
    dat$corallite[i] <- ifelse(length(temp) > 0, temp, NA)
  }
  
  if (is.na(dat$shape[i]) & dat$name[i] %in% ctdb$valid){
    temp <- ctdb %>% filter(valid == dat$name[i] & trait_name == "Growth form typical") %>%
      pull(value)
    
    dat$corallite[i] <- ifelse(length(temp) > 0, temp, NA)
  }
  
  if (dat$name[i] %in% ctdb$valid){
    temp <- ctdb %>% filter(valid == dat$name[i] & trait_name == "Depth lower") %>%
      pull(value)
    dat$max_depth[i] <- ifelse(length(temp) > 0, temp, NA)
  }
  
}

write.csv(dat, here("data", "pleist_sp.csv"), row.names = FALSE)
  
x <- read.csv(here("data", "pleist_sp.csv")) 
x <- x[,-ncol(x)]
  