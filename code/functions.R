# RTD ---------------------------------------------------------------------

clean_rtd <- function(rtd){
  n <- grep("Montastraea", rtd$genus_name)
  rtd$genus_name[n] <- "Montastrea"
  
  rtd <- rtd[!rtd$species_name == "sp.",]
  
  rtd$sp <- paste(rtd$genus_name, rtd$species_name)
  rtd$sp <- gsub("\\(\\?\\)|cf\\. |sp\\. nov\\. ", "", rtd$sp)
  rtd$sp <- gsub("  ", " ", rtd$sp)
  
  return(rtd)
}


#rescaling
scaled.new <- function (newdata, olddata) {(newdata-min(olddata, na.rm = TRUE))/(max(olddata, na.rm = TRUE)-min(olddata, na.rm = TRUE))}
unscale <- function(newdata, olddata){min(olddata, na.rm = TRUE) + newdata * (max(olddata, na.rm = TRUE) - min(olddata, na.rm=TRUE))}
scale01 <- function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}

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
    r <- range(pdp$`_yhat_`)
    
    p <- ggplot(pdp, aes(x=`_x_`, y=`_yhat_`)) + geom_line(col=col) +
      labs(x=lab, y="Predicted Value") +
      scale_y_continuous(expand=expand_scale(c(0, 0.1)))+
      geom_linerange(aes(x=`_x_`, ymin=0, ymax=r[2]/20), col="darkgrey") +
      theme_light(base_size = 15) +
      theme(axis.title = element_text(size = 12, face="bold"),
            axis.text = element_text(size = 10),
            panel.grid = element_blank(),
            legend.position = "none")
  }
  
  return(p)
}

