#obtained from http://www.marinespecies.org/aphia.php?p=webservice&type=r

#install the required packages (comment if needed)
# install.packages("jsonlite", repos="http://cran.r-project.org")
# install.packages("httr")

#Use the libraries
library(jsonlite) #https://cran.r-project.org/web/packages/jsonlite/
library(httr)

worms_api <- function (namesToMatch, verbose=FALSE){
  #Convert the namesToMatch to a valid REST-url
  urlNamesPart <- ""
  
  for (index in 1:length(namesToMatch)) {
    urlNamesPart <- sprintf("%s&scientificnames[]=%s", urlNamesPart, namesToMatch[index]);
  }
  
  #The url can contain special characters that need to be converted
  urlNamesPart <- URLencode(urlNamesPart)
  
  #The dyanmic build of the URL causes an obsolete '&' at the beginning of the string, so remove it
  urlNamesPart <- substring(urlNamesPart, 2)
  
  #Build the final REST-url
  url <- sprintf("http://www.marinespecies.org/rest/AphiaRecordsByMatchNames?%s", urlNamesPart);
  
  #Get the actual data from the URL
  matches <- fromJSON(url)
  
  valid <- c()
  #Handle the data (each requested name has an list of results)
  for (matchesindex in 1:length(namesToMatch)) {
    #Get the results for the current index
    currentResultList = matches[[matchesindex]]
    currentResultList = currentResultList[!is.na(currentResultList$valid_AphiaID),] #removes na before giving valid taxa
    
    #Handle empty data due to no matches found
    if (length(currentResultList) <1) {
      numberOfResults <- 0
    } else {
      #Get the number of list entries for the first column
      numberOfResults <- length(currentResultList[[1]])
    }
    
    if (verbose == TRUE){
      print(sprintf("%d Result(s) found for %s", numberOfResults, namesToMatch[matchesindex]))      
    }
    
    
    if (numberOfResults > 0) {
      for (listentry in 1:numberOfResults) {
        if (verbose == TRUE){
          print(sprintf("ID: %d, SCIENTIFICNAME: %s, MATCH_TYPE: %s, VALIDNAME: %s",
                        currentResultList[["AphiaID"]][listentry],
                        currentResultList[["scientificname"]][listentry],
                        currentResultList[["match_type"]][listentry],
                        currentResultList[["valid_name"]][listentry]
          ))}
        valid[matchesindex] <- currentResultList[["valid_name"]][listentry] 
      } 
    }else valid[matchesindex] <- NA 
  }
  return(valid)
}

# correcting the point recontstruction problem, wrapper around the point reconstruction funciton
match_taxa <- function(namesToMatch, chunk=50, verbose=TRUE){
  # number of coordinates
  coordNum <- length(namesToMatch)
  
  # do only when the number of coordinates is large enough
  if(coordNum>chunk){
    # batch number
    moreIndex <- rep(1:ceiling(coordNum/chunk), each=chunk)
    index<-moreIndex[1:coordNum]
    
    # new container
    newCoords <- matrix(NA, ncol=1, nrow=coordNum)
    
    # the number of batches 
    maxIndex <- max(index)
    
    # iterate - for() is easier, performance impact unlikely
    for(i in 1:maxIndex){
      cat("\r", i , "out of", maxIndex)
      # index of batch number
      bIndex <- index==i
      
      # current batch 
      current <- namesToMatch[bIndex]
      
      # do the reconstruction
      tryCatch({
        iterRes <- worms_api(current, verbose=verbose)
      },
      error=function(cond){
        stop("Query URL is too long. Round coordinates or decrease chunk size.")
      }
      ) 
      
      
      # store
      newCoords[bIndex] <- iterRes
    }
    
    # save some time by skipping this
  } else{
    tryCatch({
      newCoords <- worms_api(namesToMatch, verbose=verbose)
    }, error=function(cond){
      stop("Query URL is too long. Round coordinates or decrease chunk size.")
    }
    )
    
  }
  
  return(newCoords)
  
}