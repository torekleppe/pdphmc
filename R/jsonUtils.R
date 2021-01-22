


.json.matrix <- function(l){
  if(length(l)>0){
    for(i in 1:length(l)){
      if(is.matrix(l[[i]])){
        l[[i]] <- list(nrow=dim(l[[i]])[1],
                       vals=as.vector(l[[i]])) 
      }
    }
  }
  return(l)
}

