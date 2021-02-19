


.blocks <- c('INCLUDE_BLOCK','DATA_BLOCK','SETUP_BLOCK','VARIABLES_BLOCK','MODEL_BLOCK','ATEVENT_BLOCK')

.data.types <- function(){
  keyword <- c('DATA_SCALAR','DATA_VECTOR','DATA_MATRIX',
               'DATA_ISCALAR','DATA_IVECTOR','DATA_IMATRIX')
  cppType <- c('double','Eigen::VectorXd','Eigen::MatrixXd',
               'int','Eigen::VectorXi','Eigen::MatrixXi')
  numArgs <- c(1,1,1,
               1,1,1)
  storageType <- c("scalar","vector","matrix",
                   "scalar","vector","matrix")
  return(data.frame(keyword,cppType,numArgs,storageType))
}

.parameter.types <- function(){
  keyword <- c('PARAMETER_SCALAR','PARAMETER_VECTOR','PARAMETER_MATRIX')
  cppType <- c('var','Eigen::Matrix<var,Eigen::Dynamic,1>','Eigen::Matrix<var,Eigen::Dynamic,Eigen::Dynamic>')
  numArgs <- c(1,2,3)
  storageType <- c("scalar","vector","matrix")
  return(data.frame(keyword,cppType,numArgs,storageType))
}

.generated.types <- function(){
  keyword <- c('GENERATED_SCALAR','GENERATED_VECTOR','GENERATED_MATRIX')
  cppType <- c('var','Eigen::Matrix<var,Eigen::Dynamic,1>','Eigen::Matrix<var,Eigen::Dynamic,Eigen::Dynamic>')
  numArgs <- c(1,2,3)
  storageType <- c("scalar","vector","matrix")
  return(data.frame(keyword,cppType,numArgs,storageType))
}


.remove_comments <- function(input){
  output <- input
  # /* - */
  com.on <- FALSE
  for(i in 1:length(output)){
    if(nchar(output[i])>1){
      line <- substring(output[i],first=1,last=1)
      j <- 2
      while(j <= nchar(output[i])){
        c <- substring(output[i],first=j-1,last=j)
        
        if(c=="/*"){
          com.on <- TRUE
          if(j==2){
            line <- ""
          } else {
            line <- substring(line,first=1,last=j-2)
          } 
        }
        if(! com.on) line <- paste0(line,substring(output[i],first=j,last=j))
        if(com.on && c=="*/"){
          com.on <- FALSE
        }
        j <- j+1
      }
      output[i] <- line
    }
  }
  # //
  for(i in 1:length(output)){
    r <- gregexpr("//",output[i])[[1]]
    if(min(r)>0){
      output[i] <- substring(output[i],first=1,last=min(r)-1)
    }
  }
  return(output)
}

.block.in.braces <- function(m.file,initialBrace){
  s <- substr(m.file,initialBrace,initialBrace)
  if(! identical(s,"{")) stop("bad input")
  i <- initialBrace+1
  lastBrace <- -1L
  np <- 1
  
  while(i <= nchar(m.file)){
    s <- substr(m.file,i,i)
    if(identical(s,"{")){
      np <- np+1
    } else if(identical(s,"}")) {
      np <- np-1
    }
    
    if(np==0){
      lastBrace <- i
      break
    }
    
    i <- i+1
  }
  
  if(lastBrace==-1L) {
    warning("did not find braced block")
    blockLen <- 0
  } else {
    blockLen <- lastBrace - initialBrace - 1
  }
  return(list(blockLen = blockLen,first=initialBrace+1,last=lastBrace-1))
}


.resolve.keywords.arguments <- function(arg.string){
  arg.list <- list()
  str.len <- nchar(arg.string)
  if(str.len>0){
    np <- 0
    i <- 1
    arg <- ""
    while(i<=str.len){
      s <- substr(arg.string,i,i)
      if(identical(",",s) && np == 0){
        if(nchar(arg)>0){
          arg.list[[length(arg.list)+1]] <- arg
          arg <- ""
        } else {
          stop("empty argument for keyword")
        } 
      } else {
        arg <- paste0(arg,s)
        if(identical(s,"(")) np<-np+1
        if(identical(s,")")) np<-np-1
      }
      i <- i+1
    }
    if(nchar(arg)>0){
      arg.list[[length(arg.list)+1]] <- arg
    } else {
      stop("empty argument for keyword")
    } 
  } else {
    stop("empty argument for keyword")
  }
  return(arg.list)
}


.locate.keywords <- function(keys,m.file,from,to){
  
  out.keys <- NULL
  out.varNames <- NULL
  out.arg1 <- NULL
  out.arg2 <- NULL
  out.storageType <- NULL
  starts.at <- NULL
  
  for(i in 1:length(keys[,"keyword"])){
    loc <- gregexpr(paste0(keys[i,"keyword"],"\\s*\\(.*?\\)\\s*;"),m.file)[[1]]
    all.loc <- gregexpr(keys[i,"keyword"],m.file)[[1]]
    nfound <- sum(loc>0)
    if(nfound != sum(all.loc>0)) stop("malformed keyword expression found")
    if(nfound>0){
      for(k in 1:nfound){
        if(loc[k]>=from && loc[k]< to){
          # extract complete string, excluding ";" (hence -2)
          completeExpr <- substr(m.file,loc[k],as.integer(loc[k]+attr(loc,"match.length")[k]-2))
          # remove white space and the keyword
          arg.string <- gsub(" ","",gsub(keys[i,"keyword"],"",completeExpr))
          # remove the (outer) parenthesis
          arg.string <- substr(arg.string,2,nchar(arg.string)-1)
          # split the argument list 
          arg.list <- .resolve.keywords.arguments(arg.string)
          
          if(length(arg.list)==keys[i,"numArgs"]){
            starts.at <- c(starts.at,as.integer(loc[k]))
            out.keys <- c(out.keys,keys[i,"keyword"]) 
            out.storageType <- c(out.storageType,keys[i,"storageType"]) 
            out.varNames <- c(out.varNames,arg.list[[1]])
            if(length(arg.list)>=2) {
              out.arg1 <- c(out.arg1,arg.list[[2]])
            } else {
              out.arg1 <- c(out.arg1,"")
            }
            if(length(arg.list)==3) {
              out.arg2 <- c(out.arg2,arg.list[[3]])
            } else {
              out.arg2 <- c(out.arg2,"")
            }
            
          } else {
            stop(paste0("wrong number of arguments in : ",completeExpr))
          }
        } else {
          stop(paste0("found ",keys[i,"keyword"]," outside appropriate block"))
        } # in appropriate block?
      } # loop over found keyse
    } # found given key
  } # loop over keys
  if(! is.null(out.keys)){
    out <- data.frame(out.keys,out.varNames,out.arg1,out.arg2,out.storageType)
    # sort so that variables appear in the same order as in the code
    out <- out[order(starts.at),]
    colnames(out) <- c("keyword","varName","arg1","arg2","storageType")
    rownames(out) <- as.character(1:dim(out)[1])
    return(out)
  } else {
    return(data.frame())
  }
}

.cpp.storage.info <- function(var){
  a1 <- ifelse(identical(var[1,"arg1"],""),"1",var[1,"arg1"])
  a2 <- ifelse(identical(var[1,"arg2"],""),"1",var[1,"arg2"])
  a3 <- switch(var[1,"storageType"],
               scalar="0",
               vector="1",
               matrix="2")
  return(c(a1,a2,a3))
}

.cpp.total.dim <- function(vars){
  ret <- ""
  for(i in 1:dim(vars)[1]){
    cpp.inc <- switch(vars[i,"storageType"],
                      scalar=" 1",
                      vector=paste0(" ",vars[i,"arg1"]),
                      matrix=paste0(" (",vars[i,"arg1"],")*(",vars[i,"arg2"],")"),
                      stop("unknown storageType"))
    ret <- paste0(ret,
                  cpp.inc,
                  ifelse(i==dim(vars)[1],""," +"))
  }
  return(ret)
} 


.is.empty.df<-function(df){return(dim(df)[1]==0)}


