#
# methods for pdphmc-output objects
#
#' Obtain summary statistics from pdphmc output using rstan::monitor
#' 
#' @description Obtain summary statistics from pdphmc output using rstan::monitor
#' @param object S4 class of type `pdphmc-output` typically generated by function `run` 
#' @param print logical passed to rstan::monitor
#' @export
getMonitor<-function(object,print=TRUE){
  if(!identical(methods::is(object),"pdphmc-output")) stop("bad input")
  object@monitor <- list(rstan::monitor(sims=object@pointSamples[2:(object@samples+1),,], warmup = object@warmup, print=print))
  return(object@monitor)
}


#' Obtain summary statistics from pdphmc integrated output using rstan::monitor
#' 
#' @description Obtain summary statistics from integrated pdphmc output using rstan::monitor
#' @param object S4 class of type `pdphmc-output` typically generated by function `run` 
#' @param print logical passed to rstan::monitor
#' @export
getIntMonitor<-function(object,print=TRUE){
  if(!identical(methods::is(object),"pdphmc-output")) stop("bad input")
  if(print) {
    message("Warning: showing monitor for *integrated samples*")
    message("only mean and se_mean reflect the original target.")
  }
  object@intMonitor <- list(rstan::monitor(sims=object@intSamples[2:(object@samples+1),,], warmup = object@warmup, print=print))
  return(object@intMonitor)
}


#' trace plot from pdphmc output
#' 
#' @description trace plot from pdphmc output
#' @param object S4 class of type `pdphmc-output` typically generated by function `run` 
#' @param which.par character vector of parameter/generated names
#' @export
trace.plot<-function(object,which.par=NULL){
  if(!identical(methods::is(object),"pdphmc-output")) stop("bad input")
  if(is.null(which.par)){
    nms <- dimnames(object@pointSamples)[[3]]
  } else {
    nms <- which.par
  }
  nn <- length(nms)
  graphics::par(mfrow=c(nn,1))
  xlab <- ""
  for(i in 1:nn){
    if(i==nn) xlab <- "iteration #"
    stats::ts.plot(as.vector(object@pointSamples[2:(object@samples+1),,nms[i]]),ylab=nms[i],xlab=xlab )
    if(object@chains>1){
      for(j in 0:object@chains){
        graphics::lines(j*object@samples*c(1,1),c(-1.0e300,1.0e300),col="red")
      }
    }
  }
}


#' trace plot of integrated samples from pdphmc output
#' 
#' @description trace plot of integrated samples from pdphmc output
#' @param object S4 class of type `pdphmc-output` typically generated by function `run` 
#' @param which.par character vector of parameter/generated names
#' @export
trace.plot.int<-function(object,which.par=NULL){
  if(!identical(methods::is(object),"pdphmc-output")) stop("bad input")
  if(is.null(which.par)){
    nms <- dimnames(object@intSamples)[[3]]
  } else {
    nms <- which.par
  }
  nn <- length(nms)
  graphics::par(mfrow=c(nn,1))
  xlab <- ""
  main <- "note; integrated samples"
  for(i in 1:nn){
    if(i==nn) xlab <- "iteration #"
    stats::ts.plot(as.vector(object@intSamples[2:(object@samples+1),,nms[i]]),ylab=nms[i],xlab=xlab,main=main )
    main <- ""
    if(object@chains>1){
      for(j in 0:object@chains){
        graphics::lines(j*object@samples*c(1,1),c(-1.0e300,1.0e300),col="red")
      }
    }
  }
}

#' ACF plot from pdphmc output
#' 
#' @description autocorrelation plot from pdphmc output
#' @param object S4 class of type `pdphmc-output` typically generated by function `run` 
#' @param which.par character vector of parameter/generated names
#' @export
acf.plot<-function(object,which.par=NULL){
  if(!identical(methods::is(object),"pdphmc-output")) stop("bad input")
  if(is.null(which.par)){
    nms <- dimnames(object@pointSamples)[[3]]
  } else {
    nms <- which.par
  }
  nn <- length(nms)
  graphics::par(mfrow=c(nn,object@chains))
  
  for(i in 1:nn){
    ylab <- nms[i]
    
    for(j in 1:object@chains){
      if(i==1){
        main <- paste0("chain # ",j)
      } else {
        main <- ""
      }
      stats::acf(object@pointSamples[(object@warmup+1):(object@samples+1),j,nms[i]],ylab=ylab,
          main=main)
      ylab <- ""
    }
  }
  
}
#' ACF plot of integrated samples from pdphmc output
#' 
#' @description autocorrelation plot of integrated samples from pdphmc output
#' @param object S4 class of type `pdphmc-output` typically generated by function `run` 
#' @param which.par character vector of parameter/generated names
#' @export
acf.plot.int<-function(object,which.par=NULL){
  if(!identical(methods::is(object),"pdphmc-output")) stop("bad input")
  if(is.null(which.par)){
    nms <- dimnames(object@intSamples)[[3]]
  } else {
    nms <- which.par
  }
  nn <- length(nms)
  graphics::par(mfrow=c(nn,object@chains))
  
  for(i in 1:nn){
    ylab <- nms[i]
    
    for(j in 1:object@chains){
      if(i==1){
        main <- paste0("chain # ",j)
      } else {
        main <- ""
      }
      stats::acf(object@intSamples[(object@warmup+1):(object@samples+1),j,nms[i]],ylab=ylab,
          main=main)
      ylab <- ""
    }
  }
  
}
#' Histogram plot from pdphmc output
#' 
#' @description Histogram plot from pdphmc output
#' @param object S4 class of type `pdphmc-output` typically generated by function `run` 
#' @param which.par character vector of parameter/generated names
#' @export
histogram.plot<-function(object,which.par=NULL){
  if(!identical(methods::is(object),"pdphmc-output")) stop("bad input")
  if(is.null(which.par)){
    nms <- dimnames(object@pointSamples)[[3]]
  } else {
    nms <- which.par
  }
  nn <- length(nms)
  graphics::par(mfrow=c(nn,object@chains))
  
  for(i in 1:nn){
    ylab <- nms[i]
    
    for(j in 1:object@chains){
      if(i==1){
        main <- paste0("chain # ",j)
      } else {
        main <- ""
      }
      graphics::hist(object@pointSamples[(object@warmup+1):(object@samples+1),j,nms[i]],ylab=ylab,
           main=main,probability = TRUE,xlab="")
      ylab <- ""
    }
  }
  
}

#' extract samples from pdphmc output
#' 
#' @description extract samples of specific parameters from pdphmc output
#' @param object S4 class of type `pdphmc-output` typically generated by function `run` 
#' @param which.par character vector of parameter/generated names
#' @param include.warmup logical 
#' @export
getSample<-function(object,which.par=NULL,include.warmup=FALSE){
  if(!identical(methods::is(object),"pdphmc-output")) stop("bad input")
  if(is.null(which.par)){
    nms <- dimnames(object@pointSamples)[[3]]
  } else {
    nms <- which.par
  }
  nn <- length(nms)
  dd <- object@samples
  
  if(!include.warmup){dd <- dd - object@warmup}
  
  ret <- matrix(0.0,nrow = dd*object@chains, ncol=nn)
  
  for(i in 1:nn){
    for(j in 1:object@chains){
      if(include.warmup){
        ret[((j-1)*dd+1):(j*dd),i] = object@pointSamples[2:(object@samples+1),j,nms[i]]  
      } else {
        ret[((j-1)*dd+1):(j*dd),i] = object@pointSamples[(object@warmup+2):(object@samples+1),j,nms[i]]
      }
    }
  }
  return(ret)
  
}

#' Get CPU times from pdphmc output
#' 
#' @description CPU times from pdphmc output
#' @param object S4 class of type `pdphmc-output` typically generated by function `run` 
#' @export
get_CPU_time<-function(object){
  if(!identical(methods::is(object),"pdphmc-output")) stop("bad input")
  l <- object@CPUtime
  ret <- matrix(0.0,nrow=length(l),ncol=2)
  for(i in 1:length(l)){
    t <- l[[i]]
    ret[i,] <- t[1:2]
  }
  colnames(ret) <- c("warmup","sampling")
  return(ret)
}
