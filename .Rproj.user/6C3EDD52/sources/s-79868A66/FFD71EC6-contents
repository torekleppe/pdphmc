#' Build and compile pdphmc model
#' 
#' @description 
#' `build` parses a model file and compiles the program used by `run`
#'
#'
#' @param file path to file containing model
#' @param model.name character model name
#' @param work.folder character where to put binary and other temporary files
#' @param massMatrix character type of mass matrix used
#' @param lambda character type of event rate used
#' @param target.copies integer, number of copies of the target distribution used in extended target, currently only 1 and 2 implemented
#' @param include.flags character additional flags put at the end of compiler call.
#' @param verbose integer controls how much information is printed during parsing.
#' @param cpp.only logical if true, then only build c++ file and do not compile it.
#' @param clean logical if true, remove temporary files not used by `run`
#' @param debug logical compile with debug macro
#' @param compiler.info Information on which compiler to use
#' @export
build <- function(file,
                  model.name=strsplit(file,split=".",fixed = TRUE)[[1]][1],
                  work.folder=paste0(getwd(),"/",model.name,"_files/"),
                  massMatrix=c("diagMassISG","diagMassISG_P","diagMassVARI","identityMass","diagMassFixed"),
                  lambda=c("lambda_constant","lambda_arclength"),
                  target.copies = 1L,
                  include.flags = "",
                  verbose = 1L,
                  cpp.only = FALSE,
                  clean = TRUE,
                  debug = FALSE,
                  compiler.info=.default.compiler.info()
                  ){
  
  
  out <- methods::new("build-output") 
  out@model.name <- model.name
  
  
  
  # make folder for output files
  if(! file.exists(normalizePath(work.folder,mustWork = FALSE))){
    dir.create(normalizePath(work.folder,mustWork = FALSE))
  }
  out@work.folder <- work.folder
  
  # get time of build
  build.time <- Sys.time()
  out@build.time <- build.time
  build.ID <- as.integer(format(build.time,"%m%d%H%M%S"))
  out@build.ID <- build.ID
  
  # files in work folder will start with this name
  file.name.base <- sub("//","/",paste0(work.folder,"/",model.name,"_",build.ID))
  out@file.name.base <- file.name.base
  
  
  # read input file
  model.file <- readLines(file) 
  
  #initial cleaning of the model file
  model.file <- paste0(.remove_comments(model.file),collapse = "\n")
  
  # read template
  out.file <- paste0(readLines(normalizePath(paste0(system.file('mainTemplate', package = "pdphmc"),
                                      "/main_body_template.cpp"))),collapse = "\n")
  
  
  ################################################################################
  # build flags flags
  ################################################################################
  build.flags <- ""
  if(debug){
    build.flags <- "#define __DEBUG__ \n"
  }
  out.file <- sub(pattern="//{[[BUILD_FLAGS]]}",
                  replacement = build.flags,
                  x=out.file,fixed=TRUE)
  
  
  ################################################################################
  # Fill in sampler templates
  ################################################################################
  
  if(target.copies==1L){
    ADwrap = "ADtarget"
  } else if(target.copies==2L) {
    ADwrap = "ADtarget_cp2"
  } else {
    stop("bad value of target.copies, currently only 1 and 2 implemented")
  }
  
  
  message(paste0("ADwrapper: ", ADwrap))
  out.file <- sub(pattern="{[[AD_WRAPPER_CLASS]]}",
                  replacement = ADwrap,
                  x=out.file,fixed=TRUE)
  
  massMatrix <- match.arg(massMatrix)
  message(paste0("mass matrix: ", massMatrix))
  out.file <- sub(pattern="{[[MASS_CLASS]]}",
                  replacement = massMatrix,
                  x=out.file,fixed=TRUE)
  
  lambda <- match.arg(lambda)
  message(paste0("lambda: ",lambda))
  out.file <- sub(pattern="{[[LAMBDA_CLASS]]}",
                  replacement = lambda,
                  x=out.file,fixed=TRUE)
  
  
  ################################################################################
  # Locate the blocks
  ################################################################################
  block.start <- rep(-1,length(.blocks))
  block.end <- rep(-1,length(.blocks))
  
  for(i in 1:length(.blocks)){
    initPos <- gregexpr(pattern=paste0(.blocks[i],"\\s*\\{"),
                        text=model.file,fixed=FALSE,useBytes = TRUE)[[1]]
    
    if(length(initPos)>1){
      stop(paste0("more than one instance of block ",.blocks[i]))
    } else if(initPos>=1){
      initBrace <- as.integer(initPos+attr(initPos,"match.length")-1)
      block.loc <- .block.in.braces(m.file = model.file,initBrace)
      if(block.loc$blockLen>0){
        block.start[i] <- block.loc$first
        block.end[i] <- block.loc$last
        #print(.blocks[i])
        #cat(substr(model.file,block.start[i],block.end[i]))
      }
    }
  }
  
  if(block.start[4]==-1 || block.start[5]==-1){
    stop("model file is missing a VARIABLES_BLOCK or a MODEL_BLOCK block")
  }
  
  
  ###########################################################################
  # Fill inn model file blocks at appropriate places
  ###########################################################################
  # INCLUDE BLOCK
  out.file <- sub(pattern="//{[[INCLUDE_BLOC]]}",
                  replacement = substr(model.file, block.start[1], block.end[1]),
                  x=out.file ,fixed=TRUE)
  
  
  # DATA BLOCK
  out.file <- sub(pattern="//{[[DATA_BLOC]]}",
                  replacement = substr(model.file, block.start[2], block.end[2]),
                  x=out.file ,fixed=TRUE)
  
  # SETUP BLOCK
  out.file <- sub(pattern="//{[[SETUP_BLOCK]]}",
                  replacement = substr(model.file, block.start[3], block.end[3]),
                  x=out.file ,fixed=TRUE)
  
  # VARIABLES BLOCK
  out.file <- sub(pattern="//{[[VARIABLES_BLOCK]]}",
                  replacement = substr(model.file, block.start[4], block.end[4]),
                  x=out.file ,fixed=TRUE)
  
  
  # MODEL BLOCK
  out.file <- sub(pattern="//{[[MODEL_BLOCK]]}",
                  replacement = substr(model.file, block.start[5], block.end[5]),
                  x=out.file ,fixed=TRUE)
  
  ###########################################################################
  # Locate DATA variables
  ###########################################################################
  if(block.start[2]>0){
    data.vars <- .locate.keywords(.data.types(),model.file,block.start[2],block.end[2])
  } # model has a DATA_BLOCK
  if(verbose>=1){
    if(.is.empty.df(data.vars)){
      message("no DATA keywords found")
    } else {
      message("found the following DATA keywords:")
      print(data.vars[,c("keyword","varName")])
    }
  }
  out@data.vars <- data.vars
  
  
  ###########################################################################
  # Locate PARAMETER and GENERATED variables
  ###########################################################################
  
  par.vars <- .locate.keywords(.parameter.types(),model.file,block.start[4],block.end[4])
  if(.is.empty.df(par.vars)) stop("model does not contain any parameters !!!")
  if(verbose>=1){
    message("found the following PARAMETER keywords:")
    print(par.vars)
  }
  out@par.vars <- par.vars
  
  gen.vars <- .locate.keywords(.generated.types(),model.file,block.start[4],block.end[4])
  if(verbose>=1){
    if(.is.empty.df(gen.vars)){
      message("no GENERATED keywords found")
    } else {
      message("found the following GENERATED keywords:")
      print(gen.vars)
    }
  }
  out@gen.vars <- gen.vars
  
  ###########################################################################
  # make code dimensions of problem 
  ###########################################################################
  
  cpp.par.total.dim <- .cpp.total.dim(par.vars)
  out.file <- sub(pattern="//{[[PARAMETER_DIM]]}",
                  replacement = cpp.par.total.dim,
                  x=out.file ,fixed=TRUE)
  
  cpp.gen.total.dim <- ifelse(.is.empty.df(gen.vars),"0",.cpp.total.dim(gen.vars))
  out.file <- sub(pattern="//{[[GENERATED_DIM]]}",
                  replacement = cpp.gen.total.dim,
                  x=out.file ,fixed=TRUE)
  
  ###########################################################################
  # make code names and dimensions of problem for class __VARIABLES_dimension
  ###########################################################################
  
  cpp.par.names.push = ""
  cpp.par.info <- paste0("par_info_.resize(",dim(par.vars)[1],",5);\n")
  
  for(i in 1:dim(par.vars)[1]){
    cpp.par.names.push <- paste0(cpp.par.names.push,
                                 "par_names_.push_back(\"",
                                 par.vars[i,"varName"],"\");\n")
    tmp <- .cpp.storage.info(par.vars[i,])
    cpp.par.info <- paste0(cpp.par.info,
                           "par_info_(",i-1,",0) = ",tmp[1],";\n",
                           "par_info_(",i-1,",1) = ",tmp[2],";\n",
                           "par_info_(",i-1,",2) = ",tmp[3],";\n"
    )
  } 
  out.file <- sub(pattern="//{[[PARAMETER_NAMES_PUSH]]}",
                  replacement = cpp.par.names.push,
                  x=out.file, fixed=TRUE)
  
  out.file <- sub(pattern="//{[[PARAMETER_INFO]]}",
                  replacement = cpp.par.info,
                  x=out.file, fixed=TRUE)
  
  cpp.gen.names.push <- ""
  cpp.gen.info <- ""
  if(! .is.empty.df(gen.vars)){
    cpp.gen.info <- paste0("gen_info_.resize(",dim(gen.vars)[1],",5);\n")
    for(i in 1:dim(gen.vars)[1]){
      cpp.gen.names.push <- paste0(cpp.gen.names.push,
                                   "gen_names_.push_back(\"",
                                   gen.vars[i,"varName"],"\");\n")
      tmp <- .cpp.storage.info(gen.vars[i,])
      cpp.gen.info <- paste0(cpp.gen.info,
                             "gen_info_(",i-1,",0) = ",tmp[1],";\n",
                             "gen_info_(",i-1,",1) = ",tmp[2],";\n",
                             "gen_info_(",i-1,",2) = ",tmp[3],";\n"
      )
    } 
  }
  out.file <- sub(pattern="//{[[GENERATED_NAMES_PUSH]]}",
                  replacement = cpp.gen.names.push,
                  x=out.file, fixed=TRUE)
  
  out.file <- sub(pattern="//{[[GENERATED_INFO]]}",
                  replacement = cpp.gen.info,
                  x=out.file, fixed=TRUE)
  
  ###########################################################################
  # fill in file name of JSON file
  ###########################################################################
  
  
  json.path.in.cpp <- paste0("\""
                             ,normalizePath(paste0(file.name.base,".json"),mustWork = FALSE),
                             "\"")
  if(identical(.Platform$OS.type,"windows")){
    json.path.in.cpp <- gsub(pattern = "\\",replacement = "\\\\",x=json.path.in.cpp,fixed=TRUE)
  }
  out.file <- sub(pattern="//{[[JSON_FILE_NAME]]}",
                  replacement = json.path.in.cpp,
                  x=out.file ,fixed=TRUE)
  
  ###########################################################################
  # make code getting data from JSON file
  ###########################################################################
  
  
  cpp.load.data <- ""
  if(! .is.empty.df(data.vars)){
    for(i in 1:dim(data.vars)[1]){
      vn <- data.vars[i,"varName"]
      cpp.load.data <- paste0(cpp.load.data,
                              "if(! __jw.getNumeric(\"data\",\"",vn,"\",",vn,")){ \n",
                              "  std::cout << \"WARNING: missing / wrong type of data provided",
                              "for DATA variable ",vn," !!! \" << std::endl;\n",
                              "}\n")
    }
  }
  out.file <- sub(pattern="//{[[READDATA]]}",
                  replacement = cpp.load.data,
                  x=out.file ,fixed=TRUE)
  
  
  ###########################################################################
  # Make code copying to PARAMETER variables from the collection 
  # of parameters __par
  ###########################################################################
  
  cpp.cp.from.par <- "int __par_counter = 0;\n"
  dim.tmp <- ""
  for(i in 1:dim(par.vars)[1]){
    dim.inc <- switch(par.vars[i,"storageType"],
                      scalar="1",
                      vector=paste0(par.vars[i,"varName"],".size()"),
                      matrix=paste0(par.vars[i,"varName"],".rows()*",
                                    par.vars[i,"varName"],".cols()"),
                      stop("unknown storageType")
    )
    cpp.inc <- switch(par.vars[i,"storageType"],
                      scalar=paste0(par.vars[i,"varName"],
                                    " = __par.coeff(__par_counter); \n",
                                    "__par_counter++; \n"),
                      vector=paste0(par.vars[i,"varName"],
                                    " = __par.segment(__par_counter,",
                                    par.vars[i,"varName"],".size()); \n",
                                    "__par_counter += ",
                                    par.vars[i,"varName"],".size(); \n"),
                      matrix=paste0("for(int __i = 0; __i < ",
                                    par.vars[i,"varName"],
                                    ".cols() ; __i++){\n",
                                    "  ",par.vars[i,"varName"],
                                    ".col(__i) = __par.segment(__par_counter,",
                                    par.vars[i,"varName"],".rows());\n",
                                    "  ","__par_counter += ",
                                    par.vars[i,"varName"],".rows();\n",
                                    "}\n"
                      ),
                      stop("unknown storageType"))
    dim.tmp <- paste0(dim.tmp,dim.inc,ifelse(i==dim(par.vars)[1],"","+"))
    cpp.cp.from.par <- paste0(cpp.cp.from.par,cpp.inc,"\n")
  }
  
  cpp.cp.from.par <- paste0("if(",dim.tmp," != __par.size()){\n",
                            "  std::cout << \" parameters have been resized, please revise model, exiting \" << std::endl;\n",
                            "  return(-1);\n",
                            "}\n",
                            cpp.cp.from.par)
  
  
  
  out.file <- sub(pattern="//{[[CP_FROM_PAR]]}",
                  replacement = cpp.cp.from.par,
                  x=out.file ,fixed=TRUE)
  
  
  ###########################################################################
  # Make code copying to GENERATED variables to the collection of 
  # generated variables __gen
  ###########################################################################
  cpp.cp.to.gen <- "\nint __gen_counter = 0;\n"
  
  if(! .is.empty.df(gen.vars)){
    dim.tmp <- ""
    for(i in 1:dim(gen.vars)[1]){
      dim.inc <- switch(gen.vars[i,"storageType"],
                        scalar="1",
                        vector=paste0(gen.vars[i,"varName"],".size()"),
                        matrix=paste0(gen.vars[i,"varName"],".rows()*",
                                      gen.vars[i,"varName"],".cols()"),
                        stop("unknown storageType")
      )
      
      cpp.inc <- switch(gen.vars[i,"storageType"],
                        scalar=paste0("__gen.coeffRef(__gen_counter) = ",
                                      "doubleValue(",gen.vars[i,"varName"],");\n",
                                      "__gen_counter++;\n"
                        ),
                        vector=paste0("__gen.segment(__gen_counter,",
                                      gen.vars[i,"varName"],".size()) = ",
                                      "doubleValue(",gen.vars[i,"varName"],");\n",
                                      "__gen_counter += ",
                                      gen.vars[i,"varName"],".size();\n"
                        ),
                        matrix=paste0("for(int __i=0; __i < ",gen.vars[i,"varName"],
                                      ".cols(); __i++){\n",
                                      "  ","__gen.segment(__gen_counter,",
                                      gen.vars[i,"varName"],".rows()) = ",
                                      "doubleValue(",gen.vars[i,"varName"],".col(__i));\n",
                                      "  ","__gen_counter += ",
                                      gen.vars[i,"varName"],".rows();\n",
                                      "}\n"
                        ),
                        stop("unknown storageType")
      )
      dim.tmp <- paste0(dim.tmp,dim.inc,ifelse(i==dim(gen.vars)[1],"","+"))
      cpp.cp.to.gen <- paste0(cpp.cp.to.gen,cpp.inc,"\n")
    }
    cpp.cp.to.gen <- paste0("if(",dim.tmp," != __gen.size()){\n",
                            "  std::cout << \" generated quantities have been resized, please revise model, exiting \" << std::endl;\n",
                            "  return(-1);\n",
                            "}\n",
                            cpp.cp.to.gen)
  }
  
  out.file <- sub(pattern="//{[[CP_TO_GEN]]}",
                  replacement = cpp.cp.to.gen,
                  x=out.file ,fixed=TRUE)
  
  ###########################################################################
  # Done building the c++ file, now write it to file
  ###########################################################################
  
  out@model.cppcode <- out.file
  cpp.file.path <- normalizePath(paste0(file.name.base,".cpp"),mustWork = FALSE)
  cat(out.file,file=cpp.file.path)
  
  ###########################################################################
  # compile c++ file
  ###########################################################################
  if(! cpp.only){
    
    compile.exit.flag <- .compileCpp(out,compiler=compiler.info,
                                       include=include.flags)
    
    if(clean) file.remove(cpp.file.path)
    
    if(identical(.Platform$OS.type,"windows")){
      out@binary.path <- normalizePath(paste0(file.name.base,".exe"))
    } else {
      out@binary.path <- file.name.base
    }
  } else {
    compile.exit.flag <- NULL
  }
  
  out@compile.exit.flag <- compile.exit.flag
  return(out)
}


##########################################################################
