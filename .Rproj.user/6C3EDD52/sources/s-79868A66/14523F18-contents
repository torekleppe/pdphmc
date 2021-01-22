
.default.compiler.info<- function(){
  if(identical(.Platform$OS.type,"windows")){
    if(pkgbuild::has_build_tools()){
      rtools.path <- pkgbuild::rtools_path() # note, this is not the path to the compiler binary!!
      compiler.path <- paste0(strsplit(rtools.path,"/")[[1]][1],"\\mingw64\\bin\\g++.exe")
      return(list(compiler=normalizePath(compiler.path,mustWork = TRUE)))
    } else {
      stop("requires a working c++ compiler, get the rtools package")
    }
  } else if(identical(.Platform$OS.type,"unix")) {
    return(list(compiler="g++"))  
  } else {
    stop("Unknown OS.type")
  }
  
}

.compileCpp <- function(bo,
                        compiler.info=.default.compiler.info(),
                        flags="-O3",
                        include=""){
  
  package.includes <- paste0(
    " -I",normalizePath(system.file('include', package = "StanHeaders")),
    " -I",normalizePath(system.file('include', package = "RcppEigen")),
    " -I",normalizePath(system.file('include', package = "BH")),
    " -I",normalizePath(system.file('include', package = "RcppProgress")),
    " -I",normalizePath(system.file('include', package = "pdphmc"))
  )
  
  # compilerCall <- paste0(compiler," ",
  #                        bo@file.name.base,".cpp ",
  #                        " -o ",bo@file.name.base,
  #                        " -std=c++14 ",flags," ",
  #                        package.includes," ",
  #                        include)  
  
  compilerArgs <- paste0(bo@file.name.base,".cpp ",
                         " -o ",bo@file.name.base,
                         " -std=c++14 ",flags," ",
                         package.includes," ",
                         include)  
  
  
  ret <- system2(command=compiler.info$compiler,
                 args=compilerArgs,
                 stdout = TRUE,
                 stderr = TRUE)
  
  eflag <- attr(ret,"status")
  
  if(is.null(eflag)){
    message("compilation exited successfully")
    return(0L)
  } else {
    message("problem with compilation: ")
    cat(paste0(command=compiler.info$compiler," ",compilerArgs))
    cat(ret,file=normalizePath(paste0(bo@file.name.base,"_compiler_out.txt"),mustWork = FALSE))
    print(ret,quote = FALSE)
    return(1L)
  } 
}
