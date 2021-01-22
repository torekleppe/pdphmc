library(methods)
#' Output of the pdphmc build process
#' 
#' @description Contains path to binary and also specification of the sampler used and other auxiliary information
#' @rdname build-output
#' @slot model.name model name
#' @slot model.cppcode c++ code
#' @slot work.folder folder containing binaries etc
#' @slot build.time build.time
#' @slot build.ID identifier
#' @slot file.name.base absolute path base
#' @slot data.vars names and info on DATA variables in model
#' @slot par.vars names and info on PARAMETER variables in model
#' @slot gen.vars names and info on GENERATED variables in model
#' @slot massMatrix type of mass Matrix
#' @slot lambda type of event rate
#' @slot compile.exit.flag compile exit flag
#' @slot binary.path path to the binary
#' @exportClass build-output 
setClass(
  "build-output",
  slots = c(
    model.name = "character",
    model.cppcode = "character",
    work.folder = "character",
    build.time = "POSIXct",
    build.ID = "integer",
    file.name.base = "character",
    data.vars = "data.frame",
    par.vars = "data.frame",
    gen.vars = "data.frame",
    massMatrix = "character",
    lambda = "character",
    compile.exit.flag = "integer",
    binary.path = "character"
  )
)

#' Class containing pdphmc output and diagnostics information
#' 
#' @rdname pdphmc-output
#' @description pdphmc simulation output class
#' @slot modelname model name
#' @slot model.obj copy of the build-output object
#' @slot niter not currently in use
#' @slot chains number of trajectories
#' @slot warmup number of samples in warmup mode
#' @slot samples number of samples 
#' @slot pointSamples array containing (point) samples (sample,chain,parameter)
#' @slot intSamples array containing integrated samples (sample,chain,parameter)
#' @slot eventSamples not currently in use
#' @slot diagonstics diagnostics output for each integration leg
#' @slot monitor not in use
#' @slot intMonitor not in use
#' @slot par.vec.names names of parameter vector
#' @slot CPUtime CPUtime
#' @slot fixedMi inverse mass diagonal fixed
#' @slot lastMiDiag last inverse mass diagonal
#' @exportClass pdphmc-output
setClass(
  "pdphmc-output",
  slots = c(
    modelname = "character",
    model.obj = "build-output",
    niter = "integer",
    chains = "integer",
    warmup = "integer",
    samples = "integer",
    pointSamples = "array",
    intSamples = "array",
    eventSamples = "list",
    diagnostics = "list",
    monitor = "list",
    intMonitor = "list",
    par.vec.names = "character",
    CPUtime = "list",
    fixedMi = "numeric",
    lastMiDiag = "list"
  )
)
