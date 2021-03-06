% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.R
\name{run}
\alias{run}
\title{Run pdphmc model}
\usage{
run(
  model,
  samples = 2000L,
  Tmax = 10000,
  warmupFrac = 0.5,
  chains = 4L,
  seed = 1L,
  data = list(),
  pars = character(0),
  init = list(),
  control = list(absTol = 0.001, relTol = 0.001),
  fixedMiDiag = list(),
  out.files = list(csv.prec = 8),
  clean = TRUE,
  cores = 1L
)
}
\arguments{
\item{model}{S4 class of type `build-output` typically generated by function `build`}

\item{samples}{integer, number of discrete samples}

\item{Tmax}{positive real, (time) length of trajectories}

\item{warmupFrac}{real in [0,1] : proportion of trajectory used for warmup}

\item{chains}{positive integer number of independent realizations}

\item{seed}{integer random number seed}

\item{data}{named list passed to the model}

\item{pars}{character vector names of parameters for which samples are stored, default is all}

\item{init}{named list of initial configurations (otherwise N(0,1) is used)}

\item{control}{named list of different sampler parameters}

\item{fixedMiDiag}{named list giving fixed (i.e. not adaptively learned) parts of inverse mass matrix diagonal}

\item{out.files}{named list giving parameters related to the files returned by the binary}

\item{clean}{logical should output files be deleted after they are loaded?}

\item{cores}{positive integer number of CPU cores to be used (so far requires a unix-type system)}
}
\description{
Run a model generated by `build`
}
\examples{
# EXAMPLE 1 : an independent Gaussian target
# copy model file to current wd
file.copy(system.file("examples/test_case.cpp", package = "pdphmc"),".")
# compile the model
mdl <- pdphmc::build("test_case.cpp")
# run the model
fit <- pdphmc::run(mdl)
# get summary stats of the chains
pdphmc::getMonitor(fit)
# make a trace plot of the x(0) marginal
pdphmc::trace.plot(fit,"x(0)")
# remove files
pdphmc::clean(mdl) # will no longer work as binary is deleted...
# remove model file
file.remove("test_case.cpp")
# EXAMPLE 2 : simple Gaussian mean variance model
# copy model file to current wd
file.copy(system.file("examples/OneSampleGaussian.cpp", package = "pdphmc"),".")
# compile the model
mdl <- pdphmc::build("OneSampleGaussian.cpp")
# simulate some data
y <- rnorm(100)
# run the model
fit <- pdphmc::run(mdl,data=list(y=y))
# get summary stats from the chains
pdphmc::getMonitor(fit)
# make a trace plot of the mu marginal
pdphmc::trace.plot(fit,"mu")
# remove files
pdphmc::clean(mdl) # will no longer work as binary is deleted...
# remove model file
file.remove("OneSampleGaussian.cpp")
}
