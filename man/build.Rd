% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.R
\name{build}
\alias{build}
\title{Build and compile pdphmc model}
\usage{
build(
  file,
  model.name = strsplit(file, split = ".", fixed = TRUE)[[1]][1],
  work.folder = paste0(getwd(), "/", model.name, "_files/"),
  massMatrix = c("diagMassISG", "diagMassISG_P", "diagMassVARI", "identityMass",
    "diagMassFixed"),
  lambda = c("lambda_constant", "lambda_arclength"),
  target.copies = 1L,
  include.flags = "",
  verbose = 1L,
  cpp.only = FALSE,
  clean = TRUE,
  debug = FALSE,
  compiler.info = .default.compiler.info()
)
}
\arguments{
\item{file}{path to file containing model}

\item{model.name}{character model name}

\item{work.folder}{character where to put binary and other temporary files}

\item{massMatrix}{character type of mass matrix used}

\item{lambda}{character type of event rate used}

\item{target.copies}{integer, number of copies of the target distribution used in extended target, currently only 1 and 2 implemented}

\item{include.flags}{character additional flags put at the end of compiler call.}

\item{verbose}{integer controls how much information is printed during parsing.}

\item{cpp.only}{logical if true, then only build c++ file and do not compile it.}

\item{clean}{logical if true, remove temporary files not used by `run`}

\item{debug}{logical compile with debug macro}

\item{compiler.info}{Information on which compiler to use}
}
\description{
`build` parses a model file and compiles the program used by `run`
}
