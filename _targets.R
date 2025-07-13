# Load libraries for targets
library(targets)
library(tarchetypes)

# optional libraries for better coding
library(languageserver)
library(httpgd)

# Source the functions and sub-pipelines
tar_source(file.path("scripts", "functions", "preprocess.R"))
tar_source(file.path("scripts", "targets_preprocess.R"))

# set options
future::plan("multicore", workers = 6)
options(future.globals.maxSize = 16000 * 1024^2)

# run target pipeline
list(
    targets_preprocess
)

