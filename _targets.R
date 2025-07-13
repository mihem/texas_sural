# Load libraries for targets
library(targets)
library(tarchetypes)

# Source the functions and sub-pipelines
tar_source("scripts")

# set options
future::plan("multicore", workers = 6)
options(future.globals.maxSize = 16000 * 1024^2)

# colors from https://romanhaa.github.io/projects/scrnaseq_workflow/
colors_dutch <- c(
    '#FFC312',
    '#C4E538',
    '#12CBC4',
    '#FDA7DF',
    '#ED4C67',
    '#F79F1F',
    '#A3CB38',
    '#1289A7',
    '#D980FA',
    '#B53471',
    '#EE5A24',
    '#009432',
    '#0652DD',
    '#9980FA',
    '#833471',
    '#EA2027',
    '#006266',
    '#1B1464',
    '#5758BB',
    '#6F1E51'
)

# run target pipeline
list(
    targets_preprocess,
    targets_predict
)
