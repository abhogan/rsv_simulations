# to reinstall the package
#pak::pak("IsaacStopard/RSVsim")
library(tidyverse)
library(RSVsim)
library(foreach)
library(doParallel)
library(patchwork)
library(GGally)
library(glue)

source("data/data.r")
source("R/hosp_function.R")
source("R/set_up_variables.r")
