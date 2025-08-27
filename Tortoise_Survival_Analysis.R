## ---------------------------
##
## Tortoise Survival Analysis
##
## This script is used to run models that replicate the analysis of Hromada et. al.,
## "An integrated model improves inferences about survival in the Mojave desert tortoise"
## It requires the script "Tortoise_Survival_Functions.R" and the appropriate 
## datasets (.RDS files) in the working directory to function properly
##
## Author(s): Hromada, S.
##
## Date Created: 2025-08-15
## Date Last Modified: 2025-08-15
## 
## Email: stevehromada@gmail.com
##
## ---------------------------
##
## Notes: All location data has been randomly shifted from the original mark-recapture 
##  datasets, to obscure tortoise locations, though kept in the same geographical
##  relation within a plot to allow the spatial-capture-recapture model to work.
##
## ---------------------------

rm(list=ls())
options(scipen = 6, digits = 4)

## ---------------------------

## load up packages 

library(jagsUI)

## ---------------------------

## load up our functions into memory

source("Tortoise_Survival_Functions.R") 

## ---------------------------

# load up the datasets for each analysis
dfj_t <- readRDS(file = "surv_data_telemetry_public.RDS")
dfj_mrc <- readRDS(file = "surv_data_markrecapture_public.RDS")
dfj_im <- readRDS(file = "surv_data_integrated_public.RDS")

# Load up the JAGS file for each model
dfj_t$filename <- WriteTelFileForJAGS()
dfj_mrc$filename <- WriteMRCFileForJAGS()
dfj_im$filename <- WriteIMFileForJags()

# Run the models with the desired parameters

# Run the telemetry model
mod_tel <-  RunJags(dfj_t, burn = 2000, nits = 30000)

# run the mark-recap model
mod_mrc <- RunJags(dfj_mrc,burn = 25000, nits = 60000)

# run the integrated model
mod_im <- RunJags(dfj_im, burn = 60000, nits = 80000)

# Examine Results
summary(mod_tel)
summary(mod_mrc)
summary(mod_im)

