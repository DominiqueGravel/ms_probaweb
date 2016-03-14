##################################################################
# Script querying WorldClim to get the climate data for each location
# of the Salix dataset
# Timothée Poisot
# February 15th, 2015
##################################################################

rm(list = ls())
setwd("/Users/DGravel/Documents/Manuscripts/Inprep/ms_probaweb/analysis")

# Load packages
library(dismo)
library(raster)

# Load data
load("analysis/data/salix_int.Rdata")

# Data proprocessing
salix_int$lat = as.numeric(salix_int$lat)
salix_int$lon = as.numeric(salix_int$lon)

xl = range(as.numeric(salix_int$lon))
yl = range(as.numeric(salix_int$lat))

# Get the bioclim data
bclim = brick("analysis/data/bioclim.grd")

# Extract for coordinates
salix_climate = extract(bclim, salix_int[,c('lon', 'lat')])
salix_climate = data.frame(salix_int, salix_climate)

save(salix_climate, file="analysis/data/salix_climate.Rdata")
