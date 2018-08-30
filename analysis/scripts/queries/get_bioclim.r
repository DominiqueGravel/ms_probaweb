##################################################################
# Script querying WorldClim to get a raster of climate conditions
# across continental Europe
# Timothée Poisot
# February 15th, 2015
##################################################################

rm(list = ls())

# Load packages
library(dismo)
library(raster)

# Load data
load("analysis/data/salix_int_direct.Rdata")

# Data pre-processing
xl = range(as.numeric(salix_int$lon))
yl = range(as.numeric(salix_int$lat))

# Create the raster
bioclim = crop(getData("worldclim", var="bio", res=5), c(min(xl), max(xl), min(yl), max(yl)))
writeRaster(bioclim, filename="analysis/data/bioclim.grd", overwrite=T)
