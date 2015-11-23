##################################################################
# Script for Figure X,
# Mapping across Europe the interaction probability between pairs of species
# Dominique Gravel
# October 29th, 2015
##################################################################

rm(list=ls())

	library(raster)
library(sp)
library(rgdal)

load("analysis/data/DF_split.Rdata")
load("analysis/data/expand_data.Rdata")

IDi = data$pairs.IDi
IDj = data$pairs.IDj
Si = length(unique(IDi))
Sj = length(unique(IDj))
nobs = 374
np  = length(DF_split)

# Load the stuff to make the raster object
wrld <- readOGR("analysis/data/map", layer="level1")
europe <- wrld[1,]
exteur <- extent(europe)

# Extract the climatic variables
clim = crop(getData("worldclim", var="bio", res=10), exteur)
T = getValues(clim$bio1)
PP = getValues(clim$bio12)

############################################################
# Extract the local networks over Europe
############################################################

ncell = length(T)
map_E = data.frame(T = T, PP = PP, T2 = T^2, PP2 = PP^2)
#pair = 14038
#pair = 14456
#pair = 14458
#pair = 14476
#pair = 14618
#pair = 15626
#pair = 21395
#pair = 21415
pair = 21418

data = DF_split[[pair]]
Xi = data$Xi
Xj = data$Xj
Xij = data$Xij
Lij = data$Lij
E = data$E	

i = as.numeric(as.vector(data$IDi[1]))
j = as.numeric(as.vector(data$IDj[1]))

# Compute the models
lm_i = glm(Xi ~ T+PP+T2+PP2, family = "binomial", data = E)
lm_j = glm(Xj ~ T+PP+T2+PP2, family = "binomial", data = E)
subE = subset(E,Xij == 1)
subLij = subset(Lij, Xij == 1)
lm_L = glm(subLij ~ T+PP+T2+PP2, family = "binomial", data = subE)

summary(lm_i)
summary(lm_j)
summary(lm_L)

# Compute expected probabilities over the gradient
PXi_E = predict(lm_i, type="response", newdata=map_E)
PXj_E = predict(lm_j, type="response", newdata=map_E)
PLij_E = predict(lm_L, type="response", newdata=map_E)
PXij_E = PXi_E*PXj_E
PLijXij_E = PXij_E*PLij_E

############################################################
# Plot network properties
############################################################

quartz(height = 3.5, width = 10)
par(mfrow = c(1,3),mar = c(2,2,3,2))

# Map parameters
library("colorRamps")
library("RColorBrewer")
pal <-colorRampPalette(rev(brewer.pal(11,"RdYlBu"))) # Initialized Colors ramp palette, see here: http://colorbrewer2.org/

# MAP 1: CO-OCCURRENCE PORBABILITY
# Make the map
par(xaxs="i", yaxs="i")
plot(exteur,type="n", axes=FALSE, ann=FALSE)
rect(xmin(exteur),ymin(exteur),xmax(exteur),ymax(exteur), col="#9CDDF1")
rast <- raster(exteur, vals=(PXij_E), nrow = 276, ncol = 544)
image(rast, add=TRUE, col=pal(100))
plot(europe, border="grey25", lwd=1.2, add=TRUE)
plot(wrld[2:3,], add=TRUE, col="white")
mtext(text=expression(P(X[i],X[j])),side=3,line=0.5,adj=-0.1,cex=1.25)


# MAP 2: INTERACTION PROBABILITY
# Make the map
par(xaxs="i", yaxs="i")
plot(exteur,type="n", axes=FALSE, ann=FALSE)
rect(xmin(exteur),ymin(exteur),xmax(exteur),ymax(exteur), col="#9CDDF1")
rast <- raster(exteur, vals=PLij_E, nrow = 276, ncol = 544)
image(rast, add=TRUE, col=pal(100))
plot(europe, border="grey25", lwd=1.2, add=TRUE)
plot(wrld[2:3,], add=TRUE, col="white")
mtext(text=expression(P(L[ij])),side=3,line=0.5,adj=-0.1,cex=1.25)


# MAP 3: NET INTERACTION PROBABILITY
# Make the map
par(xaxs="i", yaxs="i")
plot(exteur,type="n", axes=FALSE, ann=FALSE)
rect(xmin(exteur),ymin(exteur),xmax(exteur),ymax(exteur), col="#9CDDF1")
rast <- raster(exteur, vals=(PLijXij_E), nrow = 276, ncol = 544)
image(rast, add=TRUE, col=pal(100))
plot(europe, border="grey25", lwd=1.2, add=TRUE)
plot(wrld[2:3,], add=TRUE, col="white")
mtext(text=expression(P(L[ij],X[i],X[j])),side=3,line=0.5,adj=-0.1,cex=1.25)

dev.copy2pdf(file = "ms/figures/map_pair.pdf")




