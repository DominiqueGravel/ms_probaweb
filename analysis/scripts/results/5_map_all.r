##################################################################
# Script for Figure 6
# Mapping network properties across Europe
# Dominique Gravel
# October 29th, 2015
##################################################################

rm(list=ls())
library(raster)
library(sp)
library(rgdal)

load("data/DF_split.Rdata")
load("data/expand_data.Rdata")
load("data/pairs.Rdata")

IDi = data$pairs.IDi
IDj = data$pairs.IDj
Si = length(unique(IDi))
Sj = length(unique(IDj))
nobs = 374
np  = length(DF_split)

# Load the stuff to make the raster object
wrld <- readOGR("data/map", layer="level1")
europe <- wrld[1,]
exteur <- extent(europe)

# Extract the climatic variables
clim = crop(getData("worldclim", var="bio", res=10), exteur)
mapT = getValues(clim$bio1)
mapPP = getValues(clim$bio12)

############################################################
# Extract the local networks over Europe
############################################################

ncell = length(T)
map_E = data.frame(T = mapT, PP = mapPP, T2 = mapT^2, PP2 = mapPP^2)

# Loop around all species pairs
expS = numeric(ncell)
expL = numeric(ncell)


#############
# Compute the expected local species richness
# Stack 
Xi = data[,7]
Xj = data[,8]
site = data[,1]
T = data[,12]
PP = data[,23]

data_i = cbind(site,IDi,Xi,T,PP)
data_j = cbind(site,IDj,Xj,T,PP)

all_pres = rbind(data_i,data_j)
site.ID = paste(all_pres[,1],all_pres[,2],sep = "-")

# Make sure there is only one site-species observation
pres = unique(all_pres)

# Loop across all species
S = length(unique(pres[,2]))
IDs = unique(pres[,2])

for(i in 1:S) {
	subdata = subset(pres,pres[,2]==IDs[i])
	Xi = subdata[,3]
	T = subdata[,4]
	PP = subdata[,5]
	T2 = T*T
	PP2 = PP*PP
	lmXi = glm(Xi ~ T+PP+T2+PP2, family = "binomial")
	predXi = predict(lmXi, type = "response", newdata = map_E)
	expS = expS + predXi
}

#############
# Compute the expected link density
for(pair in 1:np) {

	data = DF_split[[pair]]
	Xi = data$Xi
	Xj = data$Xj
	Xij = data$Xij
	Lij = data$Lij
	E = data$E	
	i = as.numeric(as.vector(data$IDi[1]))
	j = as.numeric(as.vector(data$IDj[1]))

	# Compute the models
	if(sum(Xij!=0)) {
		lm_Xij = glm(Xij ~ T+PP+T2+PP2, family = "binomial", data = E)
		expLij = sum(Lij)/sum(Xij)
	}

		else {
			expXij = 0
			expLij = 0
	}

	# Compute expected probabilities over the gradient
	if(sum(Xij!=0)) {
		expXij = predict(lm_Xij, type="response", newdata=map_E)
		# Compute species richness and link density
		expL = expL + expLij*expXij
	}

	cat(pair,'\n')
}

write.table(cbind(expS,expL),"figures/map_data.txt")

############################################################
# Plot network properties
############################################################

dev.new(height = 3.5, width = 10)
par(mfrow = c(1,3),mar = c(2,2,3,2))

# Map parameters
library(viridis)
pal <- viridis

# MAP 1: CO-OCCURRENCE PORBABILITY
# Make the map
par(xaxs="i", yaxs="i")
plot(exteur,type="n", axes=FALSE, ann=FALSE)
rect(xmin(exteur),ymin(exteur),xmax(exteur),ymax(exteur), col="#9CDDF1")
rast <- raster(exteur, vals=log(expS), nrow = 276, ncol = 544)
image(rast, add=TRUE, col=pal(100))
plot(europe, border="grey25", lwd=1.2, add=TRUE)
plot(wrld[2:3,], add=TRUE, col="white")
mtext(text="Species richness",side=3,line=0.5,adj=-0.1,cex=1.25)

# MAP 2: INTERACTION PROBABILITY
# Make the map
par(xaxs="i", yaxs="i")
plot(exteur,type="n", axes=FALSE, ann=FALSE)
rect(xmin(exteur),ymin(exteur),xmax(exteur),ymax(exteur), col="#9CDDF1")
rast <- raster(exteur, vals=log(expL), nrow = 276, ncol = 544)
image(rast, add=TRUE, col=pal(100))
plot(europe, border="grey25", lwd=1.2, add=TRUE)
plot(wrld[2:3,], add=TRUE, col="white")
mtext(text="Number of interactions",side=3,line=0.5,adj=-0.1,cex=1.25)

# MAP 3: NET INTERACTION PROBABILITY
# Make the map
par(xaxs="i", yaxs="i")
plot(exteur,type="n", axes=FALSE, ann=FALSE)
rect(xmin(exteur),ymin(exteur),xmax(exteur),ymax(exteur), col="#9CDDF1")
rast <- raster(exteur, vals=expL/expS/expS, nrow = 276, ncol = 544)
image(rast, add=TRUE, col=pal(100))
plot(europe, border="grey25", lwd=1.2, add=TRUE)
plot(wrld[2:3,], add=TRUE, col="white")
mtext(text="Connectance",side=3,line=0.5,adj=-0.1,cex=1.25)

dev.copy2pdf(file = "figures/map_connectance.pdf")



