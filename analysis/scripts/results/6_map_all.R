##################################################################
# Run the models across Europe to compute the expected species
# richness and the expected number of links
# Dominique Gravel
# October 29th, 2015
##################################################################

rm(list=ls())
library(raster)
library(sp)
library(rgdal)

setwd("/Users/DGravel/Documents/Manuscripts/Inprep/ms_probaweb")

load("analysis/data/DF_split.Rdata")
load("analysis/data/expand_data.Rdata")
load("analysis/data/pairs.Rdata")

IDi = data$pairs.from
IDj = data$pairs.to
Type = data$pairs.type
Si = length(unique(IDi))
Sj = length(unique(IDj))
nobs = 374
np  = length(DF_split)

# Load the stuff to make the raster object
wrld <- readOGR("analysis/data/map", layer="level1")
europe <- wrld[1,]
exteur <- extent(europe)
xmin(exteur) = 0
xmax(exteur) = 25
ymin(exteur) = 40
ymax(exteur) = 70

# Extract the climatic variables
clim = crop(getData("worldclim", var="bio", res=10), exteur)
mapT = getValues(clim$bio1)
mapPP = getValues(clim$bio12)
ncell = length(T)
map_E = data.frame(T = mapT, PP = mapPP, T2 = mapT^2, PP2 = mapPP^2)

############################################################
# Extract the local networks over Europe
############################################################

#############
# Compute the distribution of salix richness
data_S = subset(data[,c(1,2,5,10,21)], data$pairs.type == "SH")
ID_S = unique(data_S$pairs.from)
exp_S = numeric(nrow(map_E))
for(i in 1:length(ID_S)) {

	# Subset the data
	sub_data = unique(data_S[data_S$pairs.from == ID_S[i],])
	Xi = sub_data$Xi
	T = sub_data$climate.bio1
	PP = sub_data$climate.bio12
	T2 = T*T
	PP2 = PP*PP

	# Run GLM
	lmXi = glm(Xi ~ T + PP + T2 + PP2, family = "binomial")

	# Compile predictions
	predXi = predict(lmXi, type = "response", newdata = map_E)
	exp_S = exp_S + predXi
}

#############
# Compute the distribution of host richness
data_H = subset(data[,c(1,3,6,10,21)], data$pairs.type == "SH")
ID_H = unique(data_H$pairs.to)
exp_H = numeric(nrow(map_E))
# Loop around all species
for(i in 1:length(ID_H)) {

	# Subset the data
	sub_data = unique(data_H[data_H$pairs.to == ID_H[i],])
	Xi = sub_data$Xj
	T = sub_data$climate.bio1
	PP = sub_data$climate.bio12
	T2 = T*T
	PP2 = PP*PP

	# Run GLM
	lmXi = glm(Xi ~ T + PP + T2 + PP2, family = "binomial")

	# Compile predictions
	predXi = predict(lmXi, type = "response", newdata = map_E)
	exp_H = exp_H + predXi
}

#############
# Compute the distribution of parasitoid richness
data_P = subset(data[,c(1,3,6,10,21)], data$pairs.type == "HP")
ID_P = unique(data_P$pairs.to)
exp_P = numeric(nrow(map_E))
for(i in 1:length(ID_P)) {

	# Subset the data
	sub_data = unique(data_P[data_P$pairs.to == ID_P[i],])

	Xi = sub_data$Xj
	T = sub_data$climate.bio1
	PP = sub_data$climate.bio12
	T2 = T*T
	PP2 = PP*PP

	# Run GLM
	lmXi = glm(Xi ~ T + PP + T2 + PP2, family = "binomial")

	# Compile predictions
	predXi = predict(lmXi, type = "response", newdata = map_E)
	exp_P = exp_P + predXi
}

#############
# Compute the distribution of Salix-
data_P = subset(data[,c(1,3,6,10,21)], data$pairs.type == "HP")
ID_P = unique(data_P$pairs.to)
exp_P = numeric(nrow(map_E))
for(i in 1:length(ID_P)) {

	# Subset the data
	sub_data = unique(data_P[data_P$pairs.to == ID_P[i],])

	Xi = sub_data$Xj
	T = sub_data$climate.bio1
	PP = sub_data$climate.bio12
	T2 = T*T
	PP2 = PP*PP

	# Run GLM
	lmXi = glm(Xi ~ T + PP + T2 + PP2, family = "binomial")

	# Compile predictions
	predXi = predict(lmXi, type = "response", newdata = map_E)
	exp_P = exp_P + predXi
}





#############
# Compute the expected link density
exp_SH = numeric(ncell)
exp_HP = numeric(ncell)

for(pair in 1:np) {

	data = DF_split[[pair]]
	Xi = data$Xi
	Xj = data$Xj
	Xij = data$Xij
	Lij = data$Lij
	E = data$E	

	# Compute the models
	if(sum(Xij)!=0) {
		lm_Xij = glm(Xij ~ T+PP+T2+PP2, family = "binomial", data = E)
		expLij = sum(Lij)/sum(Xij)
	}

	else {
			expXij = 0
			expLij = 0
	}

	# Compute expected probabilities over the gradient
	if(sum(Xij)!=0) {
		expXij = predict(lm_Xij, type="response", newdata=map_E)
		# Compute species richness and link density
		if(pairs[pair,3] == "SH") exp_SH = exp_SH + expLij*expXij
		else exp_HP = exp_HP + expLij*expXij
	}

	cat(pair,'\n')
}



write.table(cbind(expS,expL),"ms/figures/map_data.txt")

############################################################
# Plot network properties
############################################################

quartz(height = 7, width = 8)
par(mfrow = c(2,3),mar = c(2,2,3,2))

# Map parameters
library("colorRamps")
library("RColorBrewer")
pal <-colorRampPalette(rev(brewer.pal(11,"RdYlBu"))) # Initialized Colors ramp palette, see here: http://colorbrewer2.org/

# MAP 1: Salix richness
# Make the map
par(xaxs="i", yaxs="i")
plot(exteur,type="n", axes=FALSE, ann=FALSE)
rect(xmin(exteur),ymin(exteur),xmax(exteur),ymax(exteur), col="#9CDDF1")
rast <- raster(exteur, vals = log(exp_S), nrow = 180, ncol = 150)
image(rast, add=TRUE, col=pal(100))
plot(europe, border="grey25", lwd=1.2, add=TRUE)
mtext(text="Salix richness",side=3,line=0.5,adj=-0.0,cex=1.25)

# MAP 2: host richness
# Make the map
par(xaxs="i", yaxs="i")
plot(exteur,type="n", axes=FALSE, ann=FALSE)
rect(xmin(exteur),ymin(exteur),xmax(exteur),ymax(exteur), col="#9CDDF1")
rast <- raster(exteur, vals = log(exp_H), nrow = 180, ncol = 150)
image(rast, add=TRUE, col=pal(100))
plot(europe, border="grey25", lwd=1.2, add=TRUE)
mtext(text="Gall richness",side=3,line=0.5,adj=-0.0,cex=1.25)

# MAP 3: parasitoid richness
# Make the map
par(xaxs="i", yaxs="i")
plot(exteur,type="n", axes=FALSE, ann=FALSE)
rect(xmin(exteur),ymin(exteur),xmax(exteur),ymax(exteur), col="#9CDDF1")
rast <- raster(exteur, vals = log(exp_P), nrow = 180, ncol = 150)
image(rast, add=TRUE, col=pal(100))
plot(europe, border="grey25", lwd=1.2, add=TRUE)
mtext(text="Parasitoid richness",side=3,line=0.5,adj=-0.0,cex=1.25)

# MAP 4: S-H link density
# Make the map
par(xaxs="i", yaxs="i")
plot(exteur,type="n", axes=FALSE, ann=FALSE)
rect(xmin(exteur),ymin(exteur),xmax(exteur),ymax(exteur), col="#9CDDF1")
rast <- raster(exteur, vals = log(exp_SH), nrow = 180, ncol = 150)
image(rast, add=TRUE, col=pal(100))
plot(europe, border="grey25", lwd=1.2, add=TRUE)
mtext(text="Link density - 
	Salix & Galls",side=3,line=0.5,adj=-0.0,cex=1.25)

# MAP 5: H-P link density
# Make the map
par(xaxs="i", yaxs="i")
plot(exteur,type="n", axes=FALSE, ann=FALSE)
rect(xmin(exteur),ymin(exteur),xmax(exteur),ymax(exteur), col="#9CDDF1")
rast <- raster(exteur, vals = log(exp_HP), nrow = 180, ncol = 150)
image(rast, add=TRUE, col=pal(100))
plot(europe, border="grey25", lwd=1.2, add=TRUE)
mtext(text="Link density - 
	Galls & Parasitoids",side=3,line=0.5,adj=-0.0,cex=1.25)

# MAP 6: S-H link density
# Make the map
par(xaxs="i", yaxs="i")
plot(exteur,type="n", axes=FALSE, ann=FALSE)
rect(xmin(exteur),ymin(exteur),xmax(exteur),ymax(exteur), col="#9CDDF1")
C = (exp_SH + exp_HP)/(exp_S*exp_H + exp_H*exp_P)
rast <- raster(exteur, vals = log(exp_P), nrow = 180, ncol = 150)
image(rast, add=TRUE, col=pal(100))
plot(europe, border="grey25", lwd=1.2, add=TRUE)
mtext(text="Connectance", side=3, line=0.5,adj=-0.0,cex=1.25)

dev.copy2pdf(file = "ms/figures/map_connectance.pdf")






