rm(list=ls())
library(raster)
library(sp)
library(rgdal)
source("scripts/functions/collect.R")
source("scripts/functions/species_models.R")
source("scripts/functions/interactions_models.R")
source("scripts/functions/get_LL.R")
source("scripts/functions/get_probs.R")
source("scripts/functions/fit_models.R")
load("data/DF_split.Rdata")
load("data/expand_data.Rdata")

original_data = data
IDi = original_data$pairs.IDi
IDj = original_data$pairs.IDj
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
T = getValues(clim$bio1)
PP = getValues(clim$bio12)


############################################################
# Extract the local networks over Europe
############################################################

ncell = length(T)
map_E = data.frame(T = T, PP = PP, T2 = T^2, PP2 = PP^2)

# Loop around all species pairs
PXi_E = list()
PXj_E = list()
PXij_E = list()
PLij_E = list()

PXi_C = list()
PXj_C = list()
PXij_C = list()
PLij_C = list()

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
	lm_i = glm(Xi ~ T*PP*T2*PP2, family = "binomial", data = E)
	lm_j = glm(Xj ~ T*PP*T2*PP2, family = "binomial", data = E)
	if(i != j) {
		if(sum(Xij!=0)) {
			subE = subset(E,Xij == 1)
			subLij = subset(Lij, Xij == 1)
			lm_L = glm(subLij ~ T*PP*T2*PP2, family = "binomial", data = subE)
		}
			else lm_L = NULL
	}
		else lm_L = NULL 

	# Compute expected probabilities over the gradient
	PXi_E[[i]] = predict(lm_i, type="response",newdata=map_E)
	PXj_E[[j]] = predict(lm_j, type="response",newdata=map_E)
	if(is.null(lm_L)) PLij_E[[pair]] = 0
		else PLij_E[[pair]] = predict(lm_L, type="response", newdata=map_E)
	PXij_E[[pair]] = PXi_E[[i]]*PXj_E[[j]]

	# Compute the constant probabilities
	PXi_C[[i]] = sum(data$Xi)/nobs 
	PXj_C[[j]] = sum(data$Xj)/nobs 
	PXij_C[[pair]] = sum(data$Xij)/nobs

	if(PXij_C[[pair]]!=0) PLij_C[[pair]] = sum(data$Lij)/PXij_C[[pair]]
	else PLij_C[[pair]] = 0

	cat(pair,'\n')
}

# Extract network caracteristics for each locality
pick = function(X,index) {
	if(is.null(X)) NULL
	else if(length(X)>1) X[index] 
	else 0
}

mat_PXi_E = matrix(nr = nsteps, nc = Si)
mat_PXj_E = matrix(nr = nsteps, nc = Sj)
mat_PXij_E = matrix(nr = nsteps, nc = np)
mat_PLij_E = matrix(nr = nsteps, nc = np)

IDi = unique(as.numeric(as.vector(IDi)))
IDj = unique(as.numeric(as.vector(IDj)))

for(x in 1:ncell) {
	mat_PXi_E[x,] = unlist(lapply(PXi_E,FUN=pick,index=x))
	mat_PXj_E[x,] = unlist(lapply(PXj_E,FUN=pick,index=x))
	mat_PXij_E[x,] = unlist(lapply(PXij_E,FUN=pick,index=x))
	mat_PLij_E[x,] = unlist(lapply(PLij_E,FUN=pick,index=x))		
	cat(x,'\n')
}

map_exp_Si_E = apply(mat_PXi_E,1,sum)
map_exp_Sj_E = apply(mat_PXj_E,1,sum)
map_exp_L_E = apply(mat_PXij_E*mat_PLij_E,1,sum)
map_Co_1 = map_exp_L_E/map_exp_Si_E/map_exp_Sj_E

write.table(map_exp_Si_E,file = "data/map_exp_Si_E.txt")
write.table(map_exp_Sj_E,file = "data/map_exp_Sj_E.txt")
write.table(map_exp_L_E,file = "data/map_exp_L_E.txt")
write.table(map_exp_Co_1,file = "data/map_exp_Co_1.txt")

############################################################
# Do the same over a temperature gradient only and different hypotheses
############################################################
gradient_T = seq(0,180,1)
mean_PP = mean(PP,na.rm=TRUE)
gradient_E = data.frame(T = gradient_T,PP = mean_PP, T2 = gradient_T^2, PP2 = mean_PP^2)
nsteps = length(gradient_T)

# Loop around all species pairs
PXi_E = list()
PXj_E = list()
PXij_E = list()
PLij_E = list()

PXi_C = list()
PXj_C = list()
PXij_C = list()
PLij_C = list()

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
	lm_i = glm(Xi ~ T*PP*T2*PP2, family = "binomial", data = E)
	lm_j = glm(Xj ~ T*PP*T2*PP2, family = "binomial", data = E)
	if(i != j) {
		if(sum(Xij!=0)) {
			subE = subset(E,Xij == 1)
			subLij = subset(Lij, Xij == 1)
			lm_L = glm(subLij ~ T*PP*T2*PP2, family = "binomial", data = subE)
		}
			else lm_L = NULL
	}
		else lm_L = NULL 

	# Compute expected probabilities over the gradient
	PXi_E[[i]] = predict(lm_i, type="response",newdata=gradient_E)
	PXj_E[[j]] = predict(lm_j, type="response",newdata=gradient_E)
	if(is.null(lm_L)) PLij_E[[pair]] = 0
		else PLij_E[[pair]] = predict(lm_L, type="response", newdata=gradient_E)
	PXij_E[[pair]] = PXi_E[[i]]*PXj_E[[j]]

	# Compute the constant probabilities
	PXi_C[[i]] = sum(data$Xi)/nobs 
	PXj_C[[j]] = sum(data$Xj)/nobs 
	PXij_C[[pair]] = sum(data$Xij)/nobs

	if(PXij_C[[pair]]!=0) PLij_C[[pair]] = sum(data$Lij)/PXij_C[[pair]]
	else PLij_C[[pair]] = 0

	cat(pair,'\n')
}

# Extract network caracteristics for each locality
pick = function(X,index) {
	if(is.null(X)) NULL
	else if(length(X)>1) X[index] 
	else 0
}

mat_PXi_E = matrix(nr = nsteps, nc = Si)
mat_PXj_E = matrix(nr = nsteps, nc = Sj)
mat_PXij_E = matrix(nr = nsteps, nc = np)
mat_PLij_E = matrix(nr = nsteps, nc = np)

IDi = unique(as.numeric(as.vector(IDi)))
IDj = unique(as.numeric(as.vector(IDj)))

for(x in 1:nsteps) {
	mat_PXi_E[x,] = unlist(lapply(PXi_E,FUN=pick,index=x))
	mat_PXj_E[x,] = unlist(lapply(PXj_E,FUN=pick,index=x))
	mat_PXij_E[x,] = unlist(lapply(PXij_E,FUN=pick,index=x))
	mat_PLij_E[x,] = unlist(lapply(PLij_E,FUN=pick,index=x))		
	cat(x,'\n')
}

grad_exp_Si_E = apply(mat_PXi_E,1,sum)
grad_exp_Sj_E = apply(mat_PXj_E,1,sum)
grad_exp_L_E = apply(mat_PXij_E*mat_PLij_E,1,sum)
grad_Co_1 = grad_exp_L_E/grad_exp_Si_E/grad_exp_Sj_E

write.table(grad_exp_Si_E,file = "data/grad_exp_Si_E.txt")
write.table(grad_exp_Sj_E,file = "data/grad_exp_Sj_E.txt")
write.table(grad_exp_L_E,file = "data/grad_exp_L_E.txt")
write.table(grad_exp_Co_1,file = "data/grad_exp_Co_1.txt")

PXi_C = unlist(PXi_C)
PXj_C = unlist(PXj_C)
PLij_C = unlist(PLij_C)
grad_exp_Si_C = sum(PXi_C)
grad_exp_Sj_C = sum(PXi_C)
grad_Co_2 = apply(PXij_E*PLij_C,1,sum)/exp_Si_E/exp_Sj_E
grad_Co_3 = apply(PXij_C*PLij_E,1,sum)/exp_Si_C/exp_Sj_C

############################################################
# Plot network properties
############################################################

dev.new(height = 3.5, width = 10)
par(mfrow = c(1,3))

# Convert network properties into a raster
rast_C <- raster(exteur, vals=map_Co_1)

# Make the map
par(xaxs="i", yaxs="i")
plot(exteur,type="n", axes=FALSE, ann=FALSE)
rect(xmin(exteur),ymin(exteur),xmax(exteur),ymax(exteur), col="#9CDDF1")
image(rast_C, add=TRUE, col=terrain.colors(100))
plot(europe, border="grey25", lwd=1.2, add=TRUE)
plot(wrld[2:3,], add=TRUE, col=c("#F39D68","#F8D762"))

# Add on the side the slices through the temperature gradient
plot(gradient_T, grad_exp_Si_E, type = "l", xlab = "Annual mean temperature", ylab = "Species richness")
lines(gradient_T,grad_exp_Sj_E, lty = 2)
legend("topright",legend = c("Victim","Ennemy"),lty = c(1,2,bty = "n")

plot(gradient_T, grad_Co_1, ylim = c(0,0.3), type = "l", xlab = "Annual mean temperature", ylab = "Connectance")
lines(gradient_T, grad_Co_2, lty = 2)
lines(gradient_T, grad_Co_3, lty = 3)
legend("topright",legend = c("Integrated","Constant co-occurrence","Constant interactions"),lty = c(1,2,3),bty = "n")






