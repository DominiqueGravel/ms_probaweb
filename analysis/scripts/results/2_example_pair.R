##################################################################
# Script for Figure X,
# Illustrating the computation of interaction probabilities for a 
# pair of species
# Dominique Gravel
# October 29th, 2015
##################################################################

rm(list = ls())

# Load some functions
source("analysis/scripts/functions/species_models.R")
source("analysis/scripts/functions/interactions_models.R")
source("analysis/scripts/functions/get_LL.R")
source("analysis/scripts/functions/get_probs.R")
source("analysis/scripts/functions/fit_models.R")

# Load the data
load("analysis/data/DF_split.Rdata")
load("analysis/data/expand_data.Rdata")

#########################################################
# Fit the different models
#########################################################

# Find the nb of links per pair
nL  = numeric(length(DF_split))
nXi = numeric(length(DF_split))
nXj = numeric(length(DF_split))
nXij = numeric(length(DF_split))
for(x in 1:length(DF_split)) {
	nL[x] = sum(DF_split[[x]]$Lij)
	nXi[x] = sum(DF_split[[x]]$Xi)
	nXj[x] = sum(DF_split[[x]]$Xj)
	nXij[x] = sum(DF_split[[x]]$Xij)
}

cbind(nL,nXi,nXj,nXij,nL,nL/nXij)[nL/nXij<0.7 & nL > 10,]

# Subset the data

#pair_index = which(nL == 13 & nXij == 22) 
pair_index = which(nL == 21 & nXij == 38) 
#pair_index = which(nL == 13 & nXij == 22) 
#pair_index = which(nL == 30 & nXij == 41) 
#pair_index = which(nL == 9 & nXij == 19) 

#pair_index = which(nL == 14 & nXij == 25)
#pair_index = which(nL == 16 & nXij == 26)
#pair_index = which(nL == 15 & nXij == 24)
#pair_index = which(nL == 10 & nXij == 24)

#pair_index = which(nL == 16 & nXij == 23)
#pair_index = which(nL == 11 & nXij == 18)[1]



data = DF_split[[pair_index]]
data$E = data.frame(T = data$E$T/12, PP = data$E$PP/1000, T2 = data$E$T2/12^2, PP2 = data$E$PP2/1000^2)
sum(data$Lij)
sum(data$Xij)

# Pick the model
models_C2_L0 = fit_models.apply(data, selection = FALSE, funC = C2, funL = L0)
models_C2_L1 = fit_models.apply(data, selection = FALSE, funC = C2, funL = L1)
models_C2_L2 = fit_models.apply(data, selection = FALSE, funC = C2, funL = L2)
models_C0_L2 = fit_models.apply(data, selection = FALSE, funC = C0, funL = L2)
models_C1_L2 = fit_models.apply(data, selection = FALSE, funC = C1, funL = L2)
models_C3_L2 = fit_models.apply(data, selection = FALSE, funC = C3, funL = L2)

# Compute the LL
LL_C2_L0 = get_LL.apply(models_C2_L0,data)
LL_C2_L1 = get_LL.apply(models_C2_L1,data)
LL_C2_L2 = get_LL.apply(models_C2_L2,data)
LL_C0_L2 = get_LL.apply(models_C0_L2,data)
LL_C1_L2 = get_LL.apply(models_C1_L2,data)
LL_C3_L2 = get_LL.apply(models_C3_L2,data)

# Collect the results
LL = c(
	LL_C2_L0[1],
	LL_C2_L1[1],
	LL_C2_L2[1],
	LL_C0_L2[1],
	LL_C1_L2[1],
	LL_C3_L2[1]	
	)

npars = c(
	LL_C2_L0[2],
	LL_C2_L1[2],
	LL_C2_L2[2],
	LL_C0_L2[2],
	LL_C1_L2[2],
	LL_C3_L2[2]	
	)

AIC = -2*LL + 2*npars
AIC

# Write the results in a table
write.table(cbind(LL,npars,AIC), file = "ms/figures/table1.txt") 


#########################################################
# Plot the results
#########################################################

# Compute predicted values for each observation
models = models_C2_L2
probs = get_probs(modelC=models$modelC, modelL=models$modelL, newE=data$E)

# Compute predicted values for the environmental space
nsteps = 250
seqT = seq(min(data$E$T,na.rm=T),max(data$E$T,na.rm=T),diff(range(data$E$T))/(nsteps-1))
seqPP = seq(min(data$E$PP,na.rm=T),max(data$E$PP,na.rm=T),diff(range(data$E$PP))/(nsteps-1))

expE = expand.grid(seqT,seqPP)
expT = expE[,1]
expPP = expE[,2]
expT2 = expT^2
expPP2 = expPP^2

newE = data.frame(T = expT,T2 = expT2, PP = expPP, PP2 = expPP2)
probs = get_probs(modelC=models$modelC, modelL=models$modelL, newE=newE)

# Plot the results
PXijmat = matrix(probs$PXij, nr = nsteps, nc = nsteps, byrow = FALSE)
PLijmat = matrix(probs$PLij, nr = nsteps, nc = nsteps, byrow = FALSE)
PLijXijmat = PXijmat*PLijmat
PLijXijmat[PLijXijmat=="NaN"] = 0

quartz(height = 3.5, width = 10)
par(mfrow = c(1,3),mar = c(5,6,2.5,1))

image(seqT,seqPP,1-PXijmat, xlab = "Annual mean temperature",ylab = "Annual precipitation",cex.axis = 1.25, cex.lab = 1.5, col = gray(seq(0.3,1,1/1000)),xlim = 1.0*range(seqT),ylim = 1.0*range(seqPP))
points(data$E$T[data$Xij==1],data$E$PP[data$Xij==1],pch = 19)
points(data$E$T[data$Xi==1 & data$Xij==0],data$E$PP[data$Xi==1 & data$Xij!=1],pch = 3)
points(data$E$T[data$Xj==1 & data$Xij==0],data$E$PP[data$Xj==1 & data$Xij!=1],pch = 3)
points(data$E$T[data$Xi==0 & data$Xj==0],data$E$PP[data$Xi==0 & data$Xj==0],pch = 1)
mtext(text=expression(P(X[i],X[j])),side=3,line=0.5,adj=-0.1,cex=1.25)

image(seqT,seqPP,1-PLijmat, xlab = "Annual mean temperature",ylab = "Annual precipitation",cex.axis = 1.25, cex.lab = 1.5, col = gray(seq(0.5,1,1/1000)),xlim = 1.0*range(seqT),ylim = 1.0*range(seqPP))
points(data$E$T[data$Lij==1],data$E$PP[data$Lij==1],pch = 19)
points(data$E$T[data$Lij==0 & data$Xij==1],data$E$PP[data$Lij==0 & data$Xij==1],pch = 1)
mtext(text=expression(P(L[ij])),side=3,line=0.5,adj=-0.1,cex=1.25)

image(seqT,seqPP,1-PLijXijmat, xlab = "Annual mean temperature",ylab = "Annual precipitation",cex.axis = 1.25, cex.lab = 1.5, col = gray(seq(0.5,1,1/1000)),xlim = 1.0*range(seqT),ylim = 1.0*range(seqPP))
points(data$E$T[data$Lij==1],data$E$PP[data$Lij==1],pch = 19)
points(data$E$T[data$Lij==0 & data$Xij==1],data$E$PP[data$Lij==0 & data$Xij==1],pch = 1)
mtext(text=expression(P(L[ij],X[i],X[j])),side=3,line=0.5,adj=-0.1,cex=1.25)

dev.copy2pdf(file = "ms/figures/example_pair.pdf")


