##################################################################
# Compute the likelihood for each model and each species pair
# Record the results in a summary table
# Dominique Gravel
# October 29th, 2015
##################################################################


rm(list = ls())
setwd("/Users/DGravel/Documents/Manuscripts/Inprep/ms_probaweb")

# Source the functions
source("analysis/scripts/functions/species_models.R")
source("analysis/scripts/functions/interactions_models.R")
source("analysis/scripts/functions/get_probs.R")
source("analysis/scripts/functions/collect.R")
source("analysis/scripts/functions/get_LL.R")
source("analysis/scripts/functions/fit_models.R")

# Load the data
load("analysis/data/expand_data.Rdata")
load("analysis/data/DF_split.Rdata")
load("analysis/data/pairs.Rdata")

IDi = data$pairs.IDi
IDj = data$pairs.IDj
Si = length(unique(IDi))
Sj = length(unique(IDj))
unique_IDi = unique(IDi)
unique_IDj = unique(IDj)
np = length(DF_split)

# Lists to store the models 
C2_L0 = matrix(nr = np, nc = 3)
C2_L1 = matrix(nr = np, nc = 3)
C2_L2 = matrix(nr = np, nc = 3)
C0_L2 = matrix(nr = np, nc = 3)
C1_L2 = matrix(nr = np, nc = 3)
C3_L2 = matrix(nr = np, nc = 3)

# Loop around all pairs of species
count = 1
for(x in 1:np) {

	sub_data = DF_split[[x]]

	IDs = as.matrix(sub_data[1,1:2])
	test = which(pairs[,1] == as.numeric(IDs[1]) & pairs[,2]==as.numeric(IDs[2]) | pairs[,1]==as.numeric(IDs[2]) & pairs[,2]==as.numeric(IDs[1]))
	if(length(test)!=0)	type[x] = as.character(pairs[test,3])
	
	models_C2_L0 = fit_models.apply(sub_data,selection = FALSE, funC = C2, funL = L0)
	models_C2_L1 = fit_models.apply(sub_data,selection = FALSE, funC = C2, funL = L1)
	models_C2_L2 = fit_models.apply(sub_data,selection = FALSE, funC = C2, funL = L2)
	models_C0_L2 = fit_models.apply(sub_data,selection = FALSE, funC = C0, funL = L2)
	models_C1_L2 = fit_models.apply(sub_data,selection = FALSE, funC = C1, funL = L2)
	models_C3_L2 = fit_models.apply(sub_data,selection = FALSE, funC = C3, funL = L2)

	C2_L0[x,] = get_LL.apply(models_C2_L0,sub_data)
	C2_L1[x,] =	get_LL.apply(models_C2_L1,sub_data)
	C2_L2[x,] =	get_LL.apply(models_C2_L2,sub_data)	
	C0_L2[x,] =	get_LL.apply(models_C0_L2,sub_data)
	C1_L2[x,] =	get_LL.apply(models_C1_L2,sub_data)
	C3_L2[x,] =	get_LL.apply(models_C3_L2,sub_data)		

	if(count >= 10) {
		cat(x,'\n')
		count = 1
	}
	else count = count + 1
}

# Sum by type of interaction
npars_SG = c(
	sum(C2_L0[type=="SG",2]),
	sum(C2_L1[type=="SG",2]),
	sum(C2_L2[type=="SG",2]),
	sum(C0_L2[type=="SG",2]),
	sum(C1_L2[type=="SG",2]),
	sum(C3_L2[type=="SG",2]))

npars_GP = c(
	sum(C2_L0[type=="GP",2]),
	sum(C2_L1[type=="GP",2]),
	sum(C2_L2[type=="GP",2]),
	sum(C0_L2[type=="GP",2]),
	sum(C1_L2[type=="GP",2]),
	sum(C3_L2[type=="GP",2]))

sumLL_SG = c(
	sum(C2_L0[type=="SG",1]),
	sum(C2_L1[type=="SG",1]),
	sum(C2_L2[type=="SG",1]),
	sum(C0_L2[type=="SG",1]),
	sum(C1_L2[type=="SG",1]),
	sum(C3_L2[type=="SG",1]))

sumLL_GP = c(
	sum(C2_L0[type=="GP",1]),
	sum(C2_L1[type=="GP",1]),
	sum(C2_L2[type=="GP",1]),
	sum(C0_L2[type=="GP",1]),
	sum(C1_L2[type=="GP",1]),
	sum(C3_L2[type=="GP",1]))

AIC_SG = -2*sumLL_SG + 2*npars_SG
AIC_GP = -2*sumLL_GP + 2*npars_GP

comparison_SG = cbind(sumLL_SG,npars_SG,AIC_SG)
comparison_GP = cbind(sumLL_GP,npars_GP,AIC_GP)

write.table(comparison_SG,"ms/figures/Table_model_comparison_SG.txt")
write.table(comparison_GP,"ms/figures/Table_model_comparison_GP.txt")







