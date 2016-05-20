##################################################################
# Compute the likelihood for each model and each species pair
# Record the results in a summary table
# Dominique Gravel
# October 29th, 2015
##################################################################

rm(list = ls())

# Source the functions
source("./scripts/functions/species_models.R")
source("./scripts/functions/interactions_models.R")
source("./scripts/functions/get_probs.R")
source("./scripts/functions/get_LL.R")
source("./scripts/functions/fit_models.R")

# Load the data
load("./data/expand_data.Rdata")
load("./data/DF_split.Rdata")
load("./data/pairs.Rdata")

IDi = data$pairs.IDi
IDj = data$pairs.IDj
Si = length(unique(IDi))
Sj = length(unique(IDj))
unique_IDi = unique(IDi)
unique_IDj = unique(IDj)
np = length(DF_split)
type = numeric(np)

# Lists to store the models 
C2_L0 = matrix(nr = np, nc = 2)
C2_L1 = matrix(nr = np, nc = 2)
C2_L2 = matrix(nr = np, nc = 2)
C0_L2 = matrix(nr = np, nc = 2)
C1_L2 = matrix(nr = np, nc = 2)
C3_L2 = matrix(nr = np, nc = 2)

# Loop around all pairs of species
count = 1
for(x in 1:np) {

	sub_data = DF_split[[x]]

	IDs = as.matrix(sub_data[1,1:2])
	test = which(pairs[,2] == as.character(IDs[1]) & pairs[,3]==as.character(IDs[2]) | pairs[,2]==as.character(IDs[2]) & pairs[,3]==as.character(IDs[1]))
	if(length(test)!=0)	type[x] = unique(as.character(pairs[test,4]))
	
	# Compute the models
	models_C2_L0 = fit_models(sub_data, selection = FALSE, funC = C2, funL = L0)
	models_C2_L1 = fit_models(sub_data, selection = FALSE, funC = C2, funL = L1)
	models_C2_L2 = fit_models(sub_data, selection = FALSE, funC = C2, funL = L2)
	models_C0_L2 = fit_models(sub_data, selection = FALSE, funC = C0, funL = L2)
	models_C1_L2 = fit_models(sub_data, selection = FALSE, funC = C1, funL = L2)
	models_C3_L2 = fit_models(sub_data, selection = FALSE, funC = C3, funL = L2)

	# Compute the likelihood
	C2_L0[x,] = get_LL(models_C2_L0, sub_data)
	C2_L1[x,] =	get_LL(models_C2_L1, sub_data)
	C2_L2[x,] =	get_LL(models_C2_L2, sub_data)	
	C0_L2[x,] =	get_LL(models_C0_L2, sub_data)
	C1_L2[x,] =	get_LL(models_C1_L2, sub_data)
	C3_L2[x,] =	get_LL(models_C3_L2, sub_data)		

	if(count >= 10) {
		cat(x,'\n')
		count = 1
	}
	else count = count + 1
}

# Sum by type of interaction
npars_SG = c(
	sum(C2_L0[type=="SH",2]),
	sum(C2_L1[type=="SH",2]),
	sum(C2_L2[type=="SH",2]),
	sum(C0_L2[type=="SH",2]),
	sum(C1_L2[type=="SH",2]),
	sum(C3_L2[type=="SH",2]))

npars_GP = c(
	sum(C2_L0[type=="HP",2]),
	sum(C2_L1[type=="HP",2]),
	sum(C2_L2[type=="HP",2]),
	sum(C0_L2[type=="HP",2]),
	sum(C1_L2[type=="HP",2]),
	sum(C3_L2[type=="HP",2]))

sumLL_SG = c(
	sum(C2_L0[type=="SH",1]),
	sum(C2_L1[type=="SH",1]),
	sum(C2_L2[type=="SH",1]),
	sum(C0_L2[type=="SH",1]),
	sum(C1_L2[type=="SH",1]),
	sum(C3_L2[type=="SH",1]))

sumLL_GP = c(
	sum(C2_L0[type=="HP",1]),
	sum(C2_L1[type=="HP",1]),
	sum(C2_L2[type=="HP",1]),
	sum(C0_L2[type=="HP",1]),
	sum(C1_L2[type=="HP",1]),
	sum(C3_L2[type=="HP",1]))

AIC_SG = -2*sumLL_SG + 2*npars_SG
AIC_GP = -2*sumLL_GP + 2*npars_GP

comparison_SG = cbind(sumLL_SG,npars_SG,AIC_SG)
comparison_GP = cbind(sumLL_GP,npars_GP,AIC_GP)

write.table(comparison_SG,"./tables/Table_model_comparison_SG.txt")
write.table(comparison_GP,"./tables/Table_model_comparison_GP.txt")