##################################################################
# Compute the likelihood for each model and each species pair
# Record the results in a summary table
# Dominique Gravel
##################################################################

rm(list = ls())

# Source the functions
source("./scripts/functions/species_models.r")
source("./scripts/functions/interactions_models.r")
source("./scripts/functions/get_probs.r")
source("./scripts/functions/get_LL.r")
source("./scripts/functions/fit_models.r")

# Load & prepare the data
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

# Matrices to store the results
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
	type[x] = unique(as.character(pairs[x,3]))
	
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
npars_SH = c(
	sum(C2_L0[type=="SH",2]),
	sum(C2_L1[type=="SH",2]),
	sum(C2_L2[type=="SH",2]),
	sum(C0_L2[type=="SH",2]),
	sum(C1_L2[type=="SH",2]),
	sum(C3_L2[type=="SH",2]))

npars_HP = c(
	sum(C2_L0[type=="HP",2]),
	sum(C2_L1[type=="HP",2]),
	sum(C2_L2[type=="HP",2]),
	sum(C0_L2[type=="HP",2]),
	sum(C1_L2[type=="HP",2]),
	sum(C3_L2[type=="HP",2]))

sumLL_SH = c(
	sum(C2_L0[type=="SH",1]),
	sum(C2_L1[type=="SH",1]),
	sum(C2_L2[type=="SH",1]),
	sum(C0_L2[type=="SH",1]),
	sum(C1_L2[type=="SH",1]),
	sum(C3_L2[type=="SH",1]))

sumLL_HP = c(
	sum(C2_L0[type=="HP",1]),
	sum(C2_L1[type=="HP",1]),
	sum(C2_L2[type=="HP",1]),
	sum(C0_L2[type=="HP",1]),
	sum(C1_L2[type=="HP",1]),
	sum(C3_L2[type=="HP",1]))

AIC_SH = -2*sumLL_SH + 2*npars_SH
AIC_HP = -2*sumLL_HP + 2*npars_HP

comparison_SH = cbind(sumLL_SH,npars_SH,AIC_SH)
comparison_HP = cbind(sumLL_HP,npars_HP,AIC_HP)

write.table(comparison_SH,"./tables/Table_model_comparison_SH.txt")
write.table(comparison_HP,"./tables/Table_model_comparison_HP.txt")