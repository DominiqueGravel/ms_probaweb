##################################################################
# Script formatting salix interaction data in a large table 
# with IDi, IDj, Xi, Xj, Xij, Lij and all of the climate variables
# Every row corresponds to a pair of species/site
# Dominique Gravel
# October 29th, 2015
##################################################################

rm(list=ls())

# Load the data
load("analysis/data/expand_data.Rdata")
load("analysis/data/salix_int.Rdata")

# Find the IDs
IDi = unique(salix_int[,5]) # i stands for hosts
IDj = unique(salix_int[,6]) # j stands for pars
IDcomm = IDj[match(IDi,IDj,nomatch=0)]
IDun = unique(c(IDi,IDj))
Si = length(IDi) 
Sj = length(IDj)

sites = salix_int[,1]
sites_ID = unique(salix_int[,1])
coords = salix_int[,2:3]
nsites = length(sites_ID)

salix_ID = unique(salix_int$to[salix_int$type == "epibiosis"])
gall_ID = unique(c(salix_int$from[salix_int$type == "epibiosis"],salix_int$to[salix_int$type != "epibiosis"]))
par_ID = unique(salix_int$from[salix_int$type != "epibiosis"])

# Pairs of interactions
SG_pairs = cbind(expand.grid(salix_ID,gall_ID),type = "SG")
GP_pairs = cbind(expand.grid(gall_ID,par_ID),type = "GP") 
pairs = rbind(SG_pairs,GP_pairs)

# Record the pair IDs
save(pairs,file="analysis/data/pairs.Rdata")

# Construct the presence/Absence matrix 
pres_mat_i = matrix(0,nr = nsites, nc = Si)
pres_mat_j = matrix(0,nr = nsites, nc = Sj)

for(i in 1:Si) {
	sites_pres_i = c(sites[which(salix_int[,5]==IDi[i])],sites[which(salix_int[,6]==IDi[i])])
	pres_mat_i[match(sites_pres_i,sites_ID),i]=1
}

for(j in 1:Sj) {
	sites_pres_j = c(sites[which(salix_int[,5]==IDj[j])],sites[which(salix_int[,6]==IDj[j])])
	pres_mat_j[match(sites_pres_j,sites_ID),j]=1	
}

# Look at co-occurrence for all pairs and record if there is an interaction
pairs0 = data.frame(
	IDi = rep(IDi, each = Sj),
	IDj = rep(IDj, times = Si),
	index_IDi = rep(c(1:Si), each = Sj),
	index_IDj = rep(c(1:Sj), times = Si)
)

pairs = data.frame(
	sites_ID = rep(sites_ID, each = Si*Sj),
	index_sites = rep(c(1:nsites), each = Si*Sj),
	IDi = rep(pairs0[,1], times = nsites),
	IDj = rep(pairs0[,2], times = nsites),
	index_IDi = rep(pairs0[,3], times = nsites),
	index_IDj = rep(pairs0[,4], times = nsites)
)

Xi = numeric(nrow(pairs))
Xj = numeric(nrow(pairs))
Lij = numeric(nrow(pairs))

expIDi = numeric(nrow(pairs))
expIDj = numeric(nrow(pairs))

link_ID = paste(salix_int[,1],salix_int[,5],salix_int[,6])
unique_link_ID = unique(link_ID)
nlink_ID = length(unique_link_ID)

pairs_ID = paste(pairs$sites_ID,pairs$IDi,pairs$IDj)

vec = numeric(length(unique_link_ID)) + 1

for(x in 1:nrow(pairs)) {
	Xi[x] = pres_mat_i[pairs$index_sites[x],pairs$index_IDi[x]]
	Xj[x] = pres_mat_j[pairs$index_sites[x],pairs$index_IDj[x]]
	if(sum(vec[unique_link_ID == pairs_ID[x]])) Lij[x] = 1
}

# Record co-occurrence for all pairs of species
Xij = Xi*Xj


###################################################
# Match with environmental data
load("analysis/data/salix_climate.Rdata")

# Subset the data to have one line per site
sub_salix_climate = salix_climate[match(sites_ID,sites),c(1,8:26)]

# Drop the NAs
sub_salix_climate = subset(sub_salix_climate,is.na(sub_salix_climate[,3])==FALSE)

# Compute the PCA
PCA = princomp(sub_salix_climate[,2:20],cor = TRUE)
summary(PCA)

# Join the scores with the climate data
sub_salix_climate = cbind(sub_salix_climate,PCA$scores)

# Match the climate data to the full table
match_pairs_climate0 = match(as.character(pairs[,1]),as.character(sub_salix_climate[,1]))
match_pairs_climate = subset(match_pairs_climate0, is.na(match_pairs_climate0) == FALSE)
expand_climate = sub_salix_climate[match_pairs_climate,]
data = data.frame(pairs = pairs[is.na(match_pairs_climate0) == FALSE,], Xi=Xi[is.na(match_pairs_climate0) == FALSE], Xj=Xj[is.na(match_pairs_climate0) == FALSE], Xij=Xij[is.na(match_pairs_climate0) == FALSE], Lij=Lij[is.na(match_pairs_climate0) == FALSE], climate = expand_climate)

save(data,file = "analysis/data/expand_data.Rdata")




