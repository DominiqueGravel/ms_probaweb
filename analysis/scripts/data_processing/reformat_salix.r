##################################################################
# Script formatting salix interaction data in a large table 
# with IDi, IDj, Xi, Xj, Xij, Lij and all of the climate variables
# Every row corresponds to a pair of species/site
# Dominique Gravel
# October 29th, 2015
##################################################################

rm(list=ls())

# Load the data
datakop <- read.csv("data/kopeleke_data.csv",dec=",")
HP = datakop[,c(9,10,11,1,19,32)]
#HP = subset(HP, HP[,6]!="no parasitoid or inquiline")
SH = datakop[,c(9,10,11,1,29,19)]
names(HP) = c("site","lat","lon","year","from","to")
names(SH) = c("site","lat","lon","year","from","to")
SH[,5] = paste("Salix",SH[,5])
salix_int = rbind(HP,SH)
save(salix_int, file = "./data/salix_int_direct.Rdata")

# Get the IDs
sites = as.character(salix_int[,1])
sites_ID = as.character(unique(salix_int[,1]))
coords = salix_int[,2:3]
n_sites = length(sites_ID)
ID_S = unique(SH[,5])
ID_H = as.character(unique(SH[,6]))
ID_P = as.character(unique(HP[,6]))
S_S = length(ID_S)
S_H = length(ID_H)
S_P = length(ID_P)

# Pairs of interactions
SH_pairs = cbind(expand.grid(ID_S, ID_H), type = "SG")
HP_pairs = cbind(expand.grid(ID_H, ID_P), type = "GP") 
pairs = rbind(SH_pairs, HP_pairs)

# Record the pair IDs
save(pairs, file="./data/pairs.Rdata")

# Look at co-occurrence for pairs of salix ang galls
pairs0_SH = data.frame(
	ID_S = rep(ID_S, each = S_H),
	ID_H = rep(ID_H, times = S_S))

pairs_SH = data.frame(
	sites_ID = rep(sites_ID, each = S_S * S_H),
	ID_S = rep(pairs0_SH[,1], times = n_sites),
	ID_H = rep(pairs0_SH[,2], times = n_sites), 
	type = rep("SH", times = n_sites*S_S*S_H)
)
rm(pairs0_SH)
XS = numeric(nrow(pairs_SH))
XH = numeric(nrow(pairs_SH))
LSH = numeric(nrow(pairs_SH))

# if(parallel == TRUE){
	
# 	# Init parrallele
# 	#0. Packages loading
# 	library(doParallel)

# 	#1. Init the cluster of SOCKETS
# 	# Get the number of cores by node
# 	nCores <- detectCores()

# 	# Open the cluster
# 	cl <- makeCluster(nCores)
# 	# We replicated the addresses nCores times.
# 	registerDoParallel(cl)

# 	# Loop around all observations of interactions
# 	id_xs <- foreach(x=1:nrow(SH),.combine=c) %dopar% {
# 		# Get the sites where the salix is present
# 		 return(which(pairs_SH$sites_ID == SH[x,1] & pairs_SH$ID_S == SH[x,5]))
# 	}
	
# 	id_xh <- foreach(x=1:nrow(SH),.combine=c) %dopar% {
# 		# Get the sites where the salix is present
# 		 return(which(pairs_SH$sites_ID == SH[x,1] & pairs_SH$ID_H == SH[x,6]))
# 	}
	
# 	id_lsh <- foreach(x=1:nrow(SH),.combine=c) %dopar% {
# 		# Get the sites where the salix is present
# 		 return(which(pairs_SH$sites_ID == SH[x,1] & pairs_SH$ID_S == SH[x,5] & pairs_SH$ID_H == SH[x,6]))
# 	}
	
# 	# convert matches to 1
# 	XS[id_xs] = 1
# 	XH[id_xh] = 1
# 	LSH[id_lsh] = 1
	
# } else{
	
	# Loop around all observations of interactions
	for(x in 1:nrow(SH)) {

		# Get the sites where the salix is present
		XS[which(pairs_SH$sites_ID == SH[x,1] & pairs_SH$ID_S == SH[x,5])] = 1

		# Get the sites where the gall is present
		XH[which(pairs_SH$sites_ID == SH[x,1] & pairs_SH$ID_H == SH[x,6])] = 1
			
		# Get the sites where the gall and the salix interact
		LSH[which(pairs_SH$sites_ID == SH[x,1] & pairs_SH$ID_S == SH[x,5] & pairs_SH$ID_H == SH[x,6])] = 1

	}
#}


data_SH = cbind(pairs_SH,XS,XH,LSH)
save(data_SH, file = "./data/SH.Rdata")

##############################################

# Look at co-occurrence for pairs of galls and parasitoids
pairs0_HP= data.frame(
	ID_H = rep(ID_H, each = S_P),
	ID_P = rep(ID_P, times = S_H))

pairs_HP = data.frame(
	sites_ID = rep(sites_ID, each = S_H * S_P),
	ID_H = rep(pairs0_HP[,1], times = n_sites),
	ID_P = rep(pairs0_HP[,2], times = n_sites), 
	type = rep("HP", n_sites*S_H*S_P))
rm(pairs0_HP)
XH2 = numeric(nrow(pairs_HP))
XP = numeric(nrow(pairs_HP))
LHP = numeric(nrow(pairs_HP))

# Loop around all observations of interactions
# if(parallel==TRUE){
	
# 	id_xh2 <- foreach(x=1:nrow(HP),.combine=c) %dopar% {
# 		# Get the sites where the salix is present
# 		 return(which(pairs_HP$sites_ID == HP[x,1] & pairs_HP$ID_H == HP[x,5]))
# 	}
	
# 	id_xp <- foreach(x=1:nrow(HP),.combine=c) %dopar% {
# 			# Get the sites where the gall is present
# 		 return(which(pairs_HP$sites_ID == HP[x,1] & pairs_HP$ID_P == HP[x,6]))
# 	}
	
# 	id_lhp <- foreach(x=1:nrow(HP),.combine=c) %dopar% {
# 		# Get the sites where the gall and the salix interact
# 		 return(which(pairs_HP$sites_ID == HP[x,1] & pairs_HP$ID_H == HP[x,5] & pairs_HP$ID_P == HP[x,6]))
# 	}
	
# 	# convert matches to 1
# 	XH2[id_xh2] = 1
# 	XP[id_xp] = 1
# 	LHP[id_lhp] = 1
	
# 	stopCluster(cl)
	
# } else{
	
	for(x in 1:nrow(HP)) {

		# Get the sites where the salix is present
		XH2[which(pairs_HP$sites_ID == HP[x,1] & pairs_HP$ID_H == HP[x,5])] = 1

		# Get the sites where the gall is present
		XP[which(pairs_HP$sites_ID == HP[x,1] & pairs_HP$ID_P == HP[x,6])] = 1
			
		# Get the sites where the gall and the salix interact
		LHP[which(pairs_HP$sites_ID == HP[x,1] & pairs_HP$ID_H == HP[x,5] & pairs_HP$ID_P == HP[x,6])] = 1

	}
#}


data_HP = data.frame(pairs_HP = pairs_HP, XH2 = XH2, XP = XP, LHP = LHP)
save(data_HP, file = "./data/HP.Rdata")

load("./data/HP.Rdata")
XH2 = data_HP$XH2
XP = data_HP$XP
LHP = data_HP$LHP

###################################################
# Combine the sets
names(pairs_SH) = c("sites_ID", "from", "to", "type")
names(pairs_HP) = c("sites_ID", "from", "to", "type")
pairs = rbind(pairs_SH, pairs_HP)
Xi = c(XS,XH2)
Xj = c(XH,XP)
Xij = Xi*Xj
Lij = c(LSH,LHP)

###################################################
# Get climatic data
library(dismo)
library(raster)
bclim = brick("./data/bioclim.grd")
salix_climate = extract(bclim, salix_int[,c('lon', 'lat')])
salix_climate = data.frame(salix_int, salix_climate)


save(salix_climate, file="./data/salix_climate.Rdata")

# Subset the data to have one line per site
sub_salix_climate = salix_climate[match(sites_ID,sites),c(1,7:25)]

# Drop NAs
sub_salix_climate = sub_salix_climate[is.na(sub_salix_climate$bio1)==FALSE,]

# Compute the PCA
PCA = princomp(na.omit(sub_salix_climate[,2:20]),cor = TRUE)
summary(PCA)

# Join the scores with the climate data
sub_salix_climate = cbind(sub_salix_climate,PCA$scores)

# Match the climate data to the full table
match_pairs_climate = match(as.character(pairs[,1]), 
as.character(sub_salix_climate[,1]))

expand_climate = sub_salix_climate[match_pairs_climate,]

# Final dataset
data = data.frame(pairs = pairs, Xi=Xi, Xj=Xj, Xij=Xij, Lij=Lij, climate = expand_climate)
data = subset(data,is.na(data$climate.bio1)==FALSE)
save(data,file = "./data/expand_data.Rdata")

