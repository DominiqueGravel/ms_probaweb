##################################################################
# Script splitting the dataset in a list of tables for each pair of
# species
# Dominique Gravel
# October 29th, 2015
##################################################################

rm(list = ls())
# Load the data
load('analysis/data/expand_data.Rdata')

# Sub select environmental variables (annual average temperature and annual preciptiations)
DF = data.frame(IDi = data$pairs.from, IDj = data$pairs.to, Xi=data$Xi, Xj=data$Xj, Xij = data$Xij, Lij = data$Lij)

T = data$climate.bio1
T2 = T^2
PP = data$climate.bio12
PP2 = PP^2

E = data.frame(T = T, T2 = T2, PP = PP, PP2 = PP2)
DF$E = E
DF$sites = data$info_pairs.sites_ID

# Split the DF into a list of pairs of species
pairs_ID = as.factor(paste(DF$IDi,DF$IDj))
DF_split = split(DF,pairs_ID)

pairs = unique(data[,c(2,3,4)])

save(DF_split,file = "analysis/data/DF_split.Rdata")
save(pairs,file = "analysis/data/pairs.Rdata")