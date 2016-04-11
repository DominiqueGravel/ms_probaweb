##################################################################
# Script for Figure X,
# Illustrating the sampling of the metaweb of interactions among salix,
# gall insects and parasitoids across Europe
# Dominique Gravel
# October 29th, 2015
##################################################################

rm(list = ls())
setwd("/Users/DGravel/Documents/Manuscripts/Inprep/ms_probaweb")

# Load the data
load("analysis/data/expand_data.Rdata")
load("analysis/data/DF_split.Rdata")
load("analysis/data/pairs.Rdata")

IDi = as.character(data$pairs.from)
IDj = as.character(data$pairs.to)
Si = length(unique(IDi))
Sj = length(unique(IDj))
unique_IDi = unique(IDi)
unique_IDj = unique(IDj)
IDcomm = IDj[match(unique_IDi,unique_IDj,nomatch=0)]
unique_ID_all = unique(c(IDi,IDj))
Sall = length(unique_ID_all)
DF = data.frame(sites = data$pairs.sites_ID, IDi = IDi, IDj = IDj, Xi=data$Xi, Xj=data$Xj, Xij = data$Xij, Lij = data$Lij)

#########################################################
# Network sampling
#########################################################

# Adjacency matrix for the metaweb
# Loop around all pairs of species to test if there are interactions
# Put everything in a square matrix
mw = data.frame(matrix(0, nr = Sall, nc = Sall))
names(mw) = unique_ID_all
row.names(mw) = unique_ID_all

for(n in 1:nrow(pairs)) {

	# Compute the number of links
	nL = sum(DF_split[[n]]$Lij)

	if(nL!=0) {

		# Get the victim
		i = which(unique_ID_all == DF_split[[n]]$IDi[1])

		# Get the ennemy index
		j = which(unique_ID_all == DF_split[[n]]$IDj[1])

		# Put the record in the mw
		mw[i,j] = 1
		mw[j,i] = 1
	}
}



#########################################
# Compute the local network for one site
lw = matrix(0, nr = Sall, nc = Sall)
DF_split_sites = split(DF,DF$sites)
site_index = 369 # Pick one that has a decent number of links

for(n in 1:nrow(pairs)) {

	# Compute the number of links
	nL = DF_split_sites[[site_index]]$Lij[n]

	if(nL!=0) {

		# Get the victim
		i = which(unique_ID_all == pairs$pairs.from[n])

		# Get the ennemy index
		j = which(unique_ID_all == pairs$pairs.to[n])

		# Put the record in the mw
		lw[i,j] = 1
		lw[j,i] = 1
	}
}

#########################################################
# Convert the metaweb to igraph
#########################################################

library(igraph)
g = graph.incidence(mw,add.names=NULL)

#########################################################
# Set the layer attribute for every node
#########################################################

salix_ID = as.character(unique(pairs[pairs$pairs.type=="SH",1]))
gall_ID = as.character(unique(pairs[pairs$pairs.type=="HP",1]))
par_ID = as.character(unique(pairs[pairs$pairs.type=="HP",2]))

V(g)$layer[V(g)$name %in% salix_ID] = 1
V(g)$layer[V(g)$name %in% gall_ID] = 2
V(g)$layer[V(g)$name %in% par_ID] = 3

#########################################################
# Plot the metaweb
#########################################################

quartz(width = 4, height = 9)

# Function to set the layout
layout.k_partite <- function(g) {
  l <- layout.sugiyama(g)$layout[,2:1]
  l[,1] <- V(g)$layer
  l[,2] <- - l[,2] + 1 + max(l[,2])
  l
}

# Color code the nodes that are present
dlw = apply(lw,1,sum)
V(g)$color[V(g)$name %in% salix_ID] = "darkgreen"
V(g)$color[V(g)$name %in% gall_ID] = "darkred"
V(g)$color[V(g)$name %in% par_ID] = "darkblue"
V(g)$color[dlw==0] = "white"

V(g)$size = 2
V(g)$size[dlw!=0] = 5

# Plot the network

plot(g, layout = layout.k_partite(g), vertex.label=NA, edge.arrow.mode=0, edge.width=0.3, margin = c(-0.075,-0.075,-0.075,-0.075), asp=0)

dev.copy2pdf(file = "ms/figures/metaweb_sampling.pdf")





