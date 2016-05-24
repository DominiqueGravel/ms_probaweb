##################################################################
# Script for Figure 1
# Illustrating the sampling of the metaweb of interactions among salix,
# gall insects and parasitoids across Europe
# Dominique Gravel
# October 29th, 2015
##################################################################

rm(list = ls())

# Load the data
load("./data/expand_data.Rdata")
load("./data/DF_split.Rdata")
load("./data/pairs.Rdata")

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
# Adjacency matrix for the metaweb
#########################################################

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
# Adjacency matrix for one site
#########################################

lw = data.frame(matrix(0, nr = Sall, nc = Sall))
names(lw) = unique_ID_all
row.names(lw) = unique_ID_all
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
# Convert the adjacency matrices to igraph
#########################################################

library(igraph)
g_mw = graph.adjacency(as.matrix(mw))
g_lw = graph.adjacency(as.matrix(lw))  

########################################################
# Set the layer attribute for every node
#########################################################

salix_ID = as.character(unique(pairs[pairs$pairs.type=="SH",1]))
gall_ID = as.character(unique(pairs[pairs$pairs.type=="HP",1]))
par_ID = as.character(unique(pairs[pairs$pairs.type=="HP",2]))

V(g_mw)$layer[V(g_mw)$name %in% salix_ID] = 1
V(g_mw)$layer[V(g_mw)$name %in% gall_ID] = 2
V(g_mw)$layer[V(g_mw)$name %in% par_ID] = 3

#########################################################
# Plot the network
#########################################################

dev.new(width = 6.5, height = 20)

# Function to set the layout
layout.k_partite <- function(g) {
  l <- layout.sugiyama(g)$layout[,2:1]
  l[,1] <- V(g)$layer
  l[,2] <- - l[,2] + 1 + max(l[,2])
  l
}

# Color code the nodes for the different layers & put in grey the ones that are absent
library(RColorBrewer)
pal <- brewer.pal(3, "Dark2")
ecol <- rgb(0.4, 0.4, 0.4, alpha=0.4)
wcol <- rgb(0.4, 0.4, 0.4)

dlw = apply(lw,1,sum)
V(g_mw)$size = 3
V(g_mw)$size[dlw!=0] = 5
V(g_mw)$color[V(g_mw)$name %in% salix_ID] = brewer.pal(5, "Dark2")[5]
V(g_mw)$color[V(g_mw)$name %in% gall_ID] = brewer.pal(5, "Dark2")[2]
V(g_mw)$color[V(g_mw)$name %in% par_ID] = brewer.pal(5, "Set1")[1]
V(g_mw)$color[dlw==0] = wcol

# Same for the edges
el_mw = apply(get.edgelist(g_mw), 1, paste, collapse = "-")
el_lw = apply(get.edgelist(g_lw), 1, paste, collapse = "-")
E(g_mw)$width = ifelse(el_mw %in% el_lw, 1, 0.3)
E(g_mw)$color = ifelse(el_mw %in% el_lw, "darkblue", ecol)

# Plot the network
plot(g_mw, layout = layout.k_partite(g_mw), vertex.label=NA, edge.arrow.mode=0, margin = c(-0.075,-0.075,-0.075,-0.075), asp=0, vertex.frame.color = NA)

# Save file
dev.copy2pdf(file = "./figures/mw_sampling.pdf")