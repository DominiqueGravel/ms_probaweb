##################################################################
# Script for Figure X,
# Illustrating the computation of interaction probabilities for a 
# pair of species
# Dominique Gravel
# October 29th, 2015
##################################################################

rm(list = ls())
setwd("/Users/DGravel/Documents/Manuscripts/Inprep/ms_probaweb")
load("analysis/data/DF_split.Rdata")
load("analysis/data/expand_data.Rdata")

IDi = as.character(data$pairs.IDi)
IDj = as.character(data$pairs.IDj)
Si = length(unique(IDi))
Sj = length(unique(IDj))
unique_IDi = unique(IDi)
unique_IDj = unique(IDj)
IDcomm = IDj[match(unique_IDi,unique_IDj,nomatch=0)]
unique_ID_all = unique(c(IDi,IDj))
Sall = length(unique_ID_all)
DF = data.frame(sites = data$pairs.sites_ID, IDi = IDi, IDj = IDj, Xi=data$Xi, Xj=data$Xj, Xij = data$Xij, Lij = data$Lij)

#########################################################
# Compile the metaweb
#########################################################

Lij= matrix(0, nr = Sall, nc = Sall)
Xij = matrix(0, nr = Sall, nc = Sall)

n = 1
for(i in 1:Si) {
	for(j in 1:Sj) {

		# Compute the number of links
		nL = sum(DF_split[[n]]$Lij)

		# Compute the number of co-occurrences
		nX = sum(DF_split[[n]]$Xij)

		if(nX!=0) {

			# Get the victim
			index_i = which(unique_ID_all == as.character(DF_split[[n]]$IDi)[1])

			# Get the ennemy index
			index_j = which(unique_ID_all == as.character(DF_split[[n]]$IDj)[1])

			# Put the record in the mw
			Lij[index_i,index_j] = nL
			Xij[index_i,index_j] = nX
		}
		n = n+1
	}
}


# Subset the two matrices to have an Si X Sj matrix
PL = Lij/Xij
subPL = PL[match(unique_IDi,unique_ID_all),match(unique_IDj,unique_ID_all)]
NAs = is.na(subPL)
subPL[NAs]=0

# Sort the matrix by generality
di = apply(subPL,1,sum,na.rm=TRUE)
dj = apply(subPL,2,sum,na.rm=TRUE)
subPL = subPL[order(di,decreasing = FALSE), order(dj, decreasing = TRUE)]

# Flip the matrices to have ennemies in colums and victims in rows
NAs = NAs[order(di,decreasing = FALSE), order(dj, decreasing = TRUE)]
subPL = t(subPL)
NAs = t(NAs)

#########################################################
# Plot the results
#########################################################
quartz(width = 8, height = 4)

par(mfrow = c(1,2),mar = c(2,3,2,0))
image(c(1:Sj),c(1:Si),subPL,col = rev(gray(seq(0,1,1/1000))),xlab = "", ylab = "",  cex.lab = 1.25,axes = FALSE)
box()
mtext("Ennemy",side = 3, cex = 1.25, line = 0.5)
mtext("Victim",side = 2, cex = 1.25, line = 0.5)

par(mar = c(2,0,2,3))
image(c(1:Sj),c(1:Si),NAs,col = c("white","black"),xlab = "", ylab = "", cex.lab = 1.25,axes = FALSE)
box()
mtext("Ennemy",side = 3, cex = 1.25, line = 0.5)
mtext("Victim",side = 4, cex = 1.25, line = 0.5)

dev.copy2pdf(file = "ms/figures/mw_holes.pdf")






