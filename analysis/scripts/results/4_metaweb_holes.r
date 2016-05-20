##################################################################
# Script for Figure 5
# Illustrating the computation of interaction probabilities for a 
# pair of species
# Dominique Gravel
# October 29th, 2015
##################################################################

rm(list = ls())
load("./data/DF_split.Rdata")
load("./data/pairs.Rdata")

ID_S = unique(as.character(pairs$pairs.from[pairs$pairs.type=="SH"]))
ID_H = unique(as.character(pairs$pairs.to[pairs$pairs.type=="SH"]))
ID_P = unique(as.character(pairs$pairs.to[pairs$pairs.type=="HP"]))

S_S = length(ID_S)
S_H = length(ID_H)
S_P = length(ID_P)

#########################################################
# Compile the metaweb
#########################################################

L_SH = matrix(0, nr = S_S, nc = S_H)
X_SH = matrix(0,nr = S_S, nc = S_H)

L_HP = matrix(0, nr = S_H, nc = S_P)
X_HP = matrix(0,nr = S_H, nc = S_P)

for(n in 1:length(DF_split)) {

	# Compute the number of links
	nL = sum(DF_split[[n]]$Lij)

	# Compute the number of co-occurrences
	nX = sum(DF_split[[n]]$Xij)

	if(nX!=0) {

		# Indices for the different species
		if(pairs$pairs.type[n]=="SH") {
			i = which(ID_S == pairs$pairs.from[n])
			j = which(ID_H == pairs$pairs.to[n])			
			L_SH[i,j] = nL/nX
			X_SH[i,j] = nX
		}

		else {
			i = which(ID_H == pairs$pairs.from[n])
			j = which(ID_P == pairs$pairs.to[n])
			L_HP[i,j] = nL/nX
			X_HP[i,j] = nX			
		}
	}
	cat(n,'\n')
}

# Sort the matrices by degree
d_S = apply(L_SH,1,sum,na.rm=TRUE)
d_Hi = apply(L_SH,2,sum,na.rm=TRUE)
d_Hj = apply(L_HP,1,sum,na.rm=TRUE)
d_P = apply(L_HP,2,sum,na.rm=TRUE)

order_S = order(d_S, decreasing = FALSE)
order_Hi = order(d_Hi,decreasing = FALSE)
order_Hj = order(d_Hj, decreasing = FALSE)
order_P = order(d_P, decreasing = FALSE)

L_SH = L_SH[order_S,order_Hi]
L_HP = L_HP[order_Hj,order_P]
X_SH = X_SH[order_S,order_Hi]
X_HP = X_HP[order_Hj,order_P]

# Flip the matrices to have ennemies in colums and victims in rows
L_SH = t(L_SH)
L_HP = t(L_HP)
X_SH = t(X_SH)
X_HP = t(X_HP)

X_SH[X_SH>0] = 1 
X_HP[X_HP>0] = 1 

L_SH[L_SH>0] = 1
L_HP[L_HP>0] = 1
L_SH[X_SH == 0] = 2
L_HP[X_HP == 0] = 2

#########################################################
# Plot the results
#########################################################
dev.new(width = 8, height = 3)

par(mfrow = c(1,2))
par(mar = c(0.5,3,2,0.5))
#image(c(1:S_H),c(1:S_S),L_SH,col = rev(gray(seq(0,1,1/1000))),xlab = "", ylab = "",  cex.lab = 1.25,axes = FALSE)
image(c(1:S_H),c(1:S_S),L_SH,col = c("lightgrey", "black", "white"),xlab = "", ylab = "",  cex.lab = 1.25,axes = FALSE)
box()
mtext("Salix",side = 2, cex = 1.25, line = 0.5)
mtext("Galls",side = 3, cex = 1.25, line = 0.5)

par(mar = c(0.5,0.5,2,3))
#image(c(1:S_P),c(1:S_H),L_HP,col = rev(gray(seq(0,1,1/1000))),xlab = "", ylab = "",  cex.lab = 1.25,axes = FALSE)
image(c(1:S_P),c(1:S_H),L_HP,col = c("lightgrey", "black","white"),xlab = "", ylab = "",  cex.lab = 1.25,axes = FALSE)
box()
mtext("Galls",side = 3, cex = 1.25, line = 0.5)
mtext("Parasitoids",side = 4, cex = 1.25, line = 0.5)

dev.copy2pdf(file = "./figures/mw_holes.pdf")






