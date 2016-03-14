# Collect the expected network from a the list of probabilities
collect_net = function(Si,Sj,probs_list) {

	npairs = Si*Sj
	L = matrix(nr = Si, nc = Sj)
	X = matrix(nr = Si, nc = Sj)
	LX = matrix(nr = Si, nc = Sj)
	L.CI = matrix(nr = Si, nc = Sj)  

	for(x in 1:npairs) {
		i = probs_list[[x]]$i
		j = probs_list[[x]]$j
		L[i,j] = probs_list[[x]]$PLij
		X[i,j] = probs_list[[x]]$PXij
		LX[i,j] = probs_list[[x]]$PLijXij
		L.CI[i,j] = probs_list[[x]]$PLijhigh - probs_list[[x]]$PLijlow
	}
	return(list(L = L, X = X, LX = LX, L.CI = L.CI))
}

# Collect the total likelihood from a list of sum log-likelihoods
collect_LL = function(LL_list, DF_split, pairs) {
	npairs = length(DF_split)
	sumLL = numeric(npairs)
	npars = numeric(npairs)
	AIC = numeric(npairs)
	type = numeric(npairs)

	for(x in 1:npairs) {
		data = DF_split[[x]]
		IDs = as.matrix(data[1,1:2])
		test =  which(pairs[,1]==as.numeric(IDs[1]) & pairs[,2]==as.numeric(IDs[2]) | pairs[,1]==as.numeric(IDs[2]) & pairs[,2]==as.numeric(IDs[1]))
		if(length(test)!=0)	type[x] = as.character(pairs[test,3])

		# Compute the likelihood if the species are not the same		
		if(as.numeric(IDs[1]) != as.numeric(IDs[2])) {

			# Test if the species are co-occurring
			if(sum(data$Xij)!=0)
				if(sum(data$Lij!=0)) {
					sumLL[x] = LL_list[[x]]$sumLL
					npars[x] = LL_list[[x]]$npars
			}
			else {

			}
		}
		cat(x,'\n')
				
	}

	sumLL[sumLL==-Inf] = 0
	return(cbind(type = type, sumLL = sumLL, npars = npars, AIC = -2*sumLL + 2*npars))
}

