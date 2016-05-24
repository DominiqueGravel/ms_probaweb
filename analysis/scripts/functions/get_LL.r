# Compute the likelihood for a set of observations, given the model
get_LL = function(model_list, data) {

	modelC = model_list$modelC
	modelL = model_list$modelL

	if(is.null(modelC)) {
		sumLL = 0
		npars = 0
	}

	else {
		# NUMBER OF PARAMETERS
		if(is.null(modelC$ij)) nparsXij = length(modelC$i$coefficients)*2	
		else nparsXij = length(modelC$ij$coefficients)
		nparsLij = length(modelL$coefficients)	
		npars = (nparsXij+nparsLij)

		# COMPUTATION OF THE PROBABILITIES
		probs = get_probs(model_list, newE = data$E)
		PLijXij = probs$PLij*probs$PXij 

		# COMPUTATION OF THE LIKELIHOOD
		LL = numeric(length(data$Xij))
		LL[data$Lij == 1] = log(PLijXij[data$Lij==1])
		LL[data$Lij == 0] = log(1 - PLijXij[data$Lij==0])	
		sumLL = sum(LL)
	}

	return(c(sumLL, npars))

}
