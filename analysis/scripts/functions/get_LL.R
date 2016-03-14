# Compute the likelihood for a set of observations, given the model
get_LL = function(data, modelC, modelL) {

	probs = get_probs(newE = data$E,modelC,modelL)

	# NUMBER OF PARAMETERS
	if(is.null(modelC$ij)) nparsXij = length(modelC$i$coefficients)*2	
	else nparsXij = length(modelC$ij$coefficients)
	nparsLij = length(modelL$coefficients)	

	# COMPUTATION OF THE LIKELIHOOD
	PLijXij = probs$PLij*probs$PXij 
	LL = numeric(length(data$Xij))
	LL[data$Lij == 1] = log(PLijXij[data$Lij==1])
	LL[data$Lij == 0] = log(1 - PLijXij[data$Lij==0])		

	return(list(LL = LL, npars = (nparsXij+nparsLij)))	
}
#test_LL = get_LL(data,modelC=test_models[[1]], modelL=test_models[[2]])


# Wrapper around the get_LL function
get_LL.apply = function(model_list,data) {

	modelC = model_list$modelC
	modelL = model_list$modelL

	if(is.null(modelC)) {
		sumLL = 0
		npars = 0
	}
	else {
		fit = get_LL(data,modelC, modelL)
		sumLL = sum(fit$LL)
		npars = fit$npars
	}

	return(c(sumLL,npars))

}
