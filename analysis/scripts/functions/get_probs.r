# Compute co-occurrence and interaction probabilities for a given environment
get_probs = function(modelC, modelL, newE) {

	# PREDICT CO-OCCURRENCE PROBABILITY
	if(is.null(modelC$ij)) 
		PXij = predict(modelC$i, type = "response", newdata = newE)*predict(modelC$j, type = "response", newdata = newE)	
	else 
		PXij = predict(modelC$ij, type = "response", newdata = newE)	

	# PREDICT INTERACTION PROBABILITY
#	PLij = predict(modelL, se.fit = TRUE, newdata = newE)
#	PLijlow = exp(PLij$fit-1.96*PLij$se.fit)/(1+exp(PLij$fit-1.96*PLij$se.fit))
#	PLijhigh = exp(PLij$fit+1.96*PLij$se.fit)/(1+exp(PLij$fit+1.96*PLij$se.fit))
#	PLijhigh[PLijhigh=="NaN"]=1
#	PLij = exp(PLij$fit)/(1+exp(PLij$fit))

	PLij = predict(modelL, type = "response", newdata = newE)
	PLij = exp(PLij)/(1+exp(PLij))
	PLij[is.na(PLij)]=1

#	return(list(PLijXij=PLij*PXij,PLij=PLij,PXij=PXij, PLijlow = PLijlow, PLijhigh = PLijhigh))	
	return(list(PLijXij=PLij*PXij,PLij=PLij,PXij=PXij))	

}

#test_get_probs = get_probs(modelC=test_models[[1]], modelL=test_models[[2]], newE=data$E)

# Wrapper around the get_probs function
get_probs.apply = function(model_list,newE) {

	modelC = model_list$modelC
	modelL = model_list$modelL

	if(is.null(modelC)) {
		PXij = NA
		PLij = NA
		PLijXij = NA
		PLijlow = NA
		PLijhigh = NA
	}
	else {
		probs = get_probs(modelC, modelL, newE)
		PXij = probs$PXij 
		PLij = probs$PLij 
		PLijXij = PXij*PLij 			
		PLijlow = probs$PLijlow
		PLijhigh = probs$PLijhigh	
	}
	return(list(PXij = PXij, PLij = PLij, PLijXij = PLijXij, PLijlow = PLijlow, PLijhigh = PLijhigh))	
} 






