# Compute co-occurrence and interaction probabilities for a given environment
get_probs = function(model_list, newE) {

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
		if(is.null(modelC$ij)) 
			PXij = predict(modelC$i, type = "response", newdata = newE)*predict(modelC$j, type = "response", newdata = newE)	
		else 
			PXij = predict(modelC$ij, type = "response", newdata = newE)	
		PLij = predict(modelL, type = "response", newdata = newE)
		PLijXij = PXij*PLij 			
	}

	return(list(PXij = PXij, PLij = PLij, PLijXij = PLijXij))	
} 