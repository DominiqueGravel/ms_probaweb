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
		probs = get_probs(modelC, modelL, newE)
		PXij = probs$PXij 
		PLij = probs$PLij 
		PLijXij = PXij*PLij 			
		PLijlow = probs$PLijlow
		PLijhigh = probs$PLijhigh	
	}

	return(list(PXij = PXij, PLij = PLij, PLijXij = PLijXij, PLijlow = PLijlow, PLijhigh = PLijhigh))	
} 