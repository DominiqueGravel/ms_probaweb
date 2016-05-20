# Compute the co-occurrence and the interaction models
fit_models = function(data, funC, funL, selection) {
	
	if(sum(data$Xij)!=0 ) {

		# MODEL CO-OCCURRENCE
		modelC = funC(data,selection)

		# MODEL INTERACTIONS
		modelL = funL(data,selection)
	}

	else {
		modelC = NULL
	 	modelL = NULL
	}

	return(list(i=as.character(data$IDi[1]),j=as.character(data$IDj[1]), modelC = modelC,modelL = modelL))
}