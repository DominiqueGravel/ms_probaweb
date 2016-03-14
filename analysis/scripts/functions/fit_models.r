
# Compute the co-occurrence and the interaction models
fit_models = function(data, selection, funC, funL) {

	# MODEL CO-OCCURRENCE
	modelC = funC(data,selection)

	# MODEL INTERACTIONS
	modelL = funL(data,selection)

	return(list(modelC=modelC,modelL=modelL))	
}

#test_fit_models = fit_models(data,selection = TRUE, funC = C2, funL = L2)

# Wrapper around the fit_models function
fit_models.apply = function(data, funC, funL, selection) {
	
	if(sum(data$Xij)!=0 ) {
		models = fit_models(data,selection,funC,funL)
	 	modelC = models$modelC
	 	modelL = models$modelL
	}

	else {
		modelC = NULL
	 	modelL = NULL
	}

	return(list(i=as.character(data$IDi[1]),j=as.character(data$IDj[1]),modelC = modelC,modelL = modelL))
}

#fit_models.apply(c(1,1), funC=C2, funL=L2, selection=TRUE, DF_split)









