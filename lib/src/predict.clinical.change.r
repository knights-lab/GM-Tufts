# Use the microbiome to predict the change in a clinical covariate (stratified by Treatment Group)
predict.clinical.change <- function(map_delta, otu, clinical.vars, n.iterations=100, otu.delta=TRUE)
{
	if(otu.delta == FALSE)
	{
		outputfile="initial-mb-predict-delta-covariate.txt"
		cat("Use the initial microbiome to predict change in clinical covariate, stratified by treatment group\n", file=outputfile)
	} else {
		outputfile="delta-mb-predict-delta-covariate.txt"
		cat("Use the change in microbiome to predict change in clinical covariate, stratified by treatment group\n", file=outputfile)	
	}
	
	samples_A <- rownames(map_delta)[map_delta$GROUP=="A"]
	samples_B <- rownames(map_delta)[map_delta$GROUP=="B"]

#	participants_A <- as.factor(map_delta[samples_A,"PARTICIPANTID"])
#	participants_B <- as.factor(map_delta[samples_B,"PARTICIPANTID"])

	for(i in 1:length(clinical.vars))
	{	
		clinical.var <- clinical.vars[i]
		delta.clinical.var <- paste("delta", clinical.var,sep="_")
#		delta.clinical.var <- clinical.var
		
		r_squared_A <- numeric()
		r_squared_B <- numeric()
		for(a in 1:n.iterations)
		{
				y <- map_delta[samples_A, delta.clinical.var,drop=F]
				samples_A_filtered <- rownames(y)[!is.na(y)]
				rf.resultA <- rf.cross.validation.classification(x=otu[samples_A_filtered,], y=map_delta[samples_A_filtered, delta.clinical.var], 
																								nfolds=-1, regression=T)  
				r_squared_A[a] <- rf.resultA$r_squared				
				
				y <- map_delta[samples_B, delta.clinical.var,drop=F]
				samples_B_filtered <- rownames(y)[!is.na(y)]
				rf.resultB <- rf.cross.validation.classification(x=otu[samples_B_filtered,], y=map_delta[samples_B_filtered, delta.clinical.var], 
																								nfolds=-1, regression=T)  
				r_squared_B[a] <- rf.resultB$r_squared				
		}
		
		cat(clinical.var, "Group A R^2:", mean(r_squared_A), "\n", file=outputfile, append=T)
		cat(clinical.var, "Group B R^2:", mean(r_squared_B), "\n", file=outputfile, append=T)
	}
}
