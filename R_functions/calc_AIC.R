calc_AIC <- function(model_names){

	aic_mat <- matrix(NA, nrow=length(model_names), ncol=4)
	colnames(aic_mat) <- c("AIC", "AICc", "deltaAIC", "deltaAICc")
	rownames(aic_mat) <- model_names

	for(i in 1:length(model_names)){
		dir <- file.path(init_dir, "output", model_names[i])
		input <- readRDS(file.path(dir, "Inputs2.rds"))
		report <- readRDS(file.path(dir, "Report.rds"))

		nll <- report$jnll
		params <- length(as.vector(unlist(input$Parameters))) - length(which(grepl("Nu_input", names(unlist(input$Parameters)))))
		# sampsize <- length(which(input$Data$I_t>0)) + length(which(input$Data$obs_per_yr>0)) + length(which(input$Data$C_t>0))
		sampsize <- length(which(input$Data$I_t>0)) + sum(input$Data$obs_per_yr) + length(which(input$Data$C_t>0))


		AIC <- 2*params + 2*nll
		AICc <- AIC + (2*params*(params + 1))/(sampsize - params - 1)

		aic_mat[i,1] <- AIC
		aic_mat[i,2] <- AICc
	}

	aic_mat[,"deltaAIC"] <- sapply(1:nrow(aic_mat), function(x) exp((min(aic_mat[,"AIC"]) -aic_mat[x,"AIC"])/2))
	aic_mat[,"deltaAICc"] <- sapply(1:nrow(aic_mat), function(x) exp((min(aic_mat[,"AICc"]) - aic_mat[x,"AICc"])/2))

	return(aic_mat)
}