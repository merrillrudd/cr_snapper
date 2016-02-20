run_lprof <- function(index, lengthfreq, catch, obs_per_yr,
	years, dat_avail=c("lengthfreq", "index_total"), meanlen,
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	rec_type, adjust_param=FALSE, adjust_val=FALSE, lprof_dir, plot=FALSE,
	lh, run){

	path <- file.path(lprof_dir, adjust_param)
	dir.create(path, showWarnings=FALSE)

	if(run==FALSE & file.exists(file.path(path, adjust_val[1], "Report.rds"))==FALSE){
		message("no lprof results - running lprof")
		run <- TRUE
	}

	nll_vec <- rep(NA, length(adjust_val))    


	if(run==TRUE){

    	for(i in 1:length(adjust_val)){    

    		dat_input <- create_inputs(param=adjust_param, val=adjust_val[i], lh_dat=lh)    

    		model_dir <- file.path(path, adjust_val[i])
    		dir.create(model_dir, showWarnings=FALSE)    

    		run_model <- run_statespace(lh=dat_input, years=years, catch=catch, 
    		  	index=index, lengthfreq=lengthfreq, obs_per_yr=obs_per_yr, 
    		  	dat_avail=dat_avail, est_params=est_params, RecType=rec_type,
    		  	model_name=NULL, model_dir=model_dir, obs_meanlen=meanlen, plot=plot)      

    		if(file.exists(file.path(model_dir, "Report.rds"))) rep <- readRDS(file.path(model_dir, "Report.rds"))
    		if(file.exists(file.path(model_dir, "Report.rds"))==FALSE) next    

    		nll_vec[i] <- rep$jnll    

    		rm(run_model)    

    	}		
	}

	if(run==FALSE & file.exists(file.path(path, adjust_val[1], "Report.rds"))){
		for(i in 1:length(adjust_val)){
			rep <- readRDS(file.path(path, adjust_val[i], "Report.rds"))
			nll_vec[i] <- rep$jnll

			rm(rep)
		}
	}


	return(nll_vec)

}