run_retro <- function(index, lengthfreq, catch, obs_per_yr,
	years, dat_avail=c("lengthfreq", "index_total"), obs_meanlen,
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	rec_type=0, adjust_param=FALSE, adjust_val=FALSE){

  dat_input <- create_inputs(simdir="retrospective", param=adjust_param, val=adjust_val, lh_dat=cr_lh)   

  tyrs <- 1:length(years)

  for(i in 1:length(tyrs)){

    model_dir <- file.path(retro_dir, years[i])
    if(file.exists(model_dir)) unlink(model_dir, TRUE)
    dir.create(model_dir, showWarnings=FALSE)   

    	iyrs <- tyrs[1]:tyrs[i]
    	if(is.null(catch)==FALSE) icatch <- catch[which(names(catch) %in% iyrs)]
    	if(is.null(catch)) icatch <- NULL
    	if(length(iyrs) > 2) iindex <- index[which(names(index) %in% iyrs)]
    	if(length(iyrs)<=2){
    		iindex <- NULL
    		dat_avail2 <- c("lengthfreq")
    		est_params2 <- c("log_F_t_input", "S50", "CV_l")
    	}
    	if(length(iyrs) > 2){
    		dat_avail2 <- dat_avail
    		est_params2 <- est_params
    	}
    	ilf <- lengthfreq[which(rownames(lengthfreq) %in% iyrs),]
    	iobs <- obs_per_yr[which(tyrs %in% iyrs)]  

    	run_model <- run_statespace(lh=dat_input, years=iyrs, catch=icatch, 
    		index=iindex, lengthfreq=ilf, obs_per_yr=iobs,
    		dat_avail=dat_avail2, est_params=est_params2, RecType=rec_type,
    		model_name=model_name, model_dir=model_dir, obs_meanlen=ml$all_gears,
    		plot=FALSE)  

    	saveRDS(run_model$df, file.path(model_dir, "df.rds"))

  }


}