get_output <- function(model_name, dat_avail, rec_type, est_params, SPR_cut=0.3,
	adjust_param=FALSE, adjust_val=FALSE, years, index, lengthfreq, catch, meanlen, obs_per_yr){

  dat_input <- create_inputs(simdir=model_name, param=adjust_param, val=adjust_val, lh_dat=cr_lh)  

  model_dir <- file.path(init_dir, "output", model_name)
  if(file.exists(model_dir)) unlink(model_dir, TRUE)
  dir.create(model_dir, showWarnings=FALSE)  

  run_model <- run_statespace(lh=dat_input, years=years, catch=catch, 
  	index=index, lengthfreq=lengthfreq, obs_per_yr=obs_per_yr, 
  	dat_avail=dat_avail, est_params=est_params, RecType=rec_type,
  	model_name=model_name, model_dir=model_dir, obs_meanlen=meanlen)  

  report <- readRDS(file.path(init_dir, "output", model_name, "Report.rds"))  

  Fref <- uniroot(calc_ref, lower=0, upper=50, Mat_a=report$Mat_a, W_a=report$W_a, 
    M=report$M, S_a=report$S_a, ref="Fref", cut=0.3)$root  

  SPR <- calc_ref(Mat_a=report$Mat_a, M=report$M, W_a=report$W_a,
  	S_a=report$S_a, ref="SPR",
  	F=mean(run_model$Ft[(nrow(run_model$Ft)-2):nrow(run_model$Ft),"Estimate"]))  

  FFref <- run_model$Ft[nrow(run_model$Ft), "Estimate"]/Fref 

  Outs <- NULL
  Outs$df <- run_model$df
  Outs$FFref <- FFref
  Outs$SPR <- SPR
  Outs$report <- report
  Outs$Ft <- run_model$Ft
  Outs$Dt <- run_model$Dt
  return(Outs) 

}