format_tmb_inputs <- function(Nyears, Nlenbins, catch, index, lengthfreq, meanlen,
	linf, vbk, t0, M, AgeMax, lbhighs, lbmids, Mat_a, lwa, lwb, 
	CV_catch, CV_length, F1, SigmaR, qcoef, R0, Sslope, S50, dome,
	RecDev_biasadj, Fpen, Dpen, Dprior, obs_per_yr, SigmaF,
	RecType, h, dat_avail, est_params){

	if("lengthfreq" %in% dat_avail){
		if(is.matrix(lengthfreq)){
          n_lc <- nrow(lengthfreq)
          LC_yrs <- as.numeric(rownames(lengthfreq))
          LF <- as.matrix(lengthfreq)
    	}
    	if(is.vector(lengthfreq)){
    	  n_lc <- 1
    	  LC_yrs <- Nyears
    	  LF <- t(as.matrix(lengthfreq))
    	}
		n_ml <- 0
		ML_yrs <- as.vector(0)
		ML_t <- as.vector(0)
	}
	if("catch_total" %in% dat_avail){
		n_c <- length(catch)
		C_yrs <- as.numeric(names(catch))
		rel_c <- 0
		C_t <- catch
	}
	if("index_total" %in% dat_avail){
		n_i <- length(index)
		I_yrs <- as.numeric(names(index))
		rel_i <- 0
		I_t <- index
	}
	if(any(grepl("index", dat_avail))==FALSE){
		n_i <- 0
		I_yrs <- as.vector(0)
		rel_i <- 0
		I_t <- as.vector(0)
	}
	if(any(grepl("catch", dat_avail))==FALSE){
		n_c <- 0
		C_yrs <- as.vector(0)
		rel_c <- 0
		C_t <- as.vector(0)
	}
	if(any(grepl("lengthfreq", dat_avail))==FALSE){
		n_lc <- 0
		LC_yrs <- as.vector(0)
		LF <- as.matrix(0)
	}
	if("meanlength" %in% dat_avail){
		n_ml <- length(meanlen)
		ML_yrs <- as.numeric(names(meanlen))
		ML_t <- meanlen
	}


	Data <- list(n_t=Nyears, n_lb=Nlenbins,
   		n_c=n_c, n_i=n_i,
   		n_lc=n_lc, n_ml=n_ml,
   		T_yrs=1:Nyears, C_yrs=C_yrs,
   		I_yrs=I_yrs,
   		LC_yrs=LC_yrs,
   		ML_yrs=ML_yrs, 
   		rel_c=rel_c, rel_i=rel_i,
   		obs_per_yr=obs_per_yr, 
   		RecType=RecType, I_t=I_t,
   		C_t=C_t, ML_t=ML_t,
   		LF=LF, linf=linf,
   		vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax,
   		lbhighs=lbhighs, lbmids=lbmids,
   		Mat_a=Mat_a, lwa=lwa, lwb=lwb,
      Fpen=Fpen, Dpen=Dpen, Dprior=Dprior,
   		RecDev_biasadj=RecDev_biasadj)


	Parameters <- list(log_F_sd=log(SigmaF), 
		log_F_t_input=log(rep(F1,Nyears)),
		log_q_I=log(qcoef),
		beta=log(R0),
		log_sigma_R=log(SigmaR), 
		Sslope=Sslope, S50=S50, dome=dome,
    	CV_c=CV_catch, CV_l=CV_length,
		Nu_input=rep(0,Nyears))

    if(RecType!=2) Random <- c("Nu_input")
    if(RecType==2) Random <- NULL

	Map <- list()

  all_params <- c("log_F_sd", "log_F_t_input", "log_q_I", 
    "beta", "log_sigma_R", "Sslope", "S50", "dome", "CV_c", 
    "CV_l")

  for(a in 1:length(all_params)){
    if(all_params[a] %in% est_params) next
    Map[[all_params[a]]] <- NA
    Map[[all_params[a]]] <- factor(Map[[all_params[a]]])
  }

    if(RecType==2){
        Map[["Nu_input"]] <- 1:length(Parameters[["Nu_input"]])
        Map[["Nu_input"]] <- rep(NA, length(Parameters[["Nu_input"]]))
        Map[["Nu_input"]] <- factor(Map[["Nu_input"]])

        Map[["log_F_t_input"]] <- 1:length(Parameters[["log_F_t_input"]])
        Map[["log_F_t_input"]][2:length(Parameters[["log_F_t_input"]])] <- NA
        Map[["log_F_t_input"]] <- factor(Map[["log_F_t_input"]])
    }

    if(length(Map)==0) Map <- NULL

    Return <- list("Parameters"=Parameters, 
    	"Data"=Data, "Random"=Random, 
    	"Map"=Map)
    return(Return)


}