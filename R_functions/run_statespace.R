run_statespace <- function(lh, catch, index, lengthfreq, 
	obs_per_yr, model_name, years, dat_avail, 
	version="lb_statespace", model_dir, obs_meanlen, est_params=c("log_F_t_input", "log_q_I", "beta", "log_sigma_R", "S50", "CV_l")){

	RecType <- 0
	if(RecType!=2) RecDev_biasadj1 <- rep(1, length(years))
	if(RecType==2) RecDev_biasadj1 <- rep(0, length(years))


TmbList1 <- format_tmb_inputs(Nyears=length(years), Nlenbins=length(lh$mids),
	catch=catch, index=index, lengthfreq=lengthfreq,
	linf=lh$linf, vbk=lh$vbk, t0=lh$t0, 
	M=lh$M, AgeMax=lh$AgeMax, lbhighs=lh$highs,
	lbmids=lh$mids, Mat_a=lh$Mat_a, lwa=lh$lwa,
	lwb=lh$lwb, CV_catch=lh$CVcatch, CV_length=lh$CVlen,
	F1=lh$F1, SigmaR=lh$SigmaR, SigmaF=lh$SigmaF, 
	qcoef=lh$qcoef, R0=lh$R0, Sslope=lh$Sslope, 
	S50=lh$S50, dome=lh$dome, RecDev_biasadj=RecDev_biasadj1, 
	Fpen=1, Dpen=0, Dprior=c(0.8, 0.2), obs_per_yr=obs_per_yr, 
	RecType=RecType, h=lh$h, dat_avail=dat_avail, est_params=est_params)

saveRDS(TmbList1, file.path(model_dir, "Inputs1.rds"))

dyn.load(paste0(run_exe, "\\", dynlib(version)))
obj1 <- MakeADFun(data=TmbList1[["Data"]], parameters=TmbList1[["Parameters"]],
                  random=TmbList1[["Random"]], map=TmbList1[["Map"]], 
                  inner.control=list(maxit=1e3), hessian=FALSE)       


## Settings
obj1$env$inner.control <- c(obj1$env$inner.control, "step.tol"=c(1e-8,1e-12)[1], "tol10"=c(1e-6,1e-8)[1], "grad.tol"=c(1e-8,1e-12)[1]) 
obj1$hessian <- FALSE 
  Upr = rep(Inf, length(obj1$par))
    Upr[match("log_sigma_R",names(obj1$par))] = log(2)
    Upr[match("S50", names(obj1$par))] = dat_input$AgeMax
    Upr[match("Sslope", names(obj1$par))] = log(5)
    Upr[which(names(obj1$par)=="log_F_t_input")] = log(2)
    Upr[match("log_F_sd", names(obj1$par))] <- log(2)
    Lwr = rep(-Inf, length(obj1$par))
    Lwr[which(names(obj1$par)=="log_F_t_input")] = log(0.001) 
    Lwr[which(names(obj1$par)=="log_sigma_R")] = log(0.001)
    Lwr[which(names(obj1$par)=="log_F_sd")] <- log(0.001)     
    Lwr[which(names(obj1$par)=="dome")] <- 0

## Run optimizer
opt1 <- tryCatch( nlminb( start=obj1$par, objective=obj1$fn, gradient=obj1$gr, upper=Upr, lower=Lwr,
  control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10) ), error=function(e) NA)
if(all(is.na(opt1))==FALSE){
  opt1[["final_gradient"]] =  obj1$gr( opt1$par )      
  df1 <- data.frame(opt1$final_gradient, names(obj1$par), opt1$par)  
}    

## Standard errors
Report1 = obj1$report()  
Sdreport1 = tryCatch( sdreport(obj1), error=function(x) NA )
ParList <- obj1$env$parList( obj1$env$last.par.best )

if(all(is.na(Sdreport1))) RecDev_biasadj <- rep(1, scenario_setup$Nyears)
if(all(is.na(Sdreport1))==FALSE){
    SD <- summary(Sdreport1)
    if(RecType!=2) RecDev_biasadj <- 1 - SD[which(rownames(SD)=="Nu_input"), "Std. Error"]^2 / Report1$sigma_R^2    
    if(RecType==2) RecDev_biasadj <- rep(0, scenario_setup$Nyears)
}

TmbList <- format_tmb_inputs(Nyears=length(years), Nlenbins=length(lh$mids),
	catch=catch, index=index, lengthfreq=lengthfreq,
	linf=lh$linf, vbk=lh$vbk, t0=lh$t0, 
	M=lh$M, AgeMax=lh$AgeMax, lbhighs=lh$highs,
	lbmids=lh$mids, Mat_a=lh$Mat_a, lwa=lh$lwa,
	lwb=lh$lwb, CV_catch=lh$CVcatch, CV_length=lh$CVlen,
	F1=lh$F1, SigmaR=lh$SigmaR, SigmaF=lh$SigmaF, 
	qcoef=lh$qcoef, R0=lh$R0, Sslope=lh$Sslope, 
	S50=lh$S50, dome=lh$dome, RecDev_biasadj=RecDev_biasadj, 
	Fpen=1, Dpen=0, Dprior=c(0.8, 0.2), obs_per_yr=obs_per_yr, 
	RecType=RecType, h=lh$h, dat_avail=dat_avail, est_params=est_params)

saveRDS(TmbList, file.path(model_dir, "Inputs2.rds"))

dyn.load(paste0(run_exe, "\\", dynlib(version)))
obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList,
                  random=TmbList[["Random"]], map=TmbList[["Map"]], 
                  inner.control=list(maxit=1e3), hessian=FALSE) 

InitVal <- obj$fn(obj$par)

## Settings
obj$env$inner.control <- c(obj$env$inner.control, "step.tol"=c(1e-8,1e-12)[1], "tol10"=c(1e-6,1e-8)[1], "grad.tol"=c(1e-8,1e-12)[1]) 
obj$hessian <- FALSE 
  Upr = rep(Inf, length(obj$par))
    Upr[match("log_sigma_R",names(obj$par))] = log(2)
    Upr[match("S50", names(obj$par))] = dat_input$AgeMax
    Upr[match("Sslope", names(obj$par))] = log(5)
    Upr[which(names(obj$par)=="log_F_t_input")] = log(2)
    Upr[match("log_F_sd", names(obj$par))] <- log(2)
    Lwr = rep(-Inf, length(obj$par))
    Lwr[which(names(obj$par)=="log_F_t_input")] = log(0.001) 
    Lwr[which(names(obj$par)=="log_sigma_R")] = log(0.001)
    Lwr[which(names(obj$par)=="log_F_sd")] <- log(0.001)
    Lwr[which(names(obj$par)=="dome")] <- 0     

## Run optimizer
opt <- tryCatch( nlminb( start=obj$par, objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr,
  control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10) ), error=function(e) NA)

      for(i in 1:5){
      	if(all(is.na(opt))){
      		obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList,
                        random=TmbList[["Random"]], map=TmbList[["Map"]], 
                        inner.control=list(maxit=1e3), hessian=FALSE) 
            opt <- tryCatch( nlminb( start= obj$env$last.par.best[-obj$env$random] + rnorm(length(obj$par),0,0.2), 
              objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr, 
              control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10)), error=function(e) NA)
      	} else{
      		break
      	}
      }
      if(all(is.na(opt))){
      	write("NAs final gradient", file.path(model_dir, "NAs_final_gradient.txt"))
		    next
	    }
      if(all(is.na(opt))==FALSE){

      	opt[["final_gradient"]] = obj$gr( opt$par )       
      	df <- data.frame(opt$final_gradient, names(obj$par), opt$par)      

      	for(i in 1:5){
          if(abs(min(opt[["final_gradient"]]))>0.01){
            obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList,
                        random=TmbList[["Random"]], map=TmbList[["Map"]], 
                        inner.control=list(maxit=1e3), hessian=FALSE) 
            opt <- tryCatch( nlminb( start= obj$env$last.par.best[-obj$env$random] + rnorm(length(obj$par),0,0.2), 
              objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr, 
              control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10)), error=function(e) NA)
            if(all(is.na(opt))==FALSE) opt[["final_gradient"]] <- obj$gr(opt$par) 
            if(all(is.na(opt))) break
          } else{
            break
          } 
        }
        for(i in 1:5){
      	  if(all(is.na(opt))){
      		obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList,
                        random=TmbList[["Random"]], map=TmbList[["Map"]], 
                        inner.control=list(maxit=1e3), hessian=FALSE) 
            opt <- tryCatch( nlminb( start= obj$env$last.par.best[-obj$env$random] + rnorm(length(obj$par),0,0.2), 
              objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr, 
              control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10)), error=function(e) NA)
      	  } else{
      		break
      	  }
        }
        if(all(is.na(opt))){
      	  write("NAs final gradient", file.path(model_dir, "NAs_final_gradient.txt"))
		  next
	    } 
        if(abs(min(opt[["final_gradient"]]))>0.01) write(opt[["final_gradient"]], file.path(model_dir, "high_final_gradient.txt"))
      }      

## Standard errors
Report = tryCatch( obj$report(), error=function(x) NA)
  saveRDS(Report, file.path(model_dir, "Report.rds"))

Sdreport = tryCatch( sdreport(obj), error=function(x) NA )
  saveRDS(Sdreport, file.path(model_dir, "Sdreport.rds"))


  if(file.exists(file.path(fig_dir, paste0(model_name, "_output.png")))) unlink(file.path(fig_dir, paste0(model_name, "_output.png")), TRUE)
	png(file=file.path(fig_dir, paste0(model_name, "_output.png")), width=9, height=6, res=200, units="in")
      	par(mfrow=c(3,3), mar=c(3,3,2,0))
        FUN = function(InputMat, log=TRUE){
          if(log==TRUE) return(c( exp(InputMat[,1]-1*InputMat[,2]), rev(exp(InputMat[,1]+1*InputMat[,2]))))
          if(log==FALSE) return(c( InputMat[,1]-1*InputMat[,2], rev(InputMat[,1]+1*InputMat[,2])))
        } 
        ## Abundance
        Mat <- cbind("Year"=years, "Est"=Report$SB_t_hat)
        ymax <- ifelse(max(Mat[,c("Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("Est")], na.rm=TRUE)*1.5)
        matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), main="Biomass", lwd=2)
        if(all(is.na(Sdreport))==FALSE) if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lSB_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Index
        Mat <- cbind("Year"=years, "Est"=Report$N_t_hat)
        ymax <- ifelse(max(Mat[,c("Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("Est")], na.rm=TRUE)*1.5)
        matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), main="Abundance", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lN_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Recruitment
        Mat <- cbind("Year"=years, "Est"=Report$R_t_hat)
        ymax <- ifelse(max(Mat[,c("Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("Est")], na.rm=TRUE)*1.5)    
        matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), main="Recruitment", lwd=2)
        if(all(is.na(Sdreport))==FALSE) if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Fishing Mortality
        Mat <- cbind("Year"=years, "Est"=Report$F_t)
        ymax <- ifelse(max(Mat[,c("Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("Est")], na.rm=TRUE)*1.5)
        matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), main="Fishing Mortality", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Depletion
        Mat <- cbind("Year"=years, "Est"=Report$Depl)
        ymax <- ifelse(max(Mat[,c("Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("Est")], na.rm=TRUE)*1.5)
        matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), main="Depletion", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lD_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Catch
        tcatch <- catch_input
        pcatch <- rep(NA, length(years))
          names(pcatch) <- years
        pcatch[which(as.numeric(names(pcatch)) %in% years[as.numeric(names(catch_input))])] <- tcatch
        Mat <- cbind("Year"=years, "Obs"=pcatch, "Est"=Report$C_t_hat)
        ymax <- ifelse(max(Mat[,c("Obs", "Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("Obs", "Est")], na.rm=TRUE)*1.5)
        matplot(y=Mat[,c("Obs", "Est")], x=Mat[,c("Year")], type=c("p","l"), col=c("black", "red"), lty="solid", pch=19, ylim=c(0, ymax), main="Catch", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lC_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Index
        tindex <- cpue_input
        pindex <- rep(NA, length(years))
          names(pindex) <- years
        pindex[which(as.numeric(names(pindex)) %in% years[as.numeric(names(cpue_input))])] <- tindex
        Mat <- cbind("Year"=years, "Obs"=pindex, "Est"=Report$I_t_hat)
        ymax <- ifelse(max(Mat[,c("Obs", "Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("Obs", "Est")], na.rm=TRUE)*1.5)
        matplot(y=Mat[,c("Obs", "Est")], x=Mat[,c("Year")], type=c("p","l"), col=c("black", "red"), lty="solid", pch=19, ylim=c(0, ymax), main="Index", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lI_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Average Length
        Mat <- cbind("Year"=years, "Obs"=obs_meanlen, "Est"=Report$L_t_hat)
        ymax <- ifelse(max(Mat[,c("Obs", "Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("Obs", "Est")], na.rm=TRUE)*1.5)
        matplot(y=Mat[,c("Obs", "Est")], x=Mat[,c("Year")], type=c("p","l"), col=c("black", "red"), lty="solid", pch=19, ylim=c(0, ymax), main="Average Length", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="L_t_hat"),], log=FALSE), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
    dev.off()

        dyn.unload( paste0(run_exe,"\\", dynlib(version)) ) 

        Fmat_l <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),]
        Dmat_l <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lD_t"),]

        Fmat <- data.frame("Estimate"=exp(Fmat_l[,"Estimate"]), "SE"=Fmat_l[,"Std. Error"])
        Dmat <- data.frame("Estimate"=exp(Dmat_l[,"Estimate"]), "SE"=Dmat_l[,"Std. Error"])

        Outs <- NULL
        Outs$df <- df
        Outs$Ft <- Fmat
        Outs$Dt <- Dmat

        return(Outs)

}