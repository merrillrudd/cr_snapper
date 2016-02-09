plot_retro <- function(param, years, retro_dir, rec_type){

   	bh_dir <- ifelse(rec_type==0, "sharedR", ifelse(rec_type==1, "bhR", stop("no BH directory created")))
    model_dir <- file.path(retro_dir, bh_dir, years[length(years)])

    full_rep <- readRDS(file.path(model_dir, "Report.rds"))
    full_sd <- readRDS(file.path(model_dir, "Sdreport.rds"))

    if(param=="D"){
    	est1 <- full_rep$Depl
    	sd <- "lD_t"
    }
    if(param=="F"){
    	est1 <- full_rep$F_t
    	sd <- "lF_t"
    }
    if(param=="R"){
    	est1 <- full_rep$R_t
    	sd <- "lR_t"
    }	

    par(mfrow=c(1,1), mar=c(4,4,2,2))
 		FUN = function(InputMat, log=TRUE){
          if(log==TRUE) return(c( exp(InputMat[,1]-1*InputMat[,2]), rev(exp(InputMat[,1]+1*InputMat[,2]))))
          if(log==FALSE) return(c( InputMat[,1]-1*InputMat[,2], rev(InputMat[,1]+1*InputMat[,2])))
        }
        Mat <- cbind("Year"=years, "Est"=est1)
        ymax <- ifelse(max(Mat[,c("Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("Est")], na.rm=TRUE)*1.5)
        matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), lwd=2,
        	xlab="Year", ylab="Fishing Mortality", xaxs="i", yaxs="i")
        if(all(is.na(full_sd))==FALSE)if( !("condition" %in% names(attributes(full_sd)))) polygon( y=FUN(summary(full_sd)[which(rownames(summary(full_sd))==sd),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)

    ryears <- rev(years)
	for(i in 2:length(ryears)){

		rep <- readRDS(file.path(retro_dir, bh_dir, ryears[i], "Report.rds"))
		if(param=="D") est2 <- rep$Depl
		if(param=="F") est2 <- rep$F_t
		if(param=="R") est2 <- rep$R_t
		Mat <- cbind("Year"=years[1]:ryears[i], "Est"=est2)
		lines(y=Mat[,"Est"], x=Mat[,"Year"], col="#00000075", lty="solid", lwd=2)

	}

}