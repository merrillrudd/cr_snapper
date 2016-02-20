plot_retro <- function(param, years, retro_dir, rec_type, uncertainty=FALSE, dat_avail){

   	bh_dir <- ifelse(rec_type==0, "sharedR", ifelse(rec_type==1, "bhR", stop("no BH directory created")))
    label <- ifelse(param=="D", "Relative Abundance", ifelse(param=="F", "Fishing Mortality", ifelse(param=="R", "Recruitment", stop("parameter not programmed"))))

    y <- NULL
    for(i in 1:length(dat_avail)) y <- paste(y, dat_avail[i], sep="_")
    dat_dir <- file.path(retro_dir, bh_dir, y)

    model_dir <- file.path(dat_dir, years[length(years)])

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
        ymax <- ifelse(param=="D", 1.5, ifelse(param=="F", 3, ifelse(param=="R", 2, "parameter not coded")))
        matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), lwd=2,
        	xlab="Year", ylab=label, xaxs="i", yaxs="i")
        if(all(is.na(full_sd))==FALSE)if( !("condition" %in% names(attributes(full_sd)))) polygon( y=FUN(summary(full_sd)[which(rownames(summary(full_sd))==sd),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)

    ryears <- rev(years)
	for(i in 2:length(ryears)){

		rep <- readRDS(file.path(dat_dir, ryears[i], "Report.rds"))
		sdrep <- readRDS(file.path(dat_dir, ryears[i], "Sdreport.rds"))

		if(param=="D") est2 <- rep$Depl
		if(param=="F") est2 <- rep$F_t
		if(param=="R") est2 <- rep$R_t
		Mat <- cbind("Year"=years[1]:ryears[i], "Est"=est2)
		lines(y=Mat[,"Est"], x=Mat[,"Year"], col="#00000075", lty="solid", lwd=2)
        if(all(is.na(sdrep))==FALSE & uncertainty==TRUE)if( !("condition" %in% names(attributes(sdrep)))) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))==sd),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(0,0,0,alpha=0.2), border=NA)

	}

}