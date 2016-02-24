model_fits <- function(Report, Sdreport, years, meanlen_obs, index_obs, catch_obs, retro=FALSE){

	if(retro==FALSE) par(mfrow=c(3,3), mar=c(3,3,2,0))
	if(retro==TRUE) par(mfrow=c(3,3), mar=c(3,3,2,0), new=TRUE)
        FUN = function(InputMat, log=TRUE){
          if(log==TRUE) return(c( exp(InputMat[,1]-1*InputMat[,2]), rev(exp(InputMat[,1]+1*InputMat[,2]))))
          if(log==FALSE) return(c( InputMat[,1]-1*InputMat[,2], rev(InputMat[,1]+1*InputMat[,2])))
        } 
        ## Abundance
        Mat <- cbind("Year"=years, "Est"=Report$SB_t_hat)
        ymax <- 500
        matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), main="Biomass", lwd=2)
        if(all(is.na(Sdreport))==FALSE) if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lSB_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Abundance
        Mat <- cbind("Year"=years, "Est"=Report$N_t_hat)
        ymax <- 2.5
        matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), main="Abundance", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lN_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Recruitment
        Mat <- cbind("Year"=years, "Est"=Report$R_t_hat)
        ymax <- 1.5
        matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), main="Recruitment", lwd=2)
        if(all(is.na(Sdreport))==FALSE) if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Fishing Mortality
        Mat <- cbind("Year"=years, "Est"=Report$F_t)
        ymax <- 1
        matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), main="Fishing Mortality", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Depletion
        Mat <- cbind("Year"=years, "Est"=Report$Depl)
        ymax <- 1.1
        matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax), main="Relative Abundance", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lD_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Catch
        tcatch <- catch_obs
        pcatch <- rep(NA, length(years))
          names(pcatch) <- years
        pcatch[which(as.numeric(names(pcatch)) %in% years[as.numeric(names(catch_obs))])] <- tcatch
        Mat <- cbind("Year"=years, "Obs"=pcatch, "Est"=Report$C_t_hat)
        ymax <- 120
        matplot(y=Mat[,c("Obs", "Est")], x=Mat[,c("Year")], type=c("p","l"), col=c("black", "red"), lty="solid", pch=19, ylim=c(0, ymax), main="Catch", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lC_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Index
        tindex <- index_obs
        pindex <- rep(NA, length(years))
          names(pindex) <- years
        pindex[which(as.numeric(names(pindex)) %in% years[as.numeric(names(index_obs))])] <- tindex
        Mat <- cbind("Year"=years, "Obs"=pindex, "Est"=Report$I_t_hat)
        ymax <- 0.04
        matplot(y=Mat[,c("Obs", "Est")], x=Mat[,c("Year")], type=c("p","l"), col=c("black", "red"), lty="solid", pch=19, ylim=c(0, ymax), main="Index", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lI_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Average Length
        tmeanlen <- meanlen_obs
        pmeanlen <- rep(NA, length(years))
            names(pmeanlen) <- years
        pmeanlen[which(as.numeric(names(pmeanlen)) %in% years[as.numeric(names(meanlen_obs))])] <- tmeanlen
        Mat <- cbind("Year"=years, "Obs"=pmeanlen, "Est"=Report$L_t_hat)
        ymax <- 55
        matplot(y=Mat[,c("Obs", "Est")], x=Mat[,c("Year")], type=c("p","l"), col=c("black", "red"), lty="solid", pch=19, ylim=c(0, ymax), main="Average Length", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="L_t_hat"),], log=FALSE), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)

}