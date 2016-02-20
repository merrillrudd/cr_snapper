rm(list=ls())

library(RColorBrewer)
#########################
## notes
#########################

## may need monthly time-steps (monthly recruitment)
## should compare to LB-SPR and other methods
## alternate size at first maturity: 
	## 34 cm (Rojas 1996 Gulf of Nicoya)
	## 23.5 cm (Colombia snappers)
	## 2/3 linf

#########################
## directories 
#########################

init_dir <- "C:\\Git_Projects\\cr_snapper"

output_dir <- file.path(init_dir, "output")
dir.create(output_dir, showWarnings=FALSE)

retro_dir <- file.path(init_dir, "retrospective")
dir.create(retro_dir, showWarnings=FALSE)

lprof_dir <- file.path(init_dir, "likelihood_profiles")
dir.create(lprof_dir, showWarnings=FALSE)

fig_dir <- file.path(init_dir, "figs")
dir.create(fig_dir, showWarnings=FALSE)

model_fits_dir <- file.path(fig_dir, "model_fits")
dir.create(model_fits_dir, showWarnings=FALSE)


setwd(init_dir)
source("R_functions\\functions.R")

#########################
## setup life history
#########################

cr_lh <- create_lh(species="CRSNAP", selex="asymptotic")

#########################
## read in data + 
## initial adjustments
#########################

data_raw <- read.csv(file.path(init_dir, "cr_snapper_database.csv"), 
	stringsAsFactors=FALSE)

data_raw$Year <- sapply(1:nrow(data_raw), function(x) adjust_date(data_raw$Date[x], out="year"))
data_raw$Month <- sapply(1:nrow(data_raw), function(x) adjust_date(data_raw$Date[x], out="month"))
data_raw$TL_cm <- as.numeric(data_raw$TL_cm)
data_raw$W_g <- as.numeric(data_raw$W_g)
	data_raw$W_g[which(data_raw$W_g==0)] <- NA
data_raw$W_kg <- data_raw$W_g/1000

#############################
## subset Lutjanus guttatus
#############################

lg <- data_raw[which(data_raw$Species=="Lutjanus guttatus" & data_raw$Year>1950),]
write.csv(lg, file.path(init_dir, "cr_snapper_filtered.csv"))

#############################
## GSI
#############################

gsi <- read.csv(file.path(init_dir, "GSI_2013_females.csv"), stringsAsFactors=FALSE)
gsi_yrs <- unique(gsi$Year)[order(unique(gsi$Year))]
avg_gsi <- sapply(1:length(gsi_yrs), function(x) mean(gsi$GSI[which(gsi$Year==gsi_yrs[x])]))
names(avg_gsi) <- gsi_yrs

plot(x=gsi$Twt, y=gsi$GSI)

#############################
## subset by gears
#############################
setwd(init_dir)
source("R_functions\\functions.R")

lg_bl <- lg[which(lg$Gear=="Bottom Longline"),]
lg_g <- lg[which(lg$Gear=="Gillnet"),]

## annual mean length
png(file=file.path(fig_dir, "mean_length_catch.png"), width=6, height=5, units="in", res=200)
ml <- mean_length(data=lg, plot=TRUE)
dev.off()


###############################
## length frequency
###############################

setwd(init_dir)
source("R_functions\\functions.R")

png(file=file.path(fig_dir, "length_frequency_catch.png"), width=9, height=7, units="in", res=200)
lf <- length_frequency(binwidth=1, linf=cr_lh$linf,
 lmat=cr_lh$Lmat, data=lg, plot=TRUE, weight=TRUE)
dev.off()


#################################
## catch and effort
#################################
setwd(init_dir)
source("R_functions\\functions.R")

fishery_data <- catch_effort(data=lg, sep_index=TRUE)
catch <- fishery_data$catch
cpue_bl <- fishery_data$cpue_bl
cpue_g <- fishery_data$cpue_g
years <- c((fishery_data$years[1]-10):(fishery_data$years[1]-1),fishery_data$years)
years_obs <- fishery_data$years
obs_high <- c(rep(0, 10), fishery_data$obs_high)
obs_low <- c(rep(0, 10), fishery_data$obs_low)


par(mfrow=c(2,1))
plot(cpue_bl, pch=19, ylim=c(0, max(cpue_bl, na.rm=TRUE)))
plot(cpue_g, pch=19, ylim=c(0, max(cpue_g, na.rm=TRUE)))

#######################
## data setup
#######################

setwd(init_dir)
source("R_functions\\functions.R")

## bottom longline cpue index
cpue_input <- cpue_bl[which(is.na(cpue_bl)==FALSE)] ## choose bottom longline cpue only
cpue_yrs <- which(years %in% names(cpue_bl)[which(is.na(cpue_bl)==FALSE)])
names(cpue_input) <- cpue_yrs

## length frequency
lf_input <- lf[which(rowSums(lf)>0),]
lf_yrs <- which(years %in% rownames(lf)[which(rowSums(lf)>0)])
rownames(lf_input) <- lf_yrs

## catch - approximation based on 300 metric tons/year
# catch_input <- 3e8/mean(as.numeric(lg$W_g), na.rm=TRUE)
# catch_yrs <- years[length(years)]
# names(catch_input) <- catch_yrs

## mean length
meanlen_input <- ml$all_gears
meanlen_yrs <- which(years %in% years_obs)
names(meanlen_input) <- meanlen_yrs

########################################
## length-based state-space assessment
########################################

library(TMB)
run_exe <- file.path(init_dir, "inst", "executables")
version <- "lb_statespace"
setwd(run_exe)
dyn.unload(paste0(version, ".dll"))
# file.remove(paste(version, c(".dll", ".o"), sep=""))
# compile(paste0(version, ".cpp"))

### base

setwd(init_dir)
source("R_functions\\functions.R")

base <- get_output(model_name="base", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	lh=cr_lh)

## low obs per year

lowObs <- get_output(model_name="lowObs", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low, 
	lh=cr_lh)

## dome + low obs per year

dome <- get_output(model_name="dome", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low, 
	lh=cr_lh, adjust_param="dome", adjust_val=0.01)

## dome + high obs per year

dome_highobs <- get_output(model_name="dome_highobs", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high, 
	lh=cr_lh, adjust_param="dome", adjust_val=0.01)

### base without beverton-holt stock recruit 

base_nbh <- get_output(model_name="base_nbh", dat_avail=c("lengthfreq", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	lh=cr_lh)

### low obs per year without beverton-holt stock recruit 

nbh_lowObs <- get_output(model_name="nbh_lowObs", dat_avail=c("lengthfreq", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	lh=cr_lh)

########################################
## compare BH and observation sizes
########################################

names <- c("lowObs", "nbh_lowObs", "dome")
aic_bhobs <- calc_AIC(names)

### continue forward with low ESS for length frequency
### beverton-holt stock-recruit curve
### asymptotic selectivity
### likely more accurately represents uncertainty


################################################
## likelihood profiles for important parameters
################################################

setwd(init_dir)
source("R_functions\\functions.R")

## steepness

h_vec <- seq(0.3, 1, by=0.05)
h_prof <- run_lprof(index=cpue_input, lengthfreq=lf_input, catch=NULL, obs_per_yr=obs_low,
	years=years, dat_avail=c("lengthfreq", "index_total"), meanlen=meanlen_input,
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	rec_type=1, adjust_param="h", adjust_val=h_vec, lprof_dir=lprof_dir, lh=cr_lh, run=FALSE)
par(mfrow=c(1,1))
plot(x=h_vec, y=h_prof, pch=19)
h_vec[which(h_prof==min(h_prof))] ## 0.6


## natural mortality

M_vec <- seq(0.05, 1, by=0.05)
M_prof <- run_lprof(index=cpue_input, lengthfreq=lf_input, catch=NULL, obs_per_yr=obs_low,
	years=years, dat_avail=c("lengthfreq", "index_total"), meanlen=meanlen_input,
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	rec_type=1, adjust_param="M", adjust_val=M_vec, lprof_dir=lprof_dir, lh=cr_lh2, run=FALSE)
par(mfrow=c(1,1))
plot(x=M_vec, y=M_prof, pch=19)
M_vec[which(M_prof==min(M_prof, na.rm=TRUE))]

## linf 
linf_vec <- seq(58, 70, by=1)
linf_prof <- run_lprof(index=cpue_input, lengthfreq=lf_input, catch=NULL, obs_per_yr=obs_low,
	years=years, dat_avail=c("lengthfreq", "index_total"), meanlen=meanlen_input,
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	rec_type=1, adjust_param="linf", adjust_val=linf_vec, lprof_dir=lprof_dir, lh=cr_lh2, run=FALSE)
par(mfrow=c(1,1))
plot(x=linf_vec, y=linf_prof, pch=19)
linf_vec[which(linf_prof==min(linf_prof, na.rm=TRUE))]

## vbk 
vbk_vec <- seq(0.15, 0.45, by=0.02)
vbk_prof <- run_lprof(index=cpue_input, lengthfreq=lf_input, catch=NULL, obs_per_yr=obs_low,
	years=years, dat_avail=c("lengthfreq", "index_total"), meanlen=meanlen_input,
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	rec_type=1, adjust_param="vbk", adjust_val=vbk_vec, lprof_dir=lprof_dir, lh=cr_lh2, run=TRUE)
par(mfrow=c(1,1))
plot(x=vbk_vec, y=vbk_prof, pch=19)
vbk_vec[which(vbk_prof==min(vbk_prof, na.rm=TRUE))]

## lmat 
lmat_vec <- seq(28, 42, by=1)
lmat_prof <- run_lprof(index=cpue_input, lengthfreq=lf_input, catch=NULL, obs_per_yr=obs_low,
	years=years, dat_avail=c("lengthfreq", "index_total"), meanlen=meanlen_input,
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	rec_type=1, adjust_param="Lmat", adjust_val=lmat_vec, lprof_dir=lprof_dir, lh=cr_lh2, run=TRUE)
par(mfrow=c(1,1))
plot(x=lmat_vec, y=lmat_prof, pch=19)
lmat_vec[which(lmat_prof==min(lmat_prof, na.rm=TRUE))]


png(file.path(fig_dir, "likelihood_profiles.png"), width=8, height=9, res=200, units="in")
par(mfrow=c(2,2))
plot(x=h_vec, y=h_prof, pch=19, xlab="Natural Mortality (M)", ylab="NLL",
	ylim=c(min(h_prof)*1.2, min(h_prof)*0.2))
abline(v=cr_lh2$h, lty=2)
text(x=1.0, y=min(h_prof)*0.25, "A", font=2, cex=1.2)
plot(x=M_vec, y=M_prof, pch=19, xlab="Steepness (h)", ylab="NLL",
	ylim=c(min(M_prof,na.rm=TRUE)*1.2, min(M_prof, na.rm=TRUE)*0.2))
abline(v=cr_lh2$M, lty=2)
text(x=1.0, y=min(M_prof,  na.rm=TRUE)*0.25, "B", font=2, cex=1.2)
plot(x=linf_vec, y=linf_prof, pch=19, xlab="Linfinity", ylab="NLL",
	ylim=c(min(linf_prof, na.rm=TRUE)*1.2, min(linf_prof, na.rm=TRUE)*0.2))
abline(v=cr_lh2$linf, lty=2)
text(x=70, y=min(linf_prof,  na.rm=TRUE)*0.25, "C", font=2, cex=1.2)
plot(x=vbk_vec, y=vbk_prof, pch=19, xlab="von Bertalanffy k", ylab="NLL",
	ylim=c(min(vbk_prof, na.rm=TRUE)*1.2, min(vbk_prof, na.rm=TRUE)*0.2))
abline(v=cr_lh2$vbk, lty=2)
text(x=0.45, y=min(vbk_prof,  na.rm=TRUE)*0.25, "D", font=2, cex=1.2)
mtext(side=2, "Negative Log Likelihood", outer=TRUE)

dev.off()

## rerun base case with lowest ll for steepness

cr_lh2 <- create_inputs(param="h", val=0.6, lh_dat=cr_lh)

setwd(init_dir)
source("R_functions\\functions.R")

base2 <- get_output(model_name="base2", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	lh=cr_lh2)

singleFR <- get_output(model_name="singleFR", dat_avail=c("meanlength", "index_total"),
	rec_type=2, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=NULL, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	lh=cr_lh2)


####################################################
## model uncertainty - mimic mean length estimator
####################################################


########################################
## choose base model
########################################

choose <- base2
choose_rep <- readRDS(file.path(output_dir, "base2", "Report.rds"))
choose_sd <- readRDS(file.path(output_dir, "base2", "Sdreport.rds"))

FUN = function(InputMat, log=TRUE){

  if(log==TRUE) return(c( exp(InputMat[,1]-1*InputMat[,2]), rev(exp(InputMat[,1]+1*InputMat[,2]))))

  if(log==FALSE) return(c( InputMat[,1]-1*InputMat[,2], rev(InputMat[,1]+1*InputMat[,2])))

} 

png(file.path(fig_dir, "base_DF.png"), width=8, height=10, res=200, units="in")
par(mfrow=c(2,1), mar=c(0,0,0,0), omi=c(1,1,1,1))
## Fishing Mortality
Mat <- cbind("Year"=years, "Est"=choose_rep$F_t)
ymax <- 1
matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax),  lwd=2, xlab="Year", xaxs="i", yaxs="i", xaxt="n", yaxt="n")
if(all(is.na(choose_sd))==FALSE)if( !("condition" %in% names(attributes(choose_sd)))) polygon( y=FUN(summary(choose_sd)[which(rownames(summary(choose_sd))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
axis(2, at=seq(0, ymax, by=0.2), las=2)
mtext(side=2, "Fishing Mortality", cex=1.2, line=3)

abline(h=F30, col="blue", lwd=2)
abline(h=F40, col="blue", lwd=2, lty=2)

text(x=1998, y=0.9, "A", font=2, cex=1.2)


## relative abundance
Mat <- cbind("Year"=years, "Est"=choose_rep$Depl)
ymax <- 1
matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax),  lwd=2, xlab="Year", xaxs="i", yaxs="i", yaxt="n")
if(all(is.na(choose_sd))==FALSE)if( !("condition" %in% names(attributes(choose_sd)))) polygon( y=FUN(summary(choose_sd)[which(rownames(summary(choose_sd))=="lD_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
axis(2, at=seq(0, ymax, by=0.2), las=2)
mtext(side=2, "Relative Abundance", cex=1.2, line=3)

abline(h=0.3, col="blue", lwd=2)
abline(h=0.4, col="blue", lwd=2, lty=2)

text(x=1998, y=0.9, "B", font=2, cex=1.2)

dev.off()





png(file.path(fig_dir, "base_F.png"), width=8, height=7, res=200, units="in")
par(mfrow=c(1,1))
## Fishing Mortality
Mat <- cbind("Year"=years, "Est"=choose_rep$F_t)
ymax <- 1
matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax),  lwd=2, xlab="Year", ylab="Estimated Fishing Mortality")
if(all(is.na(choose_sd))==FALSE)if( !("condition" %in% names(attributes(choose_sd)))) polygon( y=FUN(summary(choose_sd)[which(rownames(summary(choose_sd))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
abline(h=F30, col="blue", lwd=2)
abline(h=F40, col="blue", lwd=2, lty=2)
dev.off()

png(file.path(fig_dir, "base_D.png"), width=8, height=7, res=200, units="in")
par(mfrow=c(1,1))
## relative abundance
Mat <- cbind("Year"=years, "Est"=choose_rep$Depl)
ymax <- 1
matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax),  lwd=2, xlab="Year", ylab="Relative Abundance")
if(all(is.na(choose_sd))==FALSE)if( !("condition" %in% names(attributes(choose_sd)))) polygon( y=FUN(summary(choose_sd)[which(rownames(summary(choose_sd))=="lD_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
abline(h=0.3, col="blue", lwd=2)
abline(h=0.4, col="blue", lwd=2, lty=2)
dev.off()

########################################
## importance of data types
########################################

## remove index

rmIndex <- get_output(model_name="rmIndex", dat_avail=c("lengthfreq"),
	rec_type=1, est_params=c("log_F_t_input", "log_sigma_R", "S50", "CV_l"),
	years=years, index=NULL, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	lh=cr_lh2)

rmIndex_rep <- readRDS(file.path(output_dir, "rmIndex", "Report.rds"))
rmIndex_sd <- readRDS(file.path(output_dir, "rmIndex", "Sdreport.rds"))

## remove length comp, replace with mean length

subML <- get_output(model_name="subML", dat_avail=c("meanlength", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_sigma_R", "log_q_I", "S50", "CV_c"),
	years=years, index=cpue_input, lengthfreq=NULL, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	lh=cr_lh2)

subML_rep <- readRDS(file.path(output_dir, "subML", "Report.rds"))
subML_sd <- readRDS(file.path(output_dir, "subML", "Sdreport.rds"))




png(file.path(fig_dir, "compare_rmIndex_DF.png"), width=8, height=10, res=200, units="in")
par(mfrow=c(2,1), mar=c(0,0,0,0), omi=c(1,1,1,1))
## fishing mortality
Mat <- cbind("Year"=years, "Est"=choose_rep$F_t)
ymax <- 1
matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax),  lwd=2, xlab="Year", ylab="Fishing Mortality", xaxt="n", yaxt="n", xaxs="i", yaxs="i")
if(all(is.na(choose_sd))==FALSE)if( !("condition" %in% names(attributes(choose_sd)))) polygon( y=FUN(summary(choose_sd)[which(rownames(summary(choose_sd))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
axis(2, at=seq(0, ymax, by=0.2), las=2)
mtext(side=2, "Fishing Mortality", cex=1.2, line=3)

Mat <- cbind("Year"=years, "Est"=rmIndex_rep$F_t)
lines(y=Mat[,"Est"], x=Mat[,"Year"], col="black", lwd=2)
if(all(is.na(rmIndex_sd))==FALSE)if( !("condition" %in% names(attributes(rmIndex_sd)))) polygon( y=FUN(summary(rmIndex_sd)[which(rownames(summary(rmIndex_sd))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(0,0,0,alpha=0.2), border=NA)
abline(h=0.3, col="blue", lwd=2)
abline(h=0.4, col="blue", lwd=2, lty=2)
text(x=1998, y=0.9, "A", font=2, cex=1.2)

## relative abundance
Mat <- cbind("Year"=years, "Est"=choose_rep$Depl)
ymax <- 1
matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax),  lwd=2, xlab="Year", yaxt="n", xaxs="i", yaxs="i")
if(all(is.na(choose_sd))==FALSE)if( !("condition" %in% names(attributes(choose_sd)))) polygon( y=FUN(summary(choose_sd)[which(rownames(summary(choose_sd))=="lD_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
axis(2, at=seq(0, ymax, by=0.2), las=2)
mtext(side=2, "Relative Abundance", cex=1.2, line=3)

Mat <- cbind("Year"=years, "Est"=rmIndex_rep$Depl)
lines(y=Mat[,"Est"], x=Mat[,"Year"], col="black", lwd=2)
if(all(is.na(rmIndex_sd))==FALSE)if( !("condition" %in% names(attributes(rmIndex_sd)))) polygon( y=FUN(summary(rmIndex_sd)[which(rownames(summary(rmIndex_sd))=="lD_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(0,0,0,alpha=0.2), border=NA)
abline(h=0.3, col="blue", lwd=2)
abline(h=0.4, col="blue", lwd=2, lty=2)
text(x=1998, y=0.9, "B", font=2, cex=1.2)
mtext(side=1, "Year", cex=1.2, line=3)

dev.off()

png(file.path(fig_dir, "compare_subML_DF.png"), width=8, height=10, res=200, units="in")
par(mfrow=c(2,1), mar=c(0,0,0,0), omi=c(1,1,1,1))
## fishing mortality
Mat <- cbind("Year"=years, "Est"=choose_rep$F_t)
ymax <- 1
matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax),  lwd=2, xlab="Year", ylab="Fishing Mortality", yaxt="n", xaxt="n", xaxs="i", yaxs="i")
if(all(is.na(choose_sd))==FALSE)if( !("condition" %in% names(attributes(choose_sd)))) polygon( y=FUN(summary(choose_sd)[which(rownames(summary(choose_sd))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
axis(2, at=seq(0, ymax, by=0.2), las=2)
mtext(side=2, "Fishing Mortality", cex=1.2, line=3)

Mat <- cbind("Year"=years, "Est"=subML_rep$F_t)
lines(y=Mat[,"Est"], x=Mat[,"Year"], col="black", lwd=2)
if(all(is.na(subML_sd))==FALSE)if( !("condition" %in% names(attributes(subML_sd)))) polygon( y=FUN(summary(subML_sd)[which(rownames(summary(subML_sd))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(0,0,0,alpha=0.2), border=NA)
abline(h=0.3, col="blue", lwd=2)
abline(h=0.4, col="blue", lwd=2, lty=2)
text(x=1998, y=0.9, "A", font=2, cex=1.2)

## relative abundance
Mat <- cbind("Year"=years, "Est"=choose_rep$Depl)
ymax <- 1
matplot(y=Mat[,c("Est")], x=Mat[,c("Year")], type="l", col=c("red"), lty="solid", ylim=c(0, ymax),  lwd=2, xlab="Year", yaxt="n", xaxs="i", yaxs="i")
if(all(is.na(choose_sd))==FALSE)if( !("condition" %in% names(attributes(choose_sd)))) polygon( y=FUN(summary(choose_sd)[which(rownames(summary(choose_sd))=="lD_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
axis(2, at=seq(0, ymax, by=0.2), las=2)
mtext(side=2, "Relative Abundance", cex=1.2, line=3)

Mat <- cbind("Year"=years, "Est"=subML_rep$Depl)
lines(y=Mat[,"Est"], x=Mat[,"Year"], col="black", lwd=2)
if(all(is.na(subML_sd))==FALSE)if( !("condition" %in% names(attributes(subML_sd)))) polygon( y=FUN(summary(subML_sd)[which(rownames(summary(subML_sd))=="lD_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(0,0,0,alpha=0.2), border=NA)
abline(h=0.3, col="blue", lwd=2)
abline(h=0.4, col="blue", lwd=2, lty=2)
text(x=1998, y=0.9, "B", font=2, cex=1.2)
mtext(side=1, "Year", cex=1.2, line=3)

dev.off()


########################################
## sensitivities
########################################

## sensitivity - adjust F1 low

F1low <- get_output(model_name="F1low", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	adjust_param="F1", adjust_val=0.01,
	lh=cr_lh2)

## sensitivity - adjust F1 low - other value

F1low2 <- get_output(model_name="F1low2", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	adjust_param="F1", adjust_val=0.2,
	lh=cr_lh2)

## sensitivity - adjust F1 low - other value

F1high <- get_output(model_name="F1high", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	adjust_param="F1", adjust_val=0.45,
	lh=cr_lh2)


########################################
## retrospective
########################################

setwd(init_dir)
source("R_functions\\functions.R")

## all data
run_retro(index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	years=years, dat_avail=c("lengthfreq", "index_total"), 
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	rec_type=1, retro_dir=retro_dir, lh=cr_lh2)

## length comp only
run_retro(index=NULL, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	years=years, dat_avail=c("lengthfreq"), 
	est_params=c("log_F_t_input","log_sigma_R", "S50", "CV_l"),
	rec_type=1, retro_dir=retro_dir, lh=cr_lh2)

## cpue and mean length
run_retro(index=cpue_input, lengthfreq=NULL, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	years=years, dat_avail=c("meanlength", "index_total"), 
	est_params=c("log_F_t_input", "log_sigma_R", "log_q_I", "S50", "CV_c"),
	rec_type=1, retro_dir=retro_dir, lh=cr_lh2)


png(file.path(fig_dir, "retro_depl_lengthindex.png"), width=10, height=8, res=200, units="in")
plot_retro(param="D", years=years, retro_dir=retro_dir, rec_type=1, uncertainty=TRUE, dat_avail=c("lengthfreq", "index_total"))
dev.off()

png(file.path(fig_dir, "retro_F_lengthindex.png"), width=10, height=8, res=200, units="in")
plot_retro(param="F", years=years, retro_dir=retro_dir, rec_type=1, uncertainty=TRUE, dat_avail=c("lengthfreq", "index_total"))
dev.off()

plot_retro(param="R", years=years, retro_dir=retro_dir, rec_type=1, uncertainty=TRUE, dat_avail=c("lengthfreq", "index_total"))

png(file.path(fig_dir, "retro_depl_lengthonly.png"), width=10, height=8, res=200, units="in")
plot_retro(param="D", years=years, retro_dir=retro_dir, rec_type=1, uncertainty=TRUE, dat_avail=c("lengthfreq"))
dev.off()

plot_retro(param="F", years=years, retro_dir=retro_dir, rec_type=1, uncertainty=TRUE, dat_avail=c("lengthfreq"))
plot_retro(param="R", years=years, retro_dir=retro_dir, rec_type=1, uncertainty=TRUE, dat_avail=c("lengthfreq"))


########################################
## per recruit curves
########################################

F_vec <- seq(0,2, by=0.01)
SBPR <- YPR <- rep(NA, length(F_vec))

base_rep <- readRDS(file.path(init_dir, "output", "lowObs", "Report.rds"))
for(i in 1:length(F_vec)){
	SBPR[i] <- calc_ref(Mat_a=cr_lh2$Mat_a, W_a=cr_lh2$W_a, M=cr_lh2$M, S_a=base_rep$S_a, F=F_vec[i], ref="SBPR")
	YPR[i] <- calc_ref(Mat_a=cr_lh2$Mat_a, W_a=cr_lh2$W_a, M=cr_lh2$M, S_a=base_rep$S_a, F=F_vec[i], ref="YPR")
}	


## what is impact of changing vulnerability through regulations (e.g. size limit?)
png(file.path(fig_dir, "SPR_YPR.png"), width=6, height=5, res=200, units="in")
par(mfrow=c(1,1), mar=c(5,5,5,5))
plot(x=F_vec, y=SBPR, type="o", ylim=c(0, max(SBPR)), xaxs="i", yaxs="i", pch=19, lwd=2,
	xlab="Fishing Mortality", ylab="Spawning Biomass Per Recruit", yaxt="n", font.lab=2, xpd=NA)
axis(2, at=pretty(c(0, max(SBPR))), las=2)

F30 <- F_vec[which(round(SBPR)==round(0.30*SBPR[1]))][1]
S30 <- 0.30*SBPR[1]

F40 <- F_vec[which(round(SBPR)==round(0.4*SBPR[1]))][1]
S40 <- 0.4*SBPR[1]

F20 <- F_vec[which(round(SBPR)==round(0.2*SBPR[1]))][1]
S20 <- 0.2*SBPR[1]


# segments(x0=0, x1=F40, y0=S40, y1=S40, lty=2)
abline(h=S30, col="red", lty=2)
abline(h=S40, col="blue", lty=2)

segments(x0=F30, x1=F30, y0=0, y1=S30, lty=2, col="red")
segments(x0=F40, x1=F40, y0=0, y1=S40, lty=2, col="blue")

text(paste0("F30%\n", F30), x=F30+0.1, y=25, col="red", font=2)
text(paste0("F40%\n", F40), x=F40-0.1, y=25, col="blue", font=2)

par(new=TRUE)
plot(x=F_vec, y=YPR, type="o", col=gray(0.5), ylim=c(0, max(YPR)), xaxs="i", yaxs="i",
	pch=19, lwd=2, xaxt="n", yaxt="n", ylab="", xlab="")
axis(4, at=pretty(c(0, max(YPR))), las=2)
mtext("Yield Per Recruit", side=4, line=3, col=gray(0.5), font=2)
dev.off()

0.4*SBPR[1]
F_vec[which(round(SBPR)==round(0.4*SBPR[1]))]

########################################
## kobe plot
########################################

png(file.path(fig_dir, "base_kobe.png"), width=8, height=7, res=200, units="in")
plot(lowObs$Dt[,"Estimate"], lowObs$Ft[,"Estimate"]/F30, xlim=c(0,1.1), ylim=c(0, 1.1), 
	pch=19, col="blue", xaxs="i", yaxs="i", xpd=NA,
	xlab="Relative Abundance", ylab=expression("F/F"[F30]), type="o")
points(lowObs$Dt[length(years), "Estimate"], lowObs$Ft[length(years), "Estimate"]/F30, cex=2, col="blue")
text(lowObs$Dt[11:(length(years)-1), "Estimate"]+0.04, lowObs$Ft[11:(length(years)-1), "Estimate"]/F30+0.03, labels=years[11:(length(years)-1)])
text(lowObs$Dt[nrow(lowObs$Dt), "Estimate"], lowObs$Ft[nrow(lowObs$Ft), "Estimate"]/F30+0.03, labels="2015")
abline(h=1, lty=2)
abline(v=0.3, lty=2)
polygon(x=c(0,1.1,1.1,0), y=c(1,1,1.1,1.1), border=NA, 
	col="#AA000050")
polygon(x=c(0,0.3,0.3,0), y=c(0,0,1.1,1.1), border=NA,
	col="#AA000050")
dev.off()



########################################
## lb-spr
########################################

setwd(init_dir)
source("R_functions\\functions.R")

## stock and fleet objects
StockPars <- NULL
StockPars$MK <- cr_lh2$M/cr_lh2$vbk
StockPars$Linf <- cr_lh2$linf
StockPars$CVLinf <- cr_lh2$CVlen
StockPars$MaxSD <- 2
StockPars$L50 <- cr_lh2$Lmat
StockPars$L95 <- cr_lh2$linf*(1-exp(-cr_lh2$vbk*(7-cr_lh2$t0)))
StockPars$FecB  <- 3
StockPars$Walpha <- cr_lh2$lwa
StockPars$Wbeta <- cr_lh2$lwb
StockPars$Steepness <- cr_lh2$h
StockPars$Mpow <- 0
StockPars$AgeMax <- cr_lh2$AgeMax

SizeBins <- NULL
SizeBins$Linc <- cr_lh2$binwidth
SizeBins$ToSize <- max(cr_lh2$highs)

FleetPars <- NULL
FleetPars$SL50 <- round(cr_lh2$linf*(1-exp(-cr_lh2$vbk*(cr_lh2$S50 - cr_lh2$t0))))
FleetPars$SL95 <- round(cr_lh2$linf*(1-exp(-cr_lh2$vbk*(cr_lh2$S50+1 - cr_lh2$t0))))
FleetPars$FM <- 0.8

LBSim <- LBSPRSim(StockPars, FleetPars, SizeBins)
# LBSim$SPR

## expected
plot(LBSim$LenMids, LBSim$LCatchFished, lwd=2, type="l")
LBLenDat <- LBSim$LCatchFished*100

## observed

lbspr_out <- lapply(1:nrow(lf), function(x) DoOpt(StockPars, 
	LenDat=lf[x,], SizeBins=SizeBins, mod="LBSPR"))

lbspr_spr <- sapply(1:nrow(lf), function(x) lbspr_out[[x]]$Ests["SPR"])
lbspr_fm <- sapply(1:nrow(lf), function(x) lbspr_out[[x]]$Ests["FM"])


### kobe plot
# col <- brewer.pal(9, "Blues")
plot(lbspr_spr, lbspr_fm/0.87, xlim=c(0,1.1), ylim=c(0,25), 
	pch=19, col="blue", xaxs="i", yaxs="i", xpd=NA,
	xlab="SPR", ylab="F/M = 0.87", type="o")
points(lbspr_spr[length(years)], lbspr_fm[length(years)]/0.87, cex=2, col="blue")
text(lbspr_spr-0.04, lbspr_fm/0.87+0.3, labels=years)
# text(lbspr_spr[(length(years))]+0.01, lbspr_fm[length(years)]+0.04, labels="2014, 2015")
abline(h=1, lty=2)
abline(v=0.3, lty=2)
polygon(x=c(0,17,17,0), y=c(1,1,25,25), border=NA, 
	col="#AA000050")
polygon(x=c(0,0.3,0.3,0), y=c(0,0,25,25), border=NA,
	col="#AA000050")



## compare overfished status
png(file.path(fig_dir, "compare_models_overfishing.png"), width=8, height=7, res=200, units="in")
par(mfrow=c(2,1))
plot(x=years_obs, y=lbspr_spr, type="o", ylim=c(0,1), xlim=c(years[1], max(years)), lwd=2, pch=19, 
	xaxs="i", yaxs="i", xpd=NA, xlab="Year", ylab="Overfished Status")
lines(x=years, y=base2$Dt[,"Estimate"], type="o", col="blue", lwd=2, pch=19, xpd=NA)
abline(h=0.3, lty=2)
polygon(x=c(2007,2015,2015,2007), y=c(0,0,0.3,0.3), col="#AA000050", border=NA)

plot(x=years_obs, y=lbspr_fm/0.87, type="o", ylim=c(0,16), xlim=c(years[1], max(years)), lwd=2, pch=19, 
	xaxs="i", yaxs="i", xpd=NA, xlab="Year", ylab="Overfished Status")
lines(x=years, y=base2$Ft[,"Estimate"]/F30, type="o", col="blue", lwd=2, pch=19, xpd=NA)
abline(h=0.3, lty=2)
polygon(x=c(2007,2015,2015,2007), y=c(0,0,0.3,0.3), col="#AA000050", border=NA)
dev.off()

##################################
## decision table
##################################

## run new models 

states_of_nature <- matrix(NA, nrow=3, ncol=2)
colnames(states_of_nature) <- c("M", "h")
states_of_nature[1,] <- c(0.2, 0.5)
states_of_nature[2,] <- c(0.4, 0.6)
states_of_nature[3,] <- c(0.6, 0.8)

## project forward assuming these values have changed
## for each management strategy:
	## rerun base model for each state of nature
	## project into future with that state of nature but F30% equal to the estimated value in this assessment

dtM1S1 <- get_output(model_name="dtM1S1", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	lh=cr_lh2, adjust_param=c("M", "h"), adjust_val=c(0.2, 0.5))




########################################
## forecast
########################################

setwd(init_dir)
source("R_functions\\functions.R")

nproject <- 20
forecast_base <- forecast(currentNa=lowObs$report$N_ta[nrow(lowObs$report$N_ta),], R0=exp(lowObs$report$beta), SB0=lowObs$report$SB0, W_a=cr_lh2$W_a, Mat_a=cr_lh2$Mat_a, S_a=lowObs$report$S_a, h=cr_lh2$h, M=cr_lh2$M, F=F30, nyears=nproject)

forecast_depl <- forecast_base$SBt/lowObs$report$SB0

forecast_base_high <- forecast(currentNa=lowObs$report$N_ta[nrow(lowObs$report$N_ta),], R0=exp(lowObs$report$beta), SB0=lowObs$report$SB0, W_a=cr_lh2$W_a, Mat_a=cr_lh2$Mat_a, S_a=lowObs$report$S_a, h=cr_lh2$h, M=cr_lh2$M, F=F40, nyears=nproject)

forecast_depl_high <- forecast_base_high$SBt/lowObs$report$SB0

forecast_base_low <- forecast(currentNa=lowObs$report$N_ta[nrow(lowObs$report$N_ta),], R0=exp(lowObs$report$beta), SB0=lowObs$report$SB0, W_a=cr_lh2$W_a, Mat_a=cr_lh2$Mat_a, S_a=lowObs$report$S_a, h=cr_lh2$h, M=cr_lh2$M, F=F20, nyears=nproject)

forecast_depl_low <- forecast_base_low$SBt/lowObs$report$SB0

plot(x=years, y=lowObs$Dt[,"Estimate"], ylim=c(0, 1),
	xlim=c(years[1],years[length(years)]+nproject), type="l",
	xlab="Year", ylab="Relative Abundance")
lines(x=seq(years[length(years)], years[length(years)]+nproject-1, by=1),
	y=forecast_depl, lty=2)
lines(x=seq(years[length(years)], years[length(years)]+nproject-1, by=1),
	y=forecast_depl_high, lty=2, col="blue")
lines(x=seq(years[length(years)], years[length(years)]+nproject-1, by=1),
	y=forecast_depl_low, lty=2, col="red")



