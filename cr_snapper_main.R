rm(list=ls())

#############################################################
## stock assessment of the Costa Rican spotted rose snapper
## Lutjanus guttatus
## Conducted by Merrill Rudd for Conservation International
## mbrudd@uw.edu
## February 2016
#############################################################



#########################
## notes
#########################

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

## did not end up using, but was explored. 

# gsi <- read.csv(file.path(init_dir, "GSI_2013_females.csv"), stringsAsFactors=FALSE)
# gsi_yrs <- unique(gsi$Year)[order(unique(gsi$Year))]
# avg_gsi <- sapply(1:length(gsi_yrs), function(x) mean(gsi$GSI[which(gsi$Year==gsi_yrs[x])]))
# names(avg_gsi) <- gsi_yrs

# plot(x=gsi$Twt, y=gsi$GSI)

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
## length frequency data
###############################

setwd(init_dir)
source("R_functions\\functions.R")

png(file=file.path(fig_dir, "length_frequency_catch.png"), width=9, height=7, units="in", res=200)
lf <- length_frequency(binwidth=1, linf=cr_lh$linf,
 lmat=cr_lh$Lmat, data=lg, plot=TRUE, weight=TRUE)
dev.off()


#################################
## catch and effort data
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

#############################
## data setup for assessment
#############################

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
# file.remove(paste(version, c(".dll", ".o"), sep="")) ## un-comment to compile on a new machine
# compile(paste0(version, ".cpp")) ## un-comment to compile on a new machine

### base

setwd(init_dir)
source("R_functions\\functions.R")

highObs <- get_output(model_name="highObs", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	lh=cr_lh)

## low obs per year

lowObs <- get_output(model_name="lowObs", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low, 
	lh=cr_lh)

## dome + low obs per year

dome <- get_output(model_name="dome", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low, 
	lh=cr_lh, adjust_param="dome", adjust_val=0.01)

## dome + high obs per year

dome_highobs <- get_output(model_name="dome_highobs", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high, 
	lh=cr_lh, adjust_param="dome", adjust_val=0.01)

## dome + high obs per year + nbh

dome_highobs_nbh <- get_output(model_name="dome_highobs_nbh", dat_avail=c("lengthfreq", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high, 
	lh=cr_lh, adjust_param="dome", adjust_val=0.01)

### base without beverton-holt stock recruit 

nbh_highObs <- get_output(model_name="nbh_highObs", dat_avail=c("lengthfreq", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	lh=cr_lh)

### low obs per year without beverton-holt stock recruit 

nbh_lowObs <- get_output(model_name="nbh_lowObs", dat_avail=c("lengthfreq", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	lh=cr_lh)


## single F and R parameters - not well tested, did not use.

# singleFR <- get_output(model_name="singleFR", dat_avail=c("meanlength", "index_total"),
# 	rec_type=2, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
# 	years=years, index=cpue_input, lengthfreq=NULL, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
# 	lh=cr_lh)


########################################
## compare BH and observation sizes
########################################

names <- c("highObs", "nbh_highObs", "dome_highobs", "dome_highobs_nbh")
aic_bhobs <- calc_AIC(names)

### continue forward with low ESS for length frequency, likely more accurately/appropriately represents uncertainty
### beverton-holt stock-recruit curve
### asymptotic selectivity


################################################
## likelihood profiles for important parameters
################################################

setwd(init_dir)
source("R_functions\\functions.R")

## sigma R

sigR_vec <- seq(0.05, 0.9, length=10)
sigR_prof <- run_lprof(index=cpue_input, lengthfreq=lf_input, catch=NULL, obs_per_yr=obs_high,
	years=years, dat_avail=c("lengthfreq", "index_total"), meanlen=meanlen_input,
	est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	rec_type=0, adjust_param="SigmaR", adjust_val=sigR_vec, lprof_dir=lprof_dir, lh=cr_lh, run=FALSE)
par(mfrow=c(1,1))
plot(x=sigR_vec, y=sigR_prof, pch=19)
sigR_vec[which(sigR_prof==min(sigR_prof))] ## 0.6


## natural mortality

M_vec <- seq(0.05, 1, length=10)
M_prof <- run_lprof(index=cpue_input, lengthfreq=lf_input, catch=NULL, obs_per_yr=obs_high,
	years=years, dat_avail=c("lengthfreq", "index_total"), meanlen=meanlen_input,
	est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	rec_type=0, adjust_param="M", adjust_val=M_vec, lprof_dir=lprof_dir, lh=cr_lh, run=FALSE)
par(mfrow=c(1,1))
plot(x=M_vec, y=M_prof, pch=19)
M_vec[which(M_prof==min(M_prof, na.rm=TRUE))]

## linf 
linf_vec <- seq(58, 70, length=10)
linf_prof <- run_lprof(index=cpue_input, lengthfreq=lf_input, catch=NULL, obs_per_yr=obs_high,
	years=years, dat_avail=c("lengthfreq", "index_total"), meanlen=meanlen_input,
	est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	rec_type=0, adjust_param="linf", adjust_val=linf_vec, lprof_dir=lprof_dir, lh=cr_lh, run=FALSE)
par(mfrow=c(1,1))
plot(x=linf_vec, y=linf_prof, pch=19)
linf_vec[which(linf_prof==min(linf_prof, na.rm=TRUE))]

## vbk 

vbk_vec <- seq(0.15, 0.45, length=10)
vbk_prof <- run_lprof(index=cpue_input, lengthfreq=lf_input, catch=NULL, obs_per_yr=obs_high,
	years=years, dat_avail=c("lengthfreq", "index_total"), meanlen=meanlen_input,
	est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	rec_type=0, adjust_param="vbk", adjust_val=vbk_vec, lprof_dir=lprof_dir, lh=cr_lh, run=TRUE)
par(mfrow=c(1,1))
plot(x=vbk_vec, y=vbk_prof, pch=19)
vbk_vec[which(vbk_prof==min(vbk_prof, na.rm=TRUE))]

## lmat 

lmat_vec <- seq(28, 42, length=10)
lmat_prof <- run_lprof(index=cpue_input, lengthfreq=lf_input, catch=NULL, obs_per_yr=obs_high,
	years=years, dat_avail=c("lengthfreq", "index_total"), meanlen=meanlen_input,
	est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	rec_type=0, adjust_param="Lmat", adjust_val=lmat_vec, lprof_dir=lprof_dir, lh=cr_lh, run=TRUE)
par(mfrow=c(1,1))
plot(x=lmat_vec, y=lmat_prof, pch=19)
lmat_vec[which(lmat_prof==min(lmat_prof, na.rm=TRUE))]


## figure for likelihood profiles

png(file.path(fig_dir, "likelihood_profiles.png"), width=8, height=9, res=200, units="in")
par(mfrow=c(2,2))
plot(x=sigR_vec, y=sigR_prof, pch=19, xlab="Recruitment Variation", ylab="NLL",
	ylim=c(min(sigR_prof)*1.2, min(sigR_prof)*0.2))
abline(v=cr_lh$SigmaR, lty=2)
text(x=1.0, y=min(sigR_prof)*0.25, "A", font=2, cex=1.2)
plot(x=M_vec, y=M_prof, pch=19, xlab="Natural Mortality (M)", ylab="NLL",
	ylim=c(min(M_prof,na.rm=TRUE)*1.2, min(M_prof, na.rm=TRUE)*0.2))
abline(v=cr_lh$M, lty=2)
text(x=1.0, y=min(M_prof,  na.rm=TRUE)*0.25, "B", font=2, cex=1.2)
plot(x=linf_vec, y=linf_prof, pch=19, xlab="Linfinity", ylab="NLL",
	ylim=c(min(linf_prof, na.rm=TRUE)*1.2, min(linf_prof, na.rm=TRUE)*0.2))
abline(v=cr_lh$linf, lty=2)
text(x=70, y=min(linf_prof,  na.rm=TRUE)*0.25, "C", font=2, cex=1.2)
plot(x=vbk_vec, y=vbk_prof, pch=19, xlab="von Bertalanffy k", ylab="NLL",
	ylim=c(min(vbk_prof, na.rm=TRUE)*1.2, min(vbk_prof, na.rm=TRUE)*0.2))
abline(v=cr_lh$vbk, lty=2)
text(x=0.45, y=min(vbk_prof,  na.rm=TRUE)*0.25, "D", font=2, cex=1.2)
mtext(side=2, "Negative Log Likelihood", outer=TRUE)
dev.off()

########################################
## choose base model
########################################

## rerun chosen base case with lowest ll for steepness

setwd(init_dir)
source("R_functions\\functions.R")

base <- get_output(model_name="base", dat_avail=c("lengthfreq", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	lh=cr_lh)



choose <- base
choose_rep <- readRDS(file.path(output_dir, "base", "Report.rds"))
choose_sd <- readRDS(file.path(output_dir, "base", "Sdreport.rds"))


########################################
## per recruit analysis
########################################

F_vec <- seq(0,2, by=0.01)
SBPR <- YPR <- rep(NA, length(F_vec))

base_rep <- readRDS(file.path(init_dir, "output", "base", "Report.rds"))
for(i in 1:length(F_vec)){
	SBPR[i] <- calc_ref(Mat_a=cr_lh$Mat_a, W_a=cr_lh$W_a, M=cr_lh$M, S_a=base_rep$S_a, F=F_vec[i], ref="SBPR")
	YPR[i] <- calc_ref(Mat_a=cr_lh$Mat_a, W_a=cr_lh$W_a, M=cr_lh$M, S_a=base_rep$S_a, F=F_vec[i], ref="YPR")
}	

############################################
## calculate per recruit reference points
############################################
S30 <- 0.30*SBPR[1]
ceiling30 <- F_vec[which(round(SBPR)==ceiling(0.30*SBPR[1]))]
floor30 <- F_vec[which(round(SBPR)==floor(0.30*SBPR[1]))]
F30 <- c(ceiling30, floor30)[which(is.na(c(ceiling30, floor30))==FALSE)]

S40 <- 0.40*SBPR[1]
ceiling40 <- F_vec[which(round(SBPR)==ceiling(0.40*SBPR[1]))]
floor40 <- F_vec[which(round(SBPR)==floor(0.40*SBPR[1]))]
F40 <- c(ceiling40, floor40)[which(is.na(c(ceiling40, floor40))==FALSE)]


## figure for SBPR and YPR

png(file.path(fig_dir, "SPR_YPR.png"), width=6, height=5, res=200, units="in")
par(mfrow=c(1,1), mar=c(5,5,5,5))
plot(x=F_vec, y=SBPR, type="o", ylim=c(0, max(SBPR)), xaxs="i", yaxs="i", pch=19, lwd=2,
	xlab="Fishing Mortality", ylab="Spawning Biomass Per Recruit", yaxt="n", font.lab=2, xpd=NA)
axis(2, at=pretty(c(0, max(SBPR))), las=2)


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


### create figures to show results from base model

FUN = function(InputMat, log=TRUE){

  if(log==TRUE) return(c( exp(InputMat[,1]-1*InputMat[,2]), rev(exp(InputMat[,1]+1*InputMat[,2]))))

  if(log==FALSE) return(c( InputMat[,1]-1*InputMat[,2], rev(InputMat[,1]+1*InputMat[,2])))

} 

## relative abundance and fishing mortality results together

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



## fishing mortality base case results only

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

## relative abundance base case results only

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


##############################################
## kobe plot - overfished and/or overfishing?
##############################################

png(file.path(fig_dir, "base_kobe.png"), width=8, height=7, res=200, units="in")
plot(choose$Dt[,"Estimate"], choose$Ft[,"Estimate"]/F30, xlim=c(0,0.8), ylim=c(0, 1.1), 
	pch=19, col="blue", xaxs="i", yaxs="i", xpd=NA,
	xlab="Relative Abundance", ylab=expression("F/F"[F30]), type="o")
points(choose$Dt[length(years), "Estimate"], choose$Ft[length(years), "Estimate"]/F30, cex=2, col="blue")
text(choose$Dt[length(years), "Estimate"], (choose$Ft[length(years), "Estimate"]/F30)+0.03, "2015")
text(choose$Dt[length(years)-2, "Estimate"]+0.03, (choose$Ft[length(years)-2, "Estimate"])/F30 + 0.03, "2013")
text(choose$Dt[length(years)-4, "Estimate"], (choose$Ft[length(years)-4, "Estimate"])/F30 - 0.03, "2011")
abline(h=1, lty=2)
abline(v=0.3, lty=2)
polygon(x=c(0,1.1,1.1,0), y=c(1,1,1.1,1.1), border=NA, 
	col="#AA000050")
polygon(x=c(0,0.3,0.3,0), y=c(0,0,1.1,1.1), border=NA,
	col="#AA000050")
dev.off()



########################################
## information contained in data types
########################################

## remove index - results from length composition only

rmIndex <- get_output(model_name="rmIndex", dat_avail=c("lengthfreq"),
	rec_type=0, est_params=c("log_F_t_input", "S50", "CV_l"),
	years=years, index=NULL, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	lh=cr_lh)

rmIndex_rep <- readRDS(file.path(output_dir, "rmIndex", "Report.rds"))
rmIndex_sd <- readRDS(file.path(output_dir, "rmIndex", "Sdreport.rds"))

## remove length comp, replace with mean length, still include index

subML <- get_output(model_name="subML", dat_avail=c("meanlength", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_c"),
	years=years, index=cpue_input, lengthfreq=NULL, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	lh=cr_lh)

subML_rep <- readRDS(file.path(output_dir, "subML", "Report.rds"))
subML_sd <- readRDS(file.path(output_dir, "subML", "Sdreport.rds"))


### compare results when index is removed to base case

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

## compare results when mean length is substituted to base case

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
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	adjust_param="F1", adjust_val=0.01,
	lh=cr_lh)

## sensitivity - adjust F1 low - other value

F1low2 <- get_output(model_name="F1low2", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	adjust_param="F1", adjust_val=0.2,
	lh=cr_lh)

## sensitivity - adjust F1 low - other value

F1high <- get_output(model_name="F1high", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	adjust_param="F1", adjust_val=0.45,
	lh=cr_lh)


########################################
## retrospective
########################################

setwd(init_dir)
source("R_functions\\functions.R")

## all data
run_retro(index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	years=years, dat_avail=c("lengthfreq", "index_total"), 
	est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	rec_type=0, retro_dir=retro_dir, lh=cr_lh)

## length comp only
run_retro(index=NULL, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	years=years, dat_avail=c("lengthfreq"), 
	est_params=c("log_F_t_input" "S50", "CV_l"),
	rec_type=0, retro_dir=retro_dir, lh=cr_lh)

## cpue and mean length
run_retro(index=cpue_input, lengthfreq=NULL, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	years=years, dat_avail=c("meanlength", "index_total"), 
	est_params=c("log_F_t_input", "log_q_I", "S50", "CV_c"),
	rec_type=1, retro_dir=retro_dir, lh=cr_lh)


## figures from retrospective analysis

png(file.path(fig_dir, "retro_depl_lengthindex.png"), width=10, height=8, res=200, units="in")
plot_retro(param="D", years=years, retro_dir=retro_dir, rec_type=0, uncertainty=TRUE, dat_avail=c("lengthfreq", "index_total"))
dev.off()

png(file.path(fig_dir, "retro_F_lengthindex.png"), width=10, height=8, res=200, units="in")
plot_retro(param="F", years=years, retro_dir=retro_dir, rec_type=0, uncertainty=TRUE, dat_avail=c("lengthfreq", "index_total"))
dev.off()

plot_retro(param="R", years=years, retro_dir=retro_dir, rec_type=1, uncertainty=TRUE, dat_avail=c("lengthfreq", "index_total"))

png(file.path(fig_dir, "retro_depl_lengthonly.png"), width=10, height=8, res=200, units="in")
plot_retro(param="D", years=years, retro_dir=retro_dir, rec_type=1, uncertainty=TRUE, dat_avail=c("lengthfreq"))
dev.off()

plot_retro(param="F", years=years, retro_dir=retro_dir, rec_type=1, uncertainty=TRUE, dat_avail=c("lengthfreq"))
plot_retro(param="R", years=years, retro_dir=retro_dir, rec_type=1, uncertainty=TRUE, dat_avail=c("lengthfreq"))




########################################
## lb-spr - alternate assessment method
########################################

## from Hordyk et al. 2014
## on github

setwd(init_dir)
source("R_functions\\functions.R")

## stock and fleet objects
StockPars <- NULL
StockPars$MK <- cr_lh$M/cr_lh$vbk
StockPars$Linf <- cr_lh$linf
StockPars$CVLinf <- cr_lh$CVlen
StockPars$MaxSD <- 2
StockPars$L50 <- cr_lh$Lmat
StockPars$L95 <- cr_lh$linf*(1-exp(-cr_lh$vbk*(7-cr_lh$t0)))
StockPars$FecB  <- 3
StockPars$Walpha <- cr_lh$lwa
StockPars$Wbeta <- cr_lh$lwb
StockPars$Steepness <- cr_lh$h
StockPars$Mpow <- 0
StockPars$AgeMax <- cr_lh$AgeMax

SizeBins <- NULL
SizeBins$Linc <- cr_lh$binwidth
SizeBins$ToSize <- max(cr_lh$highs)

FleetPars <- NULL
FleetPars$SL50 <- round(cr_lh$linf*(1-exp(-cr_lh$vbk*(cr_lh$S50 - cr_lh$t0))))
FleetPars$SL95 <- round(cr_lh$linf*(1-exp(-cr_lh$vbk*(cr_lh$S50+1 - cr_lh$t0))))
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

## assuming reference point for F/M ratio is 0.87
plot(lbspr_spr, lbspr_fm/0.87, xlim=c(0,1.1), ylim=c(0,25), 
	pch=19, col="blue", xaxs="i", yaxs="i", xpd=NA,
	xlab="SPR", ylab="F/M = 0.87", type="o")
points(lbspr_spr[length(years)], lbspr_fm[length(years)]/0.87, cex=2, col="blue")
# text(lbspr_spr-0.04, lbspr_fm/0.87+0.3, labels=years)
text(lbspr_spr[(length(years))]+0.01, lbspr_fm[length(years)]+0.04, labels="2014, 2015")
abline(h=1, lty=2)
abline(v=0.3, lty=2)
polygon(x=c(0,17,17,0), y=c(1,1,25,25), border=NA, 
	col="#AA000050")
polygon(x=c(0,0.3,0.3,0), y=c(0,0,25,25), border=NA,
	col="#AA000050")

png(file.path(fig_dir, "compare_methods_overfished.png"), width=8, height=7, res=200, units="in")
plot(x=years_obs, y=lbspr_spr, type="o", ylim=c(0,1), xlim=c(years[1], max(years)), lwd=2, pch=19, 
	xaxs="i", yaxs="i", xpd=NA, xlab="Year", ylab="Overfished Status")
lines(x=years, y=choose$Dt[,"Estimate"], type="o", col="blue", lwd=2, pch=19, xpd=NA)
abline(h=0.3, lty=2)
polygon(x=c(years[1],years[length(years)],years[length(years)],years[1]), y=c(0,0,0.3,0.3), col="#AA000050", border=NA)
dev.off()


## compare overfished status between assessment methods
png(file.path(fig_dir, "compare_models_overfishing.png"), width=8, height=7, res=200, units="in")
plot(x=years_obs, y=lbspr_fm/0.87, type="o", ylim=c(0,16), xlim=c(years[1], max(years)), lwd=2, pch=19, 
	xaxs="i", yaxs="i", xpd=NA, xlab="Year", ylab="Overfished Status")
lines(x=years, y=choose$Ft[,"Estimate"]/F30, type="o", col="blue", lwd=2, pch=19, xpd=NA)
abline(h=0.3, lty=2)
polygon(x=c(years[1],years[length(years)],years[length(years)],years[1]), y=c(0,0,0.3,0.3), col="#AA000050", border=NA)
dev.off()

##################################
## decision table
##################################

## run new models 

Madj <- c(0.2, 0.4, 0.6)
sigRadj <- c(0.2, 0.6, 0.9)

### low sigma R
M1R1 <- get_output(model_name="M1R1", dat_avail=c("lengthfreq", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	lh=cr_lh, adjust_param=c("M", "SigmaR"), adjust_val=c(Madj[1], sigRadj[1]))

M2R1 <- get_output(model_name="M2R1", dat_avail=c("lengthfreq", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	lh=cr_lh, adjust_param=c("M", "SigmaR"), adjust_val=c(Madj[2], sigRadj[1]))

M3R1 <- get_output(model_name="M3R1", dat_avail=c("lengthfreq", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	lh=cr_lh, adjust_param=c("M", "SigmaR"), adjust_val=c(Madj[3], sigRadj[1]))

### med sigma R
M1R2 <- get_output(model_name="M1R2", dat_avail=c("lengthfreq", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	lh=cr_lh, adjust_param=c("M", "SigmaR"), adjust_val=c(Madj[1], sigRadj[2]))

M2R2 <- get_output(model_name="M2R2", dat_avail=c("lengthfreq", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	lh=cr_lh, adjust_param=c("M", "SigmaR"), adjust_val=c(Madj[2], sigRadj[2]))

M3R2 <- get_output(model_name="M3R2", dat_avail=c("lengthfreq", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	lh=cr_lh, adjust_param=c("M", "SigmaR"), adjust_val=c(Madj[3], sigRadj[2]))


### high sigma R
M1R3 <- get_output(model_name="M1R3", dat_avail=c("lengthfreq", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	lh=cr_lh, adjust_param=c("M", "SigmaR"), adjust_val=c(Madj[1], sigRadj[3]))

M2R3 <- get_output(model_name="M2R3", dat_avail=c("lengthfreq", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	lh=cr_lh, adjust_param=c("M", "SigmaR"), adjust_val=c(Madj[2], sigRadj[3]))

M3R3 <- get_output(model_name="M3R3", dat_avail=c("lengthfreq", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	lh=cr_lh, adjust_param=c("M", "SigmaR"), adjust_val=c(Madj[3], sigRadj[3]))

########################################
## stochastic forecast
########################################

setwd(init_dir)
source("R_functions\\functions.R")

nproject <- 20
# forecast_base <- forecast(currentNa=lowObs$report$N_ta[nrow(lowObs$report$N_ta),], R0=exp(lowObs$report$beta), SB0=lowObs$report$SB0, W_a=cr_lh$W_a, Mat_a=cr_lh$Mat_a, S_a=lowObs$report$S_a, h=cr_lh$h, M=cr_lh$M, F=F30, nyears=nproject)
# forecast_depl <- forecast_base$SBt/lowObs$report$SB0


#### harvest at currently-estimated F30 for different states of nature
F1M1R1 <- forecast(currentNa=M1R1$report$N_ta[nrow(M1R1$report$N_ta),], R0=exp(M1R1$report$beta), SB0=M1R1$report$SB0, W_a=M1R1$report$W_a, Mat_a=M1R1$report$Mat_a, S_a=M1R1$report$S_a, 
	M=Madj[1], SigmaR=sigRadj[1], F=F30, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

F1M2R1 <- forecast(currentNa=M2R1$report$N_ta[nrow(M2R1$report$N_ta),], R0=exp(M2R1$report$beta), SB0=M2R1$report$SB0, W_a=M2R1$report$W_a, Mat_a=M2R1$report$Mat_a, S_a=M2R1$report$S_a, 
	M=Madj[2], SigmaR=sigRadj[1], F=F30, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

F1M3R1 <- forecast(currentNa=M3R1$report$N_ta[nrow(M3R1$report$N_ta),], R0=exp(M3R1$report$beta), SB0=M3R1$report$SB0, W_a=M3R1$report$W_a, Mat_a=M3R1$report$Mat_a, S_a=M3R1$report$S_a, 
	M=Madj[3], SigmaR=sigRadj[1], F=F30, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

F1M1R2 <- forecast(currentNa=M1R2$report$N_ta[nrow(M1R2$report$N_ta),], R0=exp(M1R2$report$beta), SB0=M1R2$report$SB0, W_a=M1R2$report$W_a, Mat_a=M1R2$report$Mat_a, S_a=M1R2$report$S_a, 
	M=Madj[1], SigmaR=sigRadj[2], F=F30, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

F1M2R2 <- forecast(currentNa=M2R2$report$N_ta[nrow(M2R2$report$N_ta),], R0=exp(M2R2$report$beta), SB0=M2R2$report$SB0, W_a=M2R2$report$W_a, Mat_a=M2R2$report$Mat_a, S_a=M2R2$report$S_a, 
	M=Madj[2], SigmaR=sigRadj[2], F=F30, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

F1M3R2 <- forecast(currentNa=M3R2$report$N_ta[nrow(M3R2$report$N_ta),], R0=exp(M3R2$report$beta), SB0=M3R2$report$SB0, W_a=M3R2$report$W_a, Mat_a=M3R2$report$Mat_a, S_a=M3R2$report$S_a, 
	M=Madj[3], SigmaR=sigRadj[2], F=F30, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

F1M1R3 <- forecast(currentNa=M1R3$report$N_ta[nrow(M1R3$report$N_ta),], R0=exp(M1R3$report$beta), SB0=M1R3$report$SB0, W_a=M1R3$report$W_a, Mat_a=M1R3$report$Mat_a, S_a=M1R3$report$S_a, 
	M=Madj[1], SigmaR=sigRadj[3], F=F30, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

F1M2R3 <- forecast(currentNa=M2R3$report$N_ta[nrow(M2R3$report$N_ta),], R0=exp(M2R3$report$beta), SB0=M2R3$report$SB0, W_a=M2R3$report$W_a, Mat_a=M2R3$report$Mat_a, S_a=M2R3$report$S_a, 
	M=Madj[2], SigmaR=sigRadj[3], F=F30, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

F1M3R3 <- forecast(currentNa=M3R3$report$N_ta[nrow(M3R3$report$N_ta),], R0=exp(M3R3$report$beta), SB0=M3R3$report$SB0, W_a=M3R3$report$W_a, Mat_a=M3R3$report$Mat_a, S_a=M3R3$report$S_a, 
	M=Madj[3], SigmaR=sigRadj[3], F=F30, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)


#### harvest at currently-estimated F40 for different states of nature
F2M1R1 <- forecast(currentNa=M1R1$report$N_ta[nrow(M1R1$report$N_ta),], R0=exp(M1R1$report$beta), SB0=M1R1$report$SB0, W_a=M1R1$report$W_a, Mat_a=M1R1$report$Mat_a, S_a=M1R1$report$S_a, 
	M=Madj[1], SigmaR=sigRadj[1], F=F40, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

F2M2R1 <- forecast(currentNa=M2R1$report$N_ta[nrow(M2R1$report$N_ta),], R0=exp(M2R1$report$beta), SB0=M2R1$report$SB0, W_a=M2R1$report$W_a, Mat_a=M2R1$report$Mat_a, S_a=M2R1$report$S_a, 
	M=Madj[2], SigmaR=sigRadj[1], F=F40, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

F2M3R1 <- forecast(currentNa=M3R1$report$N_ta[nrow(M3R1$report$N_ta),], R0=exp(M3R1$report$beta), SB0=M3R1$report$SB0, W_a=M3R1$report$W_a, Mat_a=M3R1$report$Mat_a, S_a=M3R1$report$S_a, 
	M=Madj[3], SigmaR=sigRadj[1], F=F40, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

F2M1R2 <- forecast(currentNa=M1R2$report$N_ta[nrow(M1R2$report$N_ta),], R0=exp(M1R2$report$beta), SB0=M1R2$report$SB0, W_a=M1R2$report$W_a, Mat_a=M1R2$report$Mat_a, S_a=M1R2$report$S_a, 
	M=Madj[1], SigmaR=sigRadj[2], F=F40, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

F2M2R2 <- forecast(currentNa=M2R2$report$N_ta[nrow(M2R2$report$N_ta),], R0=exp(M2R2$report$beta), SB0=M2R2$report$SB0, W_a=M2R2$report$W_a, Mat_a=M2R2$report$Mat_a, S_a=M2R2$report$S_a, 
	M=Madj[2], SigmaR=sigRadj[2], F=F40, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

F2M3R2 <- forecast(currentNa=M3R2$report$N_ta[nrow(M3R2$report$N_ta),], R0=exp(M3R2$report$beta), SB0=M3R2$report$SB0, W_a=M3R2$report$W_a, Mat_a=M3R2$report$Mat_a, S_a=M3R2$report$S_a, 
	M=Madj[3], SigmaR=sigRadj[2], F=F40, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

F2M1R3 <- forecast(currentNa=M1R3$report$N_ta[nrow(M1R3$report$N_ta),], R0=exp(M1R3$report$beta), SB0=M1R3$report$SB0, W_a=M1R3$report$W_a, Mat_a=M1R3$report$Mat_a, S_a=M1R3$report$S_a, 
	M=Madj[1], SigmaR=sigRadj[3], F=F40, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

F2M2R3 <- forecast(currentNa=M2R3$report$N_ta[nrow(M2R3$report$N_ta),], R0=exp(M2R3$report$beta), SB0=M2R3$report$SB0, W_a=M2R3$report$W_a, Mat_a=M2R3$report$Mat_a, S_a=M2R3$report$S_a, 
	M=Madj[2], SigmaR=sigRadj[3], F=F40, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

F2M3R3 <- forecast(currentNa=M3R3$report$N_ta[nrow(M3R3$report$N_ta),], R0=exp(M3R3$report$beta), SB0=M3R3$report$SB0, W_a=M3R3$report$W_a, Mat_a=M3R3$report$Mat_a, S_a=M3R3$report$S_a, 
	M=Madj[3], SigmaR=sigRadj[3], F=F40, 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)


#### harvest at status quo with minimum length
S1M1R1 <- forecast(currentNa=M1R1$report$N_ta[nrow(M1R1$report$N_ta),], R0=exp(M1R1$report$beta), SB0=M1R1$report$SB0, W_a=M1R1$report$W_a, Mat_a=M1R1$report$Mat_a, S_a=M1R1$report$S_a, 
	M=Madj[1], SigmaR=sigRadj[1], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=FALSE, lh=cr_lh)

S1M2R1 <- forecast(currentNa=M2R1$report$N_ta[nrow(M2R1$report$N_ta),], R0=exp(M2R1$report$beta), SB0=M2R1$report$SB0, W_a=M2R1$report$W_a, Mat_a=M2R1$report$Mat_a, S_a=M2R1$report$S_a, 
	M=Madj[2], SigmaR=sigRadj[1], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=FALSE, lh=cr_lh)

S1M3R1 <- forecast(currentNa=M3R1$report$N_ta[nrow(M3R1$report$N_ta),], R0=exp(M3R1$report$beta), SB0=M3R1$report$SB0, W_a=M3R1$report$W_a, Mat_a=M3R1$report$Mat_a, S_a=M3R1$report$S_a, 
	M=Madj[3], SigmaR=sigRadj[1], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=FALSE, lh=cr_lh)

S1M1R2 <- forecast(currentNa=M1R2$report$N_ta[nrow(M1R2$report$N_ta),], R0=exp(M1R2$report$beta), SB0=M1R2$report$SB0, W_a=M1R2$report$W_a, Mat_a=M1R2$report$Mat_a, S_a=M1R2$report$S_a, 
	M=Madj[1], SigmaR=sigRadj[2], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=FALSE, lh=cr_lh)

S1M2R2 <- forecast(currentNa=M2R2$report$N_ta[nrow(M2R2$report$N_ta),], R0=exp(M2R2$report$beta), SB0=M2R2$report$SB0, W_a=M2R2$report$W_a, Mat_a=M2R2$report$Mat_a, S_a=M2R2$report$S_a, 
	M=Madj[2], SigmaR=sigRadj[2], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=FALSE, lh=cr_lh)

S1M3R2 <- forecast(currentNa=M3R2$report$N_ta[nrow(M3R2$report$N_ta),], R0=exp(M3R2$report$beta), SB0=M3R2$report$SB0, W_a=M3R2$report$W_a, Mat_a=M3R2$report$Mat_a, S_a=M3R2$report$S_a, 
	M=Madj[3], SigmaR=sigRadj[2], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=FALSE, lh=cr_lh)

S1M1R3 <- forecast(currentNa=M1R3$report$N_ta[nrow(M1R3$report$N_ta),], R0=exp(M1R3$report$beta), SB0=M1R3$report$SB0, W_a=M1R3$report$W_a, Mat_a=M1R3$report$Mat_a, S_a=M1R3$report$S_a, 
	M=Madj[1], SigmaR=sigRadj[3], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=FALSE, lh=cr_lh)

S1M2R3 <- forecast(currentNa=M2R3$report$N_ta[nrow(M2R3$report$N_ta),], R0=exp(M2R3$report$beta), SB0=M2R3$report$SB0, W_a=M2R3$report$W_a, Mat_a=M2R3$report$Mat_a, S_a=M2R3$report$S_a, 
	M=Madj[2], SigmaR=sigRadj[3], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=FALSE, lh=cr_lh)

S1M3R3 <- forecast(currentNa=M3R3$report$N_ta[nrow(M3R3$report$N_ta),], R0=exp(M3R3$report$beta), SB0=M3R3$report$SB0, W_a=M3R3$report$W_a, Mat_a=M3R3$report$Mat_a, S_a=M3R3$report$S_a, 
	M=Madj[3], SigmaR=sigRadj[3], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=FALSE, lh=cr_lh)



#### harvest at status quo with slot limit
S2M1R1 <- forecast(currentNa=M1R1$report$N_ta[nrow(M1R1$report$N_ta),], R0=exp(M1R1$report$beta), SB0=M1R1$report$SB0, W_a=M1R1$report$W_a, Mat_a=M1R1$report$Mat_a, S_a=M1R1$report$S_a, 
	M=Madj[1], SigmaR=sigRadj[1], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=50, lh=cr_lh)

S2M2R1 <- forecast(currentNa=M2R1$report$N_ta[nrow(M2R1$report$N_ta),], R0=exp(M2R1$report$beta), SB0=M2R1$report$SB0, W_a=M2R1$report$W_a, Mat_a=M2R1$report$Mat_a, S_a=M2R1$report$S_a, 
	M=Madj[2], SigmaR=sigRadj[1], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=50, lh=cr_lh)

S2M3R1 <- forecast(currentNa=M3R1$report$N_ta[nrow(M3R1$report$N_ta),], R0=exp(M3R1$report$beta), SB0=M3R1$report$SB0, W_a=M3R1$report$W_a, Mat_a=M3R1$report$Mat_a, S_a=M3R1$report$S_a, 
	M=Madj[3], SigmaR=sigRadj[1], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=50, lh=cr_lh)

S2M1R2 <- forecast(currentNa=M1R2$report$N_ta[nrow(M1R2$report$N_ta),], R0=exp(M1R2$report$beta), SB0=M1R2$report$SB0, W_a=M1R2$report$W_a, Mat_a=M1R2$report$Mat_a, S_a=M1R2$report$S_a, 
	M=Madj[1], SigmaR=sigRadj[2], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=50, lh=cr_lh)

S2M2R2 <- forecast(currentNa=M2R2$report$N_ta[nrow(M2R2$report$N_ta),], R0=exp(M2R2$report$beta), SB0=M2R2$report$SB0, W_a=M2R2$report$W_a, Mat_a=M2R2$report$Mat_a, S_a=M2R2$report$S_a, 
	M=Madj[2], SigmaR=sigRadj[2], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=50, lh=cr_lh)

S2M3R2 <- forecast(currentNa=M3R2$report$N_ta[nrow(M3R2$report$N_ta),], R0=exp(M3R2$report$beta), SB0=M3R2$report$SB0, W_a=M3R2$report$W_a, Mat_a=M3R2$report$Mat_a, S_a=M3R2$report$S_a, 
	M=Madj[3], SigmaR=sigRadj[2], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=50, lh=cr_lh)

S2M1R3 <- forecast(currentNa=M1R3$report$N_ta[nrow(M1R3$report$N_ta),], R0=exp(M1R3$report$beta), SB0=M1R3$report$SB0, W_a=M1R3$report$W_a, Mat_a=M1R3$report$Mat_a, S_a=M1R3$report$S_a, 
	M=Madj[1], SigmaR=sigRadj[3], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=50, lh=cr_lh)

S2M2R3 <- forecast(currentNa=M2R3$report$N_ta[nrow(M2R3$report$N_ta),], R0=exp(M2R3$report$beta), SB0=M2R3$report$SB0, W_a=M2R3$report$W_a, Mat_a=M2R3$report$Mat_a, S_a=M2R3$report$S_a, 
	M=Madj[2], SigmaR=sigRadj[3], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=50, lh=cr_lh)

S2M3R3 <- forecast(currentNa=M3R3$report$N_ta[nrow(M3R3$report$N_ta),], R0=exp(M3R3$report$beta), SB0=M3R3$report$SB0, W_a=M3R3$report$W_a, Mat_a=M3R3$report$Mat_a, S_a=M3R3$report$S_a, 
	M=Madj[3], SigmaR=sigRadj[3], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=40, maxlength=50, lh=cr_lh)


#### harvest at status quo 
SqM1R1 <- forecast(currentNa=M1R1$report$N_ta[nrow(M1R1$report$N_ta),], R0=exp(M1R1$report$beta), SB0=M1R1$report$SB0, W_a=M1R1$report$W_a, Mat_a=M1R1$report$Mat_a, S_a=M1R1$report$S_a, 
	M=Madj[1], SigmaR=sigRadj[1], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

SqM2R1 <- forecast(currentNa=M2R1$report$N_ta[nrow(M2R1$report$N_ta),], R0=exp(M2R1$report$beta), SB0=M2R1$report$SB0, W_a=M2R1$report$W_a, Mat_a=M2R1$report$Mat_a, S_a=M2R1$report$S_a, 
	M=Madj[2], SigmaR=sigRadj[1], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

SqM3R1 <- forecast(currentNa=M3R1$report$N_ta[nrow(M3R1$report$N_ta),], R0=exp(M3R1$report$beta), SB0=M3R1$report$SB0, W_a=M3R1$report$W_a, Mat_a=M3R1$report$Mat_a, S_a=M3R1$report$S_a, 
	M=Madj[3], SigmaR=sigRadj[1], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

SqM1R2 <- forecast(currentNa=M1R2$report$N_ta[nrow(M1R2$report$N_ta),], R0=exp(M1R2$report$beta), SB0=M1R2$report$SB0, W_a=M1R2$report$W_a, Mat_a=M1R2$report$Mat_a, S_a=M1R2$report$S_a, 
	M=Madj[1], SigmaR=sigRadj[2], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

SqM2R2 <- forecast(currentNa=M2R2$report$N_ta[nrow(M2R2$report$N_ta),], R0=exp(M2R2$report$beta), SB0=M2R2$report$SB0, W_a=M2R2$report$W_a, Mat_a=M2R2$report$Mat_a, S_a=M2R2$report$S_a, 
	M=Madj[2], SigmaR=sigRadj[2], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

SqM3R2 <- forecast(currentNa=M3R2$report$N_ta[nrow(M3R2$report$N_ta),], R0=exp(M3R2$report$beta), SB0=M3R2$report$SB0, W_a=M3R2$report$W_a, Mat_a=M3R2$report$Mat_a, S_a=M3R2$report$S_a, 
	M=Madj[3], SigmaR=sigRadj[2], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

SqM1R3 <- forecast(currentNa=M1R3$report$N_ta[nrow(M1R3$report$N_ta),], R0=exp(M1R3$report$beta), SB0=M1R3$report$SB0, W_a=M1R3$report$W_a, Mat_a=M1R3$report$Mat_a, S_a=M1R3$report$S_a, 
	M=Madj[1], SigmaR=sigRadj[3], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

SqM2R3 <- forecast(currentNa=M2R3$report$N_ta[nrow(M2R3$report$N_ta),], R0=exp(M2R3$report$beta), SB0=M2R3$report$SB0, W_a=M2R3$report$W_a, Mat_a=M2R3$report$Mat_a, S_a=M2R3$report$S_a, 
	M=Madj[2], SigmaR=sigRadj[3], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)

SqM3R3 <- forecast(currentNa=M3R3$report$N_ta[nrow(M3R3$report$N_ta),], R0=exp(M3R3$report$beta), SB0=M3R3$report$SB0, W_a=M3R3$report$W_a, Mat_a=M3R3$report$Mat_a, S_a=M3R3$report$S_a, 
	M=Madj[3], SigmaR=sigRadj[3], F=choose$Ft[nrow(choose$Ft),"Estimate"], 
	nyears=nproject, iter=100, ref=0.3, 
	minlength=FALSE, maxlength=FALSE, lh=cr_lh)



