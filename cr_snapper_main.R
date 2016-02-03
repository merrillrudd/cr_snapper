rm(list=ls())

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

fig_dir <- file.path(init_dir, "figs")
dir.create(fig_dir, showWarnings=FALSE)


setwd(init_dir)
source("R_functions\\functions.R")

#########################
## setup life history
#########################

cr_lh <- create_lh(species="CRSNAP", selex="asymptotic")
cr_lh_dome <- create_lh(species="CRSNAP", selex="dome")

#########################
## read in data + 
## initial adjustments
#########################

data_raw <- read.csv(file.path(init_dir, "cr_snapper_database.csv"), 
	stringsAsFactors=FALSE)

data_raw$Year <- sapply(1:nrow(data_raw), function(x) adjust_date(data_raw$Date[x]))
data_raw$TL_cm <- as.numeric(data_raw$TL_cm)

#############################
## subset Lutjanus guttatus
#############################

lg <- data_raw[which(data_raw$Species=="Lutjanus guttatus" & data_raw$Year>1950),]
write.csv(lg, file.path(init_dir, "cr_snapper_filtered.csv"))

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
obs_per_yr <- fishery_data$obs_per_yr
years <- fishery_data$years

par(mfrow=c(2,1))
plot(cpue_bl, pch=19, ylim=c(0, max(cpue_bl, na.rm=TRUE)))
plot(cpue_g, pch=19, ylim=c(0, max(cpue_g, na.rm=TRUE)))

#######################
## data setup
#######################

setwd(init_dir)
source("R_functions\\functions.R")

cpue_yrs <- which(is.na(cpue_bl)==FALSE)
lf_years <- which(years %in% rownames(lf)[which(rowSums(lf)>0)])
catch_yrs <- years

catch_input <- (rep(300, length(years))*907185)/mean(as.numeric(lg$W_g), na.rm=TRUE)
cpue_input <- cpue_bl[which(is.na(cpue_bl)==FALSE)] ## choose bottom longline cpue only
lf_input <- lf[which(rowSums(lf)>0),]
obs_input <- obs_per_yr

names(catch_input) <- catch_yrs
names(cpue_input) <- cpue_yrs
rownames(lf_input) <- lf_years

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

##########
## base
##########

model_name <- "base"
dat_input <- create_inputs(simdir=model_name, param=FALSE, val=FALSE, lh_dat=cr_lh)

model_dir <- file.path(init_dir, model_name)
if(file.exists(model_dir)) unlink(model_dir, TRUE)
dir.create(model_dir, showWarnings=FALSE)

# lh <- dat_input
# catch <- catch_input
# index <- cpue_input
# lengthfreq <- lf_input
# obs_per_yr <- obs_input
# dat_avail <- c("lengthfreq_multiple", "catch_total", "index_total")
# version <- "lb_statespace"
# obs_meanlen=ml$all_gears
# RecType <- 0
# est_params=c("log_F_t_input", "log_q_I", "beta", "log_sigma_R", "S50", "CV_l", "CV_c")

setwd(init_dir)
source("R_functions\\functions.R")

base <- run_statespace(lh=dat_input, years=years, catch=NULL, 
	index=cpue_input, lengthfreq=lf_input, obs_per_yr=obs_input, 
	dat_avail=c("lengthfreq_multiple", "index_total"),
	model_name=model_name, model_dir=model_dir, obs_meanlen=ml$all_gears, RecType=0,
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"))

base_rep <- readRDS(file.path(init_dir, model_name, "Report.rds"))


## sensitivity - adjust F1 high

model_name <- "sensitivity_R0low"
dat_input <- create_inputs(simdir=model_name, param="F1", val=0.5, lh_dat=cr_lh)

model_dir <- file.path(init_dir, model_name)
if(file.exists(model_dir)) unlink(model_dir, TRUE)
dir.create(model_dir, showWarnings=FALSE)

sens_R0low <- run_statespace(lh=dat_input, years=years, catch=NULL, 
	index=cpue_input, lengthfreq=lf_input, obs_per_yr=obs_input, 
	dat_avail=c("lengthfreq_multiple", "index_total"),
	model_name=model_name, model_dir=model_dir, obs_meanlen=ml$all_gears, RecType=0,
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"))

sens_R0low_rep <- readRDS(file.path(init_dir, model_name, "Report.rds"))


## sensitivity - adjust F1 low

model_name <- "sensitivity_F1low"
dat_input <- create_inputs(simdir=model_name, param="F1", val=0.25, lh_dat=cr_lh)

model_dir <- file.path(init_dir, model_name)
if(file.exists(model_dir)) unlink(model_dir, TRUE)
dir.create(model_dir, showWarnings=FALSE)

sens_F1low <- run_statespace(lh=dat_input, years=years, catch=NULL, 
	index=cpue_input, lengthfreq=lf_input, obs_per_yr=obs_input, 
	dat_avail=c("lengthfreq_multiple", "index_total"),
	model_name=model_name, model_dir=model_dir, obs_meanlen=ml$all_gears, RecType=0,
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"))

sens_F1low_rep <- readRDS(file.path(init_dir, model_name, "Report.rds"))


## sensitivity - adjust R0 higher

model_name <- "sensitivity_R0high"
dat_input <- create_inputs(simdir=model_name, param="R0", val=1e10, lh_dat=cr_lh)

model_dir <- file.path(init_dir, model_name)
if(file.exists(model_dir)) unlink(model_dir, TRUE)
dir.create(model_dir, showWarnings=FALSE)

sens_R0high <- run_statespace(lh=dat_input, years=years, catch=NULL, 
	index=cpue_input, lengthfreq=lf_input, obs_per_yr=obs_input, 
	dat_avail=c("lengthfreq_multiple", "index_total"),
	model_name=model_name, model_dir=model_dir, obs_meanlen=ml$all_gears, RecType=0,
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"))

sens_R0high_rep <- readRDS(file.path(init_dir, model_name, "Report.rds"))

## sensitivity - adjust R0 lower

model_name <- "sensitivity_R0low"
dat_input <- create_inputs(simdir=model_name, param="R0", val=1e6, lh_dat=cr_lh)

model_dir <- file.path(init_dir, model_name)
if(file.exists(model_dir)) unlink(model_dir, TRUE)
dir.create(model_dir, showWarnings=FALSE)

sens_R0low <- run_statespace(lh=dat_input, years=years, catch=NULL, 
	index=cpue_input, lengthfreq=lf_input, obs_per_yr=obs_input, 
	dat_avail=c("lengthfreq_multiple", "index_total"),
	model_name=model_name, model_dir=model_dir, obs_meanlen=ml$all_gears, RecType=0,
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"))

sens_R0low_rep <- readRDS(file.path(init_dir, model_name, "Report.rds"))


##########
## base assuming dome-shaped selex
##########

model_name <- "base_dome"
dat_input <- create_inputs(simdir=model_name, param=FALSE, val=FALSE, lh_dat=cr_lh_dome)

model_dir <- file.path(init_dir, model_name)
if(file.exists(model_dir)) unlink(model_dir, TRUE)
dir.create(model_dir, showWarnings=FALSE)

base_dome <- run_statespace(lh=dat_input, years=years, catch=NULL, 
	index=cpue_input, lengthfreq=lf_input, obs_per_yr=obs_input, 
	dat_avail=c("lengthfreq_multiple", "index_total"),
	model_name=model_name, model_dir=model_dir, obs_meanlen=ml$all_gears, RecType=0,
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"))

base_dome_rep <- readRDS(file.path(init_dir, model_name, "Report.rds"))








## base with BH recruitment

setwd(init_dir)
source("R_functions\\functions.R")

model_name <- "base_bh"
model_dir <- file.path(init_dir, model_name)
if(file.exists(model_dir)) unlink(model_dir, TRUE)
dir.create(model_dir, showWarnings=FALSE)

base_bh <- run_statespace(lh=dat_input, years=years, catch=catch_input, 
	index=cpue_input, lengthfreq=lf_input, obs_per_yr=obs_input, 
	dat_avail=c("lengthfreq_multiple", "catch_total", "index_total"),
	model_name=model_name, model_dir=model_dir, obs_meanlen=ml$all_gears, RecType=1,
	est_params=c("log_F_t_input", "log_q_I", "beta", "log_sigma_R", "S50", "CV_l", "CV_c"))

base_bh_rep <- readRDS(file.path(init_dir, model_name, "Report.rds"))

## base with BH recruitment fix CVc

setwd(init_dir)
source("R_functions\\functions.R")

model_name <- "base_bh_fixCVc"
model_dir <- file.path(init_dir, model_name)
if(file.exists(model_dir)) unlink(model_dir, TRUE)
dir.create(model_dir, showWarnings=FALSE)

base_bh_fixCVc <- run_statespace(lh=dat_input, years=years, catch=catch_input, 
	index=cpue_input, lengthfreq=lf_input, obs_per_yr=obs_input, 
	dat_avail=c("lengthfreq_multiple", "catch_total", "index_total"),
	model_name=model_name, model_dir=model_dir, obs_meanlen=ml$all_gears, RecType=1,
	est_params=c("log_F_t_input", "log_q_I", "beta", "log_sigma_R", "S50", "CV_l"))

base_bh_fixCVc_rep <- readRDS(file.path(init_dir, model_name, "Report.rds"))


## remove catch

setwd(init_dir)
source("R_functions\\functions.R")

model_name <- "rm_catch"
dat_input <- create_inputs(simdir=model_name, param=FALSE, val=FALSE, lh_dat=cr_lh)

model_dir <- file.path(init_dir, model_name)
if(file.exists(model_dir)) unlink(model_dir, TRUE)
dir.create(model_dir, showWarnings=FALSE)

setwd(init_dir)
source("R_functions\\functions.R")

rm_catch <- run_statespace(lh=dat_input, years=years, catch=NULL, 
	index=cpue_input, lengthfreq=lf_input, obs_per_yr=obs_input, 
	dat_avail=c("lengthfreq_multiple", "index_total"),
	model_name=model_name, model_dir=model_dir, obs_meanlen=ml$all_gears, RecType=0,
	est_params=c("log_F_t_input", "log_q_I", "beta", "log_sigma_R", "S50", "CV_l", "CV_c"))

rm_catch_rep <- readRDS(file.path(init_dir, model_name, "Report.rds"))

