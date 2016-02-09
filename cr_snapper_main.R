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

output_dir <- file.path(init_dir, "output")
dir.create(output_dir, showWarnings=FALSE)

retro_dir <- file.path(init_dir, "retrospective")
dir.create(retro_dir, showWarnings=FALSE)

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
obs_high <- fishery_data$obs_high
obs_low <- fishery_data$obs_low
years <- fishery_data$years

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
cpue_yrs <- which(is.na(cpue_bl)==FALSE)
names(cpue_input) <- cpue_yrs

## length frequency
lf_input <- lf[which(rowSums(lf)>0),]
lf_yrs <- which(years %in% rownames(lf)[which(rowSums(lf)>0)])
rownames(lf_input) <- lf_yrs

## catch - approximation based on 300 metric tons/year
catch_input <- 3e8/mean(as.numeric(lg$W_g), na.rm=TRUE)
catch_yrs <- years[length(years)]
names(catch_input) <- catch_yrs

## mean length
meanlen_input <- ml$all_gears
meanlen_yrs <- 1:length(years)
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
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high)

## low obs per year

lowObs <- get_output(model_name="lowObs", dat_avail=c("lengthfreq", "index_total"),
	rec_type=0, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low)

### base with beverton-holt stock recruit 

base_bh <- get_output(model_name="base_bh", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high)


### low obs per year with beverton-holt stock recruit 

bh_lowObs <- get_output(model_name="bh_lowObs", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low)

## remove index

rmIndex <- get_output(model_name="rmIndex", dat_avail=c("lengthfreq"),
	rec_type=1, est_params=c("log_F_t_input", "log_sigma_R", "S50", "CV_l"),
	years=years, index=NULL, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high)

## remove length comp

subML <- get_output(model_name="subML", dat_avail=c("meanlength", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_sigma_R", "log_q_I", "S50", "CV_c"),
	years=years, index=cpue_input, lengthfreq=NULL, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high)


## sensitivity - adjust F1 low

F1low <- get_output(model_name="F1low", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	adjust_param="F1", adjust_val=0.01)

## sensitivity - adjust F1 high

F1high <- get_output(model_name="F1high", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	adjust_param="F1", adjust_val=0.4)

## dome selex

dome <- get_output(model_name="dome", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	adjust_param="dome", adjust_val=0.01)

## fix CVc

fixCVc <- get_output(model_name="fixCVc", dat_avail=c("lengthfreq", "index_total"),
	rec_type=1, est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l"),
	years=years, index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high)


## retrospective

setwd(init_dir)
source("R_functions\\functions.R")

  retro_dir_high <- file.path(retro_dir, "obs_high")
	dir.create(retro_dir_high, showWarnings=FALSE)
  retro_dir_low <- file.path(retro_dir, "obs_low")
	dir.create(retro_dir_low, showWarnings=FALSE)


run_retro(index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	years=years, dat_avail=c("lengthfreq", "index_total"), 
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	rec_type=0, retro_dir=retro_dir_high)

run_retro(index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	years=years, dat_avail=c("lengthfreq", "index_total"), 
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	rec_type=0, retro_dir=retro_dir_low)


run_retro(index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_high,
	years=years, dat_avail=c("lengthfreq", "index_total"), 
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	rec_type=1, retro_dir=retro_dir_high)

run_retro(index=cpue_input, lengthfreq=lf_input, meanlen=meanlen_input, catch=NULL, obs_per_yr=obs_low,
	years=years, dat_avail=c("lengthfreq", "index_total"), 
	est_params=c("log_F_t_input", "log_q_I", "log_sigma_R", "S50", "CV_l", "CV_c"),
	rec_type=1, retro_dir=retro_dir_low)



#### using Beverton-Holt stock recruit curve helps estimate recruitment variation with fewer years of data

plot_retro(param="D", years=years, retro_dir=retro_dir_high, rec_type=0)
plot_retro(param="F", years=years, retro_dir=retro_dir_high, rec_type=0)
plot_retro(param="R", years=years, retro_dir=retro_dir_high, rec_type=0)

plot_retro(param="D", years=years, retro_dir=retro_dir_low, rec_type=0)
plot_retro(param="F", years=years, retro_dir=retro_dir_low, rec_type=0)
plot_retro(param="R", years=years, retro_dir=retro_dir_low, rec_type=0)


plot_retro(param="D", years=years, retro_dir=retro_dir_high, rec_type=1)
plot_retro(param="F", years=years, retro_dir=retro_dir_high, rec_type=1)
plot_retro(param="R", years=years, retro_dir=retro_dir_high, rec_type=1)

plot_retro(param="D", years=years, retro_dir=retro_dir_low, rec_type=1)
plot_retro(param="F", years=years, retro_dir=retro_dir_low, rec_type=1)
plot_retro(param="R", years=years, retro_dir=retro_dir_low, rec_type=1)

