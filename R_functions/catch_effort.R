catch_effort <- function(data, sep_index){


years <- unique(as.numeric(data$Year))[order(unique(as.numeric(data$Year)))]
catch <- sapply(1:length(years), function(x) nrow(data[which(data$Year==years[x]),]))
	names(catch) <- years
# plot(x=years, y=catch, type="l", ylim=c(0, max(catch)*1.2), lwd=3)

data$Soak_time_hr <- sapply(1:nrow(data), function(x) calc_soak_time(time1=data$Calado_end_time[x], time2=data$Virado_start_time[x], hr=TRUE))
data$Soak_time_min <- sapply(1:nrow(data), function(x) calc_soak_time(time1=data$Calado_end_time[x], time2=data$Virado_start_time[x], hr=FALSE))

if(sep_index==TRUE){
 data$Effort_index <- NA
 gillnet_index <- which(data$Gear=="Gillnet")
   data$Effort_index[gillnet_index] <- data$Soak_time_hr[gillnet_index]
 longline_index <- which(data$Gear=="Bottom Longline")
   data$Effort_index[longline_index] <- data$Soak_time_hr[longline_index] * data$Line_n_hooks[longline_index]
}
if(sep_index==FALSE){
	data$Effort_index <- data$Soak_time_hr
}

## check no efforts are negative
if(length(which(data$Effort_index[gillnet_index] < 0))!=0) stop("gillnet effort indices < 0")
if(length(which(data$Effort_index[longline_index] < 0))!=0) stop("longline effort indices < 0")

obs_high <- obs_low <- rep(NA, length(years))
for(i in 1:length(years)){
	subyr <- data[which(data$Year==years[i]),]
	obs_high[i] <- length(unique(subyr$Date))
	obs_low[i] <- length(unique(subyr$Month))
}


if(sep_index==FALSE){
    cpue <- rep(NA, length(years))
    for(i in 1:length(years)){    

    	subyr <- data[which(data$Year==years[i]),]

    	dates <- unique(subyr$Date)
    	catch_pass <- effort_pass <- NULL
    	for(j in 1:length(dates)){
    		subdate <- subyr[which(subyr$Date==dates[j] & subyr$Observation_type=="Onboard" & is.na(subyr$Effort_index)==FALSE),]
    		passes <- unique(subdate$Pass)
    		for(k in 1:length(passes)){
    			subpass <- subdate[which(subdate$Pass==passes[k] & is.na(subdate$Effort_index)==FALSE),]
    			if(length(unique(subpass$Effort_index))>1) stop("effort varies within pass")
    			if(length(unique(subpass$Effort_index))==1){
    				effort_pass <- c(effort_pass, unique(subpass$Effort_index))
    				catch_pass <- c(catch_pass, nrow(subpass))
    			}
    		}
    	}
    	
    	cpue[i] <- mean(catch_pass/effort_pass, na.rm=TRUE)    

    	rm(catch_pass)
    	rm(effort_pass)
    }
    names(cpue) <- years
}

if(sep_index==TRUE){
	cpue_bl <- rep(NA, length(years))
for(i in 1:length(years)){

	subyr <- data[which(data$Year==years[i] & data$Gear=="Bottom Longline"),]

	dates <- unique(subyr$Date)
	catch_pass <- effort_pass <- NULL
	for(j in 1:length(dates)){
		subdate <- subyr[which(subyr$Date==dates[j] & subyr$Observation_type=="Onboard" & is.na(subyr$Effort_index)==FALSE),]
		passes <- unique(subdate$Pass)
		for(k in 1:length(passes)){
			subpass <- subdate[which(subdate$Pass==passes[k] & is.na(subdate$Effort_index)==FALSE),]
			if(length(unique(subpass$Effort_index))>1) stop("effort varies within pass")
			if(length(unique(subpass$Effort_index))==1){
				effort_pass <- c(effort_pass, unique(subpass$Effort_index))
				catch_pass <- c(catch_pass, nrow(subpass))
			}
		}
	}
	
	cpue_bl[i] <- mean(catch_pass/effort_pass, na.rm=TRUE)

	rm(catch_pass)
	rm(effort_pass)
}
names(cpue_bl) <- years

cpue_g <- rep(NA, length(years))
for(i in 1:length(years)){

	subyr <- data[which(data$Year==years[i] & data$Gear=="Gillnet"),]

	dates <- unique(subyr$Date)
	catch_pass <- effort_pass <- NULL
	for(j in 1:length(dates)){
		subdate <- subyr[which(subyr$Date==dates[j] & subyr$Observation_type=="Onboard" & is.na(subyr$Effort_index)==FALSE),]
		passes <- unique(subdate$Pass)
		for(k in 1:length(passes)){
			subpass <- subdate[which(subdate$Pass==passes[k] & is.na(subdate$Effort_index)==FALSE),]
			if(length(unique(subpass$Effort_index))>1) stop("effort varies within pass")
			if(length(unique(subpass$Effort_index))==1){
				effort_pass <- c(effort_pass, unique(subpass$Effort_index))
				catch_pass <- c(catch_pass, nrow(subpass))
			}
		}
	}
	
	cpue_g[i] <- mean(catch_pass/effort_pass, na.rm=TRUE)

	rm(catch_pass)
	rm(effort_pass)
}
names(cpue_g) <- years

}


	Outs <- NULL
	Outs$catch <- catch
	if(sep_index==FALSE) Outs$cpue <- cpue
	Outs$obs_high <- obs_high
	Outs$obs_low <- obs_low
	Outs$years <- years
	if(sep_index==TRUE){
		Outs$cpue_bl <- cpue_bl
		Outs$cpue_g <- cpue_g
	}
	return(Outs)


}