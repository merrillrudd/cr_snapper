catch_effort <- function(data, plot){


years <- unique(as.numeric(data$Year))[order(unique(as.numeric(data$Year)))]
catch <- sapply(1:length(years), function(x) nrow(data[which(data$Year==years[x]),]))
	names(catch) <- years
# plot(x=years, y=catch, type="l", ylim=c(0, max(catch)*1.2), lwd=3)

data$Soak_time_hr <- sapply(1:nrow(data), function(x) calc_soak_time(time1=data$Calado_end_time[x], time2=data$Virado_start_time[x], hr=TRUE))
data$Soak_time_min <- sapply(1:nrow(data), function(x) calc_soak_time(time1=data$Calado_end_time[x], time2=data$Virado_start_time[x], hr=FALSE))

data$Effort_index <- data$Soak_time_hr #* data$Line_n_hooks

## check no efforts are negative
length(which(data$Effort_index < 0))==0

cpue <- rep(NA, length(years))
obs_per_yr <- rep(NA, length(years))
for(i in 1:length(years)){

	subyr <- data[which(data$Year==years[i]),]
	obs_per_yr[i] <- length(unique(subyr$Date))

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

if(plot==TRUE){

	par(mfrow=c(2,1))
	plot(x=years, y=catch, pch=19, ylim=c(0, max(catch, na.rm=TRUE)*1.2))
	plot(x=years, y=cpue, pch=19, ylim=c(0, max(cpue, na.rm=TRUE)*1.2))
}

	Outs <- NULL
	Outs$catch <- catch
	Outs$cpue <- cpue
	Outs$obs_per_yr <- obs_per_yr
	Outs$years <- years
	return(Outs)


}