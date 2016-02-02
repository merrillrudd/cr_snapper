## calculates the amount of time the gear was in the water

calc_soak_time <- function(time1, time2, hr=TRUE){
	if(time1==""|time2=="") return(NA)
	hr1 <- as.numeric(strsplit(time1, ":")[[1]][1])
	min1 <- as.numeric(strsplit(time1, ":")[[1]][2])

	hr2 <- as.numeric(strsplit(time2, ":")[[1]][1])
	min2 <- as.numeric(strsplit(time2, ":")[[1]][1])
	if(hr2 < 12) hr2_v2 <- 24 + hr2
	if(hr2 >= 12) hr2_v2 <- hr2

	hr_elapsed <- hr2_v2 - hr1
	min_elapsed <- min2 - min1

	min_return <- hr_elapsed*60 + min_elapsed
	hr_return <- round(min_return/60,2)

	if(hr==TRUE) return(hr_return)
	if(hr==FALSE) return(min_return)
}