### makes sure years are all in the same format

adjust_date <- function(date, out){
	date_type <- grepl("/", date)
	if(date_type){
		yr <- strsplit(date, "/")[[1]][3]
		mo <- strsplit(date, "/")[[1]][1]
		if(nchar(yr)<4) yr <- paste0("20", yr)
	}
	if(date_type==FALSE){
		yr <- date
		mo <- NA
	}

	if(out=="year") return(yr)
	if(out=="month") return(mo)
}