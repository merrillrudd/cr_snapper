### makes sure years are all in the same format

adjust_date <- function(date){
	date_type <- grepl("/", date)
	if(date_type){
		yr <- strsplit(date, "/")[[1]][3]
		if(nchar(yr)<4) yr <- paste0("20", yr)
	}
	if(date_type==FALSE) yr <- date

	return(yr)
}