rm(list=ls())

#########################
## directories 
#########################

init_dir <- "C:\\Git_Projects\\cr_snapper"

#########################
## read in data + 
## initial adjustments
#########################

data_raw <- read.csv(file.path(init_dir, "cr_snapper_database.csv"), 
	stringsAsFactors=FALSE)

adjust_date <- function(date){
	date_type <- grepl("/", date)
	if(date_type){
		yr <- strsplit(date, "/")[[1]][3]
		if(nchar(yr)<4) yr <- paste0("20", yr)
	}
	if(date_type==FALSE) yr <- date

	return(yr)
}

data_raw$Year <- sapply(1:nrow(data_raw), function(x) adjust_date(data_raw$Date[x]))
data_raw$TL_cm <- as.numeric(data_raw$TL_cm)

#############################
## subset Lutjanus guttatus
#############################

lg <- data_raw[which(data_raw$Species=="Lutjanus guttatus"),]

## category vectors
years <- unique(as.numeric(lg$Year))[order(unique(as.numeric(lg$Year)))]
gears <- unique(lg$Gear)

## annual mean length
ml_lg <- sapply(1:length(years), function(x) mean(lg$TL_cm[which(as.numeric(lg$Year)==years[x])], na.rm=TRUE))
plot(x=years, y=ml_lg, ylim=c(0,max(ml_lg)*1.1), xlab="Year", ylab="TL (cm)", pch=19)

## percentage of data by gear
ng_lg <- sapply(1:length(gears), function(x) nrow(lg[which(lg$Gear==gears[x]),])/nrow(lg))
names(ng_lg) <- gears