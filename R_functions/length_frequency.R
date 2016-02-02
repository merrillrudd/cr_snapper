length_frequency <- function(binwidth, linf, lmat, data, plot, weight){

	bins <- seq(1, linf*1.2, by=binwidth)

	years <- unique(as.numeric(data$Year))[order(unique(as.numeric(data$Year)))]

	data_sub <- data[which(data$Gear %in% c("Bottom Longline", "Gillnet")),]
	data_sub$TL_round <- binwidth*round(data_sub$TL_cm/binwidth, 0)
	data$TL_round <- binwidth*round(data$TL_cm/binwidth,0)

	lf <- lf_bl <- lf_g <- matrix(nrow=length(years), ncol=length(bins))
	rownames(lf_bl) <- rownames(lf_g) <- years
	colnames(lf_bl) <- colnames(lf_g) <- bins

	## length frequency
	for(i in 1:length(years)){
		lf[i,] <- sapply(1:length(bins), function(x) length(which(data$TL_round==bins[x] & data$Year==years[i])))
		lf_bl[i,] <- sapply(1:length(bins), function(x) length(which(data_sub$TL_round==bins[x] & data_sub$Year==years[i] & data_sub$Gear=="Bottom Longline")))
		lf_g[i,] <- sapply(1:length(bins), function(x) length(which(data_sub$TL_round==bins[x] & data_sub$Year==years[i] & data_sub$Gear=="Gillnet")))
	}

	## length frequency proportions
	lf_prop <- t(apply(lf, 1, function(x) x/sum(x)))
	lf_bl_prop <- t(apply(lf_bl, 1, function(x) x/sum(x)))
	lf_g_prop <- t(apply(lf_g, 1, function(x) x/sum(x)))

	lf_bl_prop[which(is.na(lf_bl_prop))] <- 0
	lf_g_prop[which(is.na(lf_g_prop))] <- 0

	## catch by gear in each year
	catch_bl <- sapply(1:length(years), function(x) nrow(data_sub[which(data_sub$Gear=="Bottom Longline" & data_sub$Year==years[x]),]))
	catch_g <- sapply(1:length(years), function(x) nrow(data_sub[which(data_sub$Gear=="Gillnet" & data_sub$Year==years[x]),]))

	lf_bl_catch <- t(sapply(1:length(years), function(x) lf_bl_prop[x,]*catch_bl[x]))
	lf_g_catch <- t(sapply(1:length(years), function(x) lf_g_prop[x,]*catch_g[x]))

	## sum length bins by gears
	lf_total <- lf_bl_catch + lf_g_catch
	lf_total_prop <- t(sapply(1:length(years), function(x) lf_total[x,]/sum(lf_total[x,])))
		rownames(lf_total_prop) <- years

	ylab_seq <- seq(1,length(years), by=3)
	xlab_seq <- (length(years)-2):length(years)

	if(plot==TRUE){
		par(mfrow=c(3,3), mar=c(0,0,0,0), omi=c(1,1,0.5,0.5))
		for(i in 1:length(years)){
			plot(lf_prop[i,], type="h", xlim=c(0, max(bins)), 
				ylim=c(0, 0.1), xaxs="i", yaxs="i", xaxt="n", yaxt="n", lwd=4)
			par(new=TRUE)
			plot(lf_bl_prop[i,], type="h", xlim=c(0, max(bins)),
				ylim=c(0, 0.1), xaxs="i", yaxs="i", xaxt="n", yaxt="n", lwd=4, col="goldenrod")
			abline(v=lmat, col="blue", lwd=2)
			if(i %in% ylab_seq) axis(2, seq(0, 0.1, by=0.02), las=2, cex.axis=1.2)
			if(i %in% xlab_seq) axis(1, seq(0, 80, by=20), cex.axis=1.2)
			text(x=10, y=0.09, years[i], font=2, cex=1.2)
		}
			mtext(side=3, outer=TRUE, "Length Frequency of Catch", line=1, cex=1.4)
			mtext(side=1, outer=TRUE, "Total Length (cm)", line=3.5, cex=1.2)
			mtext(side=2, outer=TRUE, "Proportion in Catch", line=3.5, cex=1.2)
			legend("topright", legend=c("all gears", "longline", "maturity"), 
				col=c("black", "goldenrod", "blue"), lty="solid", cex=1.3, lwd=3)
	}

	if(weight==TRUE) return(lf_total_prop)
	if(weight==FALSE) return(lf_prop)

}