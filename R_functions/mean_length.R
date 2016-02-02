mean_length <- function(data, plot){

	years <- unique(as.numeric(data$Year))[order(unique(as.numeric(data$Year)))]
	ml <- sapply(1:length(years), function(x) mean(data$TL_cm[which(as.numeric(data$Year)==years[x])], na.rm=TRUE))
	  names(ml) <- years

	bl_sub <- data[which(data$Gear=="Bottom Longline"),]
	ml_bl <- sapply(1:length(years), function(x) mean(bl_sub$TL_cm[which(as.numeric(bl_sub$Year)==years[x])], na.rm=TRUE))
	  names(ml_bl) <- years

	g_sub <- data[which(data$Gear=="Gillnet"),]
	ml_g <- sapply(1:length(years), function(x) mean(g_sub$TL_cm[which(as.numeric(g_sub$Year)==years[x])], na.rm=TRUE))
	  names(ml_g) <- years

	exp_sub <- data[which(data$Gear=="Bottom Longline Experimental"),]
	ml_exp <- sapply(1:length(years), function(x) mean(exp_sub$TL_cm[which(as.numeric(exp_sub$Year)==years[x])], na.rm=TRUE))

	## percentage of data by gear
	  gears <- unique(lg$Gear)
		ng_lg <- sapply(1:length(gears), function(x) nrow(lg[which(lg$Gear==gears[x]),])/nrow(lg))
		names(ng_lg) <- gears

	if(plot==TRUE){
		par(mfrow=c(1,1))
		plot(x=years, y=ml, ylim=c(0, max(ml)*1.2), xlab="Year", ylab="TL (cm)", pch=19, type="o", cex=1.2, lwd=2, main="Mean Length in Catch by Gear")
		points(x=years, y=ml_bl, type="o", col="blue", cex=1.4, lwd=2)
		points(x=years, y=ml_g, type="o", col="red", cex=1.4, lwd=2)
		points(x=years, y=ml_exp, type="o", col="forestgreen", cex=1.4, lwd=2)
		legend("bottomright", 
			legend=c("all gears (all observations)", 
				paste0("bottom longline (", round(ng_lg[which(names(ng_lg)=="Bottom Longline")],2), ")"), 
				paste0("gillnet (", round(ng_lg[which(names(ng_lg)=="Gillnet")],2), ")"),
				paste0("experimental longline (", round(ng_lg[which(names(ng_lg)=="Bottom Longline Experimental")],2), ")")),
			col=c("black", "blue", "red", "forestgreen"), cex=1.2, pch=c(19,1,1,1), lwd=2)
	}

	outs <- NULL
	outs$all_gears <- ml
	outs$bottom_longline <- ml_bl
	outs$gillnet <- ml_g
	outs$experimental <- ml_exp
	return(outs)
}