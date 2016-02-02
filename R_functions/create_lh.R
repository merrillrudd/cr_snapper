create_lh <- function(species, sens=FALSE, val=FALSE, selex){

	if(species=="CRSNAP"){

		## life history - from Bystrom thesis
		vbk <- 0.21 ## alternate: 0.37 
		linf <- 64.58 ## alternate: 63.20 
		t0 <- 0 
		M <- 0.43 ## based on vbk
		lwa <- 0.0245
		lwb <- 2.790
		AgeMax <- 23
		F1 <- 0.34


		## growth to estimate
		CVlen <- 0.2

		## recruitment to estimate
		R0 <- 1e6
		h <- 0.7

		## fishery related to estimate
		CVcatch <- 0.1
		qcoef <- 1e-5
		SigmaF <- 0.1

		## selectivity
		S50 <- 4 ## 4.43 fully recruited (39.1 cm TL) (from Bystrom thesis)
		Sslope <- 4
		if(selex=="asymptotic") dome <- 0
		if(selex=="dome") dome <- 0.01
		S_a <- exp(dome*Sslope*(S50-0:AgeMax))/(1+exp(-Sslope*(0:AgeMax - S50)))

		## length bins
		binwidth <- 1
		mids <- seq((binwidth/2), linf*1.2, by=binwidth)
		highs <- mids + (binwidth/2)
		lows <- mids - (binwidth)/2
		
		## related values
    	Lmat <- 34 ## Rojas 2006 Gulf of Nicoya
    	Amat <- round(t0-log(1-(Lmat/linf))/vbk)
    	L_a <- linf*(1-exp(-vbk*(0:AgeMax - t0)))
    	W_a <- lwa*L_a^lwb
    	Mat_a <- 1 / (1 + exp(Amat - 0:AgeMax)) 
	}


	Outs <- NULL
	Outs$vbk <- vbk
	Outs$linf <- linf
	Outs$t0 <- t0
	Outs$M <- M
	Outs$lwa <- lwa
	Outs$lwb <- lwb
	Outs$Lmat <- Lmat
	Outs$F1 <- F1
	Outs$CVlen <- CVlen
	Outs$R0 <- R0
	Outs$h <- h
	Outs$CVcatch <- CVcatch
	Outs$qcoef <- qcoef
	Outs$SigmaF <- SigmaF
	Outs$S50 <- S50
	Outs$Sslope <- Sslope
	Outs$dome <- dome
	Outs$S_a <- S_a
	Outs$AgeMax <- AgeMax
	Outs$Amat <- Amat
	Outs$L_a <- L_a
	Outs$W_a <- W_a
	Outs$Mat_a <- Mat_a


	Outs$binwidth <- binwidth
	Outs$mids <- mids
	Outs$highs <- highs
	Outs$lows <- lows

	return(Outs)
}