forecast <- function(currentNa, R0, SB0, W_a, Mat_a, S_a, 
	M, SigmaR, iter, F, nyears, minlength=FALSE, maxlength=FALSE, lh,
	ref){

ra <- rep(NA, iter)

for(ii in 1:iter){

	set.seed(ii)
	RecDev <- rnorm(nyears, mean=-SigmaR^2/2, sd=SigmaR)

	Nat <- matrix(NA, nrow=nyears, ncol=length(S_a))
	SBt <- rep(NA, nyears)
	Nat[1,] <- currentNa
	SBt[1] <- sum(Nat[1,-1]*W_a[-1]*Mat_a[-1])

	if(minlength==FALSE & maxlength==FALSE) S_a_proj <- S_a
	if(minlength!=FALSE | maxlength!=FALSE){
		S_a_proj <- rep(0, length(S_a))

		if(minlength!=FALSE){
			Smin <- round(lh$t0-log(1-(minlength/lh$linf))/lh$vbk)
			S_a_proj[Smin:length(S_a_proj)] <- 1
		}
		if(maxlength!=FALSE){
			Smax <- round(lh$t0-log(1-(maxlength/lh$linf))/lh$vbk)
			S_a_proj[(Smax+1):length(S_a_proj)] <- 0
		}
	
	}
	for(y in 2:nyears){
		Nat[y,1] <- R0*exp(RecDev[y])
		for(a in 2:length(W_a)){
			Nat[y,a] <- Nat[y-1,a-1]*exp(-M-S_a_proj[a-1]*F)
		}
		SBt[y] <- sum(Nat[y,-1]*W_a[-1]*Mat_a[-1])
	}

	Nt <- rowSums(Nat[,-1])

	ra[ii] <- SBt[nyears]/SB0
}

	prop_of <- length(which(ra < ref))/iter


	Outs <- NULL
	Outs$ra <- ra
	Outs$prop_of <- prop_of
	return(Outs)

}