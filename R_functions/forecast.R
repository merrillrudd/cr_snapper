forecast <- function(currentNa, R0, SB0, W_a, Mat_a, S_a, h, M, F, nyears){

	Nat <- matrix(NA, nrow=nyears, ncol=length(S_a))
	SBt <- rep(NA, nyears)
	Nat[1,] <- currentNa
	SBt[1] <- sum(Nat[1,-1]*W_a[-1]*Mat_a[-1])
	
	for(i in 2:nyears){
		Nat[i,1] <- (4*h*R0*SBt[i-1])/(SB0*(1-h)+SBt[i-1]*(5*h-1))
		for(a in 2:length(W_a)){
			Nat[i,a] <- Nat[i-1,a-1]*exp(-M-S_a[a-1]*F)
		}
		SBt[i] <- sum(Nat[i,-1]*W_a[-1]*Mat_a[-1])
	}

	Nt <- rowSums(Nat[,-1])

	Outs <- NULL
	Outs$Nt <- Nt
	Outs$SBt <- SBt
	return(Outs)

}