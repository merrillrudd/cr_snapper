LBSPRSim <- function(StockPars, FleetPars, SizeBins=NULL, P=0.001, Nage=201) {

  MK <- StockPars$MK 
  Linf <- StockPars$Linf
  CVLinf <- StockPars$CVLinf 
  L50 <- StockPars$L50 
  L95 <- StockPars$L95 
  Beta <- StockPars$FecB 
  MaxSD <- StockPars$MaxSD
  
  SDLinf <- CVLinf * Linf # Standard Deviation of Length-at-Age # Assumed constant CV here
  if (is.null(SizeBins)) {
    SizeBins$Linc <- 5
	 SizeBins$ToSize <- Linf + MaxSD * SDLinf
  }
  if (is.null(SizeBins$ToSize)) SizeBins$ToSize <- Linf + MaxSD * SDLinf
  Linc <- SizeBins$Linc 
  ToSize <- SizeBins$ToSize
  
  FM <- FleetPars$FM 
  SL50 <- FleetPars$SL50 
  SL95 <- FleetPars$SL95 
  
  LenBins <- seq(from=0, by=Linc, to=ToSize)	
  LenMids <- seq(from=0.5*Linc, by=Linc, length.out=length(LenBins)-1)
  x <- seq(from=0, to=1, length.out=Nage) # relative age vector
  EL <- (1-P^(x/MK)) * Linf # length at relative age 
  rLens <- EL/Linf # relative length 
  SDL <- EL * CVLinf # standard deviation of length-at-age
  
  Nlen <- length(LenMids) 
  Prob <- matrix(NA, nrow=Nage, ncol=Nlen)
  Prob[,1] <- pnorm((LenBins[2] - EL)/SDL, 0, 1) # probablility of length-at-age
  for (i in 2:(Nlen-1)) {
    Prob[,i] <- pnorm((LenBins[i+1] - EL)/SDL, 0, 1) - 
		pnorm((LenBins[i] - EL)/SDL, 0, 1)
  }
  Prob[,Nlen] <- 1 - pnorm((LenBins[Nlen] - EL)/SDL, 0, 1)
  
  # Truncate normal dist at MaxSD 
  mat <- array(1, dim=dim(Prob))
  for (X in 1:Nage) {
    ind <- which(abs((LenMids - EL[X]) /SDL[X]) >= MaxSD)
    mat[X,ind] <- 0
  }
  
  Prob <- Prob * mat

  SL <- 1/(1+exp(-log(19)*(LenMids-SL50)/(SL95-SL50))) # Selectivity at length
  Sx <- apply(t(Prob) * SL, 2, sum) # Selectivity at relative age 
  MSX <- cumsum(Sx) / seq_along(Sx) # Mean cumulative selectivity for each age 
  Ns <- (1-rLens)^(MK+(MK*FM)*MSX) # number at relative age in population
  
  Cx <- t(t(Prob) * SL) # Conditional catch length-at-age probablilities  
  Nc <- apply(Ns * Cx, 2, sum) # 
  Pop <- apply(Ns * Prob, 2, sum)
  
  Ml <- 1/(1+exp(-log(19)*(LenMids-L50)/(L95-L50))) # Maturity at length
  Ma <-  apply(t(Prob) * Ml, 2, sum) # Maturity at relative age 
  
  N0 <- (1-rLens)^MK # Unfished numbers-at-age 
  SPR <- sum(Ma * Ns * rLens^Beta)/sum(Ma * N0 * rLens^Beta)
  
  Output <- NULL 
  Output$SPR <- SPR 
  Output$LenMids <- LenMids
  Output$PropLen <- Nc/sum(Nc)
  Output$Pop <- Pop
  
  Output$LCatchFished <- Nc/sum(Nc)
  Output$LPopFished <- Pop
  Output$LCatchUnfished <- apply(N0 * Cx, 2, sum)

  return(Output)
}  
