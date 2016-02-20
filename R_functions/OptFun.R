OptFun <- function(tryFleetPars, LenDat, StockPars, SizeBins=NULL, 
	mod=c("GTG", "LBSPR")) {
  Fleet <- NULL
  Fleet$SL50 <- exp(tryFleetPars[1]) * StockPars$Linf
  Fleet$SL95 <- Fleet$SL50  + (exp(tryFleetPars[2]) * StockPars$Linf)
  Fleet$MLLKnife <- NA
  Fleet$FM <- exp(tryFleetPars[3])
  
  if (mod == "GTG") runMod <-  GTGLBSPRSim(StockPars, Fleet, SizeBins, Nage=(StockPars$AgeMax+1))
  if (mod == "LBSPR") runMod <- LBSPRSim(StockPars, Fleet, SizeBins, Nage=(StockPars$AgeMax+1))
  
  LenDat <- LenDat + 1E-15 # add tiny constant for zero catches
  LenProb <- LenDat/sum(LenDat)
  predProb <- runMod$LCatchFished 
  predProb <- predProb + 1E-15 # add tiny constant for zero catches
  NLL <- -sum(LenDat * log(predProb/LenProb))
  
  # add penalty for SL50 
  trySL50 <- exp(tryFleetPars[1])
  PenVal <- NLL
  Pen <- dbeta(trySL50, shape1=5, shape2=0.01) * PenVal
  if (Pen == 0) Pen <- PenVal * trySL50
  
  # plot(xx, dbeta(xx, shape1=5, shape2=0.01) )
  
  NLL <- NLL+Pen 

  return(NLL)
}