DoOpt <- function(StockPars, LenDat, SizeBins=NULL, mod=c("GTG", "LBSPR")) {
  
  SDLinf <- StockPars$CVLinf * StockPars$Linf
  if (is.null(SizeBins)) {
    SizeBins$Linc <- 5
	SizeBins$ToSize <- StockPars$Linf + StockPars$MaxSD * SDLinf
  }
  if (is.null(SizeBins$ToSize)) 
	SizeBins$ToSize <- StockPars$Linf + StockPars$MaxSD * SDLinf
  
  Linc <- SizeBins$Linc 
  ToSize <- SizeBins$ToSize
 
  LenBins <- seq(from=0, by=Linc, to=ToSize)	
  LenMids <- seq(from=0.5*Linc, by=Linc, length.out=length(LenBins)-1)
  
  sSL50 <- LenMids[which.max(LenDat)]/StockPars$Linf # Starting guesses
  sDel <- 0.2 * LenMids[which.max(LenDat)]/StockPars$Linf
  sFM <- 0.5 
  Start <- log(c(sSL50, sDel, sFM))
  opt <- nlminb(Start, OptFun, LenDat=LenDat, StockPars=StockPars, 
	SizeBins=SizeBins, mod=mod, 
	control= list(iter.max=300, eval.max=400, abs.tol=1E-20))
  
  newFleet <- NULL 
  newFleet$FM <- exp(opt$par[3])
  newFleet$SL50 <- exp(opt$par[1]) * StockPars$Linf 
  newFleet$SL95 <- newFleet$SL50 + exp(opt$par[2]) * StockPars$Linf

  if (mod == "GTG") runMod <-  GTGLBSPRSim(StockPars, newFleet, SizeBins, Nage=(StockPars$AgeMax+1))
  if (mod == "LBSPR") runMod <- LBSPRSim(StockPars, newFleet, SizeBins, Nage=(StockPars$AgeMax+1))
  
  Out <- NULL 
  Out$Ests <- c(FM=newFleet$FM, SL50=newFleet$SL50, SL95=newFleet$SL95, 
	SPR=runMod$SPR)
  Out$PredLen <- runMod$LCatchFished * sum(LenDat)
  return(Out)
}