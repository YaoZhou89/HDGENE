`HDGENE.getBLUP` <- function(Y,K){
  # Objects: get the BLUP residual
  # Input:
  #     K is a list, and length could be 1;
  #     GD: genotype data, big.matrix
  #     model: AA, only calculate A by A interaction
  #     model: full, calculate all epistatic interactions: AA,AD,DA,DD
  #
  # Output: BLUPresidual
  # Last updated: 02-05-2018
  # Authors: Yao Zhou
  KL = length(K)
  ETA = list()
  for (i in 1:KL){
    ETA[[i]] = list(Ki = K[i][[1]],model="RKHS")
  }

  dirFile = paste(getwd(),"/HDGENE_ploygenetic_BGLR_",sep="")
  fmA <- BGLR(y=Y,
              ETA=ETA,
              nIter=100000,
              burnIn=10000,
              saveAt=dirFile,
              verbose=FALSE)
  yhat = fmA$yHat
  return(yhat)
}