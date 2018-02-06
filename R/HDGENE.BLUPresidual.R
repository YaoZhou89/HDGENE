`HDGENE.BLUPresidual` <- function(Y = NULL,K = NULL,GD,model="AA"){
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
  index.y = !is.na(Y[,2])
  if(is.null(K)){
    print("Building kinship matrix...")
    if (nrow(GD) > 100000) GD = GD[sample(1:100000,replace = F),index.y]
    GD = t(as.matrix(GD))
    Ka = A.mat(GD)
    Kaa = Ka*Ka
    if(model == "full"){
      GDD = 1 - abs(GD - 1)
      Kd = A.mat(GDD)
      Kad = Ka*Kd
      Kdd = Kd*Kd
      K = list(A = Ka, D = Kd, AA = Kaa, AD = Kad, DD = Kdd)
    }else if(model == "AA"){
      K = list(A = Ka, AA = Kaa)
    }else{
      stop("Model error: model could only be \"AA\" or \"full\"!")
    }
  }
  LK = length(K)
  if(LK == 1) stop("K must more than 2!")
  yblup = matrix(NA,nrow(Y),LK+1)
  yblup = as.data.frame(yblup)
  yblup[,1] = Y[,1]
  colnames(yblup)= c(colnames(Y)[1],paste(colnames(Y)[2],"_K",seq(1:LK),sep=""))
  print("Estimating BLUP residual...")
  yhat1 = HDGENE.getBLUP(Y[index.y,2],K);
  for (i in 1:LK){
    yhat2 = HDGENE.getBLUP(Y[index.y,2],K[-i])
    yblup[index.y,i+1] = yhat1 - yhat2;
  }
  return(yblup)
}
