`HDGENE.putback`<-function(seqQTN = NULL,GD = NULL,GM=NULL, LD = 0.7){
  # Objects: get all SNPs that are correlated with seqQTN
  # Input:
  #   seqQTN: psudo QTN detected by BLINK
  #   GD: genotype, big.matrix type, m by n. m: # of marker, n: # of individual
  #   LD: same as in the LD remove part
  # Output:
  #   GDL: genotype for Bayesian LASSO model
  # Author: Yao Zhou
  # Last update: 5/3/2017

  GDi = GD[seqQTN,]
  GDi = t(as.matrix(as.data.frame(GDi)))
 # index.can= apply(GD1[1:1000,],1,function(x) sum(cor(GDi,x)>0.7)>0)
  GD1 = as.matrix(as.data.frame(GD))
  c = cor(GDi,t(GD1))
  index.can = apply(c,2,function(x) sum(x>LD)>0)
  index.can = which(index.can)
  GDL = GD[index.can,]
  GDL = t(as.matrix(as.data.frame(GDL)))
  GML = GM[index.can,]
  rm(GDi,c,GD1)
  return(list(GD=GDL,GM=GML))
}



