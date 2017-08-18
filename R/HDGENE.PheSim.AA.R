`HDGENE.PheSim.AA` <- function(GD=NULL,H2=0.75,alpha=.95,NQTN.A=10,NQTN.AA = 5,IN = 5,distribution="norm"){
  #Objects: simulate Additive and Dominance effects
  #Input:
  #	GD: genotype, nxm, with taxa 
  #	h2: heritability
  #	NQTN.A: QTN number with additive effect
  # NQTN.AA:  number QTNs interact with the one QTN
  # IN: interaction number, number of QTNs with interaction
  #	d: the effect of Dominance. d=0 : No dominance effect; 0<d<1 : incomplete dominance effect; d=1 : complete dominance effect; d >1: overdominance effect
  #	dr: the ratio of dominance effect, dr=1: all QTNs have dominance effect; dr=0: None of QTNs has dominance effect
  #	AddEff: if FALSE, only simulate dominance effect
  #   position: "same" or not, if "same", sampling Dominance position from additive position, or sampling from different position
  #   orientation: "col" or "row"
  #Output:
  #	addeffect: additive effect
  # 	effect: genetic effect=addeffect + epistatic effects
  #	y: simulated phenotype (nx1)
  #	QTN.position: position of QTNs (or additive if position!="same")
  #	QTN.position.AA: position for epistatic QTNs (pair-wise)
  # Author: Yao Zhou
  # Created date: May 31, 2017
  # Last update: May 31, 2017
  n = ncol(GD)
  m = nrow(GD)
  QTN.position.A = sample(m,NQTN.A,replace=F)
  GDA = t(GD[QTN.position.A,])
  AAinfo = matrix(NA,IN*NQTN.AA,4)
  for (i in 1:IN){
    QTN.position.AA = sample(setdiff(seq(1:m),QTN.position.A),NQTN.AA,replace=F)
    for(j in 1:NQTN.AA){
      AAinfo[(i-1)*NQTN.AA+j,1] = QTN.position.A[i]
      AAinfo[(i-1)*NQTN.AA+j,2] = QTN.position.AA[j]
    }
    GDAA = t(GD[QTN.position.AA,])
    if(i==1){
      AA = GDA[,i]*GDAA
    }else{
      AA = cbind(AA,GDA[,i]*GDAA)
    }
  }
  colnames(AAinfo) = c("QTN1","QTN2","Effect","Variance")
  #QTN effects
  if(distribution=="norm"){
    addeffect = rnorm(NQTN.A,0,1.5)
    aaeffect = rnorm(ncol(AA),0,1)
  }else if (distribution=="geometry"){
    addeffect = alpha^(1:NQTN.A)
    aaeffect = alpha^(1:ncol(AA))
  }else if(distribution=="even"){
    addeffect=rep(1,NQTN.A)
    aaeffect=rep(1,ncol(AA))
  }
  AAinfo[,3] = aaeffect
  #Simulate phenotype
  
  sum.table = matrix(0, NQTN.A, 3)
  sum.table[,1] = QTN.position.A	
  Ae = GDA %*% addeffect
  for (i in 1: NQTN.A){
    sum.table[i,3] = var(as.matrix(GDA[,i]) %*% addeffect[i])	
  }
  colnames(sum.table)=c("QTNs","Effect","Variance")
  sum.table[,2] = addeffect
  VA = var(Ae)
  AAe = AA %*% aaeffect	
  for( i in 1:ncol(AA)){
    AAinfo[i,4] = var(as.matrix(AA[,i]) %*% aaeffect[i])
  }
  VI = var(AAe)
  
  blup = Ae + AAe

  effectvar = VA + VI
  residualvar = (effectvar-H2*effectvar)/H2
  residual = rnorm(n,0,sqrt(residualvar))
  phe = as.matrix(blup+residual)
  y = cbind(NA,phe)
  return(list(add = sum.table, eps = AAinfo, VA = VA,VI = VI,VG = effectvar,ve = residualvar,y = y))
}

