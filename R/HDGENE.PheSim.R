`HDGENE.PheSim` <- function(GD=NULL,h2=.75,alpha=.95,NQTN=10,AddEff=TRUE,position="same",distribution="norm",d = 1, dr = 1,orientation="row"){
#Objects: simulate Additive and Dominance effects
#Input:
#	GD: genotype, nxm, with taxa 
#	h2: heritability
#	NQTN: QTN number
#	d: the effect of Dominance. d=0 : No dominance effect; 0<d<1 : incomplete dominance effect; d=1 : complete dominance effect; d >1: overdominance effect
#	dr: the ratio of dominance effect, dr=1: all QTNs have dominance effect; dr=0: None of QTNs has dominance effect
#	AddEff: if FALSE, only simulate dominance effect
#   position: "same" or not, if "same", sampling Dominance position from additive position, or sampling from different position
#   orientation: "col" or "row"
#Output:
#	addeffect: additive effect
# 	effect: genetic effect=addeffect + domeffect
#	y: simulated phenotype (nx1)
#	QTN.position: position of QTNs (or additive if position!="same")
#	QTN.position.D: position of dominance QTNs
#Author: Yao Zhou
#Created date: May 2nd, 2016
#Last update: Oct 31, 2016

	n=ncol(GD)
	m=nrow(GD)
#Sampling QTN
	QTN.position=sample(m,NQTN,replace=F)
	NumD=floor(NQTN*dr)
	if(position=="same"){
		QTN.position.D=QTN.position[1:NumD]
	}else{
		QTN.position.D=sample(m,NumD,replace=F)
	}
	if(NumD==0) QTN.position.D=NULL
	SNPQ = t(as.matrix(GD[QTN.position,]))
	SNPQ.D = 1-abs(t(as.matrix(GD[QTN.position.D,]))-1)	
	
#QTN effects
	if(distribution=="norm"){
		addeffect=rnorm(NQTN,0,1)
		if(dr!=0){
			domeffect=addeffect[1:NumD]*d
		}
	}else if (distribution=="geometry"){
		addeffect=alpha^(1:NQTN)
		if(dr!=0){
			domeffect=addeffect[1:NumD]*d
		}
	}else if(distribution=="even"){
		addeffect=rep(1,NQTN)
		if(dr!=0){
			domeffect=addeffect[1:NumD]*d
		}
	}
#Simulate phenotype
	VA = 0;
	VD = 0;
	A = 0;
	D = 0;
	sum.table = NULL
	lenDom = length(QTN.position.D)
	sum.table = matrix(0, NQTN, 5)
	sum.table[,1] = QTN.position	
	sum.table[1:lenDom,2] = QTN.position.D
	if(AddEff){	
		A = SNPQ %*% addeffect
		for (i in 1: NQTN){
			sum.table[i,3] = var(as.matrix(SNPQ[,i]) %*% addeffect[i])	
		}
		#VA = sum(sum.table[,3])
		VA = var(A)
		if(dr>0){
			D = SNPQ.D%*%domeffect	
			for( i in 1:lenDom){
				sum.table[i,4] = var(as.matrix(SNPQ.D[,i]) %*% domeffect[i])
			}
			#VD = sum(sum.table[,4])
			VD = var(D)
		}
	}else{
		D = SNPQ.D%*%domeffect
		for( i in 1:lenDom){
			sum.table[i,4] = var(as.matrix(SNPQ.D[,i]) %*% domeffect[i])
		}
		#VD = sum(sum.table[,4])
		VD = var(D)
	}
	
	effect = A + D
	sum.table[,5] = sum.table[,3] + sum.table[,4]	
	
	effectvar = VA + VD
	residualvar = (effectvar-h2*effectvar)/h2
	residual = rnorm(n,0,sqrt(residualvar))
	phe = as.matrix(effect+residual)
	y = cbind(NA,phe)
	sum.table = as.data.frame(sum.table)
	colnames(sum.table) = c("QTN.position.A","QTN.position.D", "VA", "VD", "VG")
	return(list(summary = sum.table, VE = residualvar, VA = VA, VD = VD, VG = effectvar, y = y, effect = effect, residual = residual, QTN.position = QTN.position, QTN.position.D = QTN.position.D))
}

