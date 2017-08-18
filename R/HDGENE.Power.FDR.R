`HDGENE.Power.FDR` <- function(GM = NULL, seqQTN = NULL, GWAS = NULL, WS = 1 ){	
#Objects: diagnose the result of GWAS
#Input: ws, the bin size considering as one marker; i
#		info, the output from HDGENE.PheSim
#		GWAS, the 
#Output:
#Author: Yao Zhou
#Last update: Nov 4, 2016
	GWAS[is.na(GWAS[,4]),4] = 1
	m = max(GM[,3])
	GM[,1] = as.numeric(GM[,2]) * m + as.numeric(GM[,3])
	GM.seq = GM[seqQTN,]
	for ( i in 1:(nrow(GM.seq))){
		seq.down = GM.seq[i,1]- WS
		if (seq.down < 1 ) seq.down = 1
		seq.up = GM.seq[i,1] + WS
		index = which(((GM[,1] > seq.down) & ((GM[,1]) < seq.up)))
		
		#GM.seq[i,2] = index[which(GWAS[index,4] == min(GWAS[index,4]))][1]
		p = min(GWAS[index,4])
		seqQTN[i] = index[which(GWAS[index,4] == min(GWAS[index,4]))][1]
		GWAS[index,4] =1
		GWAS[seqQTN[i],4] = p
	}
	
	P.order = order(GWAS[seqQTN,4],na.last = T)
	Porder = order(GWAS[,4],na.last = T)
	result = matrix(0,nrow(GM.seq),3)
	# distance

	
	# P value
	i = 0
	b.f = -1
	delta = 0
	while ( i < nrow(GM.seq)){
		i = i + 1
		a = GWAS[seqQTN[P.order[i]],4] - GWAS[,4]
		b = sum(a > 0)
		if(b == b.f){
			delta = delta + 1
		}else{
			delta = 0
		}
		result[i,1] = (i/nrow(GM.seq))
		result[i,2] = (b + 1 + delta - i)/ (b + 1 + delta)
		result[i,3] = (b + 1 + delta - i)/(nrow(GM)-nrow(GM.seq))
		b.f = b
	}
	return(result = result)
}

