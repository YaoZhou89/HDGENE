`HDGENE_AA_full` <- function(CV = NULL,GD = NULL,GM = NULL,Y = NULL,p.threshold = NULL,seqQTN = NULL,position.pruning=NULL,LD = 0.7){
  # Objects: doing G by G using GLM, full screen
  # inputs:
  #		CV: n by p , p fixed effects
  #		GD: m by n , big matrix
  #		Y: n by 2 matrix
  #   seqQTN: significant QTNs from Additive model
  #   GM: m by 3 matrix
  #   p.threshold: threshold for calculating interaction
  # Outputs:
  #		GWAS: 5 columns, (SNP1, SNP2, Pvalue_1, Pvalue_2, Pvalue_1*2)
  # Author: Yao Zhou
  # Last Update: July 19, 2017

  ## prepare input
  index.y = which(!is.na(Y[,2]))
  if(length(seqQTN)==0){
    NeedDo = FALSE
    myGWAS = NULL
    print("No marginal effect!")
    GDA = NULL
    CV.BIC = CV
    CV.BIC = CV.BIC[index.y,]
  }else{
    if(length(seqQTN)==1){
      GDA = matrix(as.numeric(GD[seqQTN,]),length(GD[seqQTN,]),1)
    }else{
      GDA = t(as.matrix(GD[seqQTN,]))
    }
    CV.BIC = cbind(CV,GDA)
    NeedDo = TRUE
    CV = CV[index.y,]
    CV.BIC = CV.BIC[index.y,]
  }

  if(!is.null(position.pruning)){
    GD1 = deepcopy(GD, rows = position.pruning)
    GD = GD1
    GM = GM[position.pruning,]
  }

  GD1 = deepcopy(GD,cols=index.y)
  n = length(index.y)
  m = nrow(GM)

  if(is.null(p.threshold)) p.threshold = 0.01/nrow(GM)
  if(is.null(CV.BIC)){
    k = 3
  }else{
    if(is.null(ncol(CV.BIC))|length(seqQTN)==1){
      k = 4
    }else {
      k = ncol(CV.BIC) + 3 ## degree of freedom
    }
  }
  Eps = FALSE
  index.AA = NULL
  indicator.s = 0
  iter = 0
  CV1 = CV
  t1 = proc.time()
  print("Testing interaction....")
  if(m > 10e4){
    aa_cor = HDGENE.cor.new(Y = Y[index.y,2], GD = GD1, w = CV.BIC, ms = 100000)
    aa_order = order(aa_cor,decreasing = T,na.last = T)
    aa_length = sum(!is.na(aa_cor))
    aa_log = n/log(n)
    aa_lim = min(aa_log,aa_length)
    index.aa = aa_order[1:aa_lim]
  }else{
    aa_cor = HDGENE.cor(Y = Y[index.y,2], GD = GD1, w = CV.BIC, ms = 100000)
    p = Blink.rtop(r=aa_cor,df = n - k)
    index.aa = which(p < p.threshold)
  }
  t2 = proc.time()
  t3 = t2 - t1
  cat("Interaction testing finished in",round(as.numeric(t3)[3],2), "....","\n")

  ## transform correlation to pvalue
  cat(length(index.aa),"interactions passed the p value threshold...","\n")
  if(length(index.aa)>0){
    cor_order = index.aa[order(aa_cor[index.aa],decreasing = T,na.last = T)]
    SNP2 = cor_order%%m
    SNP2[SNP2==0]=m
    SNP1 = trunc(cor_order/m)
    SNP1[SNP2!=m] = SNP1[SNP2!=m] +1 ## max correlation is at aa_cor[SNP2[1],SNP1[1]]
    index.snp1 = SNP1
    index.snp2 = SNP2
    GD2 = GD[index.snp2,]
    if(length(index.snp2)==1){
      GD2 = as.matrix(GD2)
      GD1 = as.matrix(GD[index.snp1,])
    }else{
      GD2 = t(as.matrix(GD2))
      GD1 = t(as.matrix(GD[index.snp1,]))
    }
    GDI = GD1 * GD2
    cat(ncol(GDI),"QTNs for LDRemove","\n")
    if(ncol(GDI)>1){
      seqQTN.index = Blink.LDRemove(GDneo=t(GDI),Porder = seq(1:ncol(GDI)),LD = LD,orientation = "row",LD.num = 100)
      cat(length(seqQTN.index),"QTNs for BIC selection","\n")
      GDI = GDI[,seqQTN.index]
      SNP1 = SNP1[seqQTN.index]
      SNP2 = SNP2[seqQTN.index]
      if(length(seqQTN.index)>1){
        myBIC = Blink.BICselection(CV = CV.BIC, Psort=seq(1:length(seqQTN.index)),GD=GDI[index.y,],Y=Y[index.y,],orientation="col",BIC.method="naive")
        seqQTN.index = myBIC$seqQTN
        cat(length(seqQTN.index),"QTNs detected with interactions!","\n")
        GDI = GDI[,seqQTN.index]
        SNP1 = SNP1[seqQTN.index]
        SNP2 = SNP2[seqQTN.index]
      }
    }
    GDAI = cbind(GDA,GDI)
    GMAI = matrix(NA,ncol(GDAI),3)
    GMAI[,1]=paste("rs",seq(1:ncol(GDAI)),sep="")
    print(paste(ncol(GDAI),"markers are fitted in the EM-BLASSO model ..."))
    if(ncol(GDAI)>1){
      myEM = EM_LASSO(CV=CV1,GD = GDAI[index.y,],y=Y[index.y,2],GM=GMAI)
      GWAS = myEM$GWAS
      index.s = match(myEM$GWAS[,1],GMAI[,1])
    }else{
      index.s = c(1)
      GWAS = data.frame(rs=c("rs1"),chr=c(1),pos=c(1),u=c(1),LOD=c(4),Variance=c(10))
    }


    SNPa = index.s[which(index.s<=length(seqQTN))]
    SNPaa = setdiff(index.s,SNPa)
    if(is.null(ncol(GDA))){
      NGDA = 0
    }else{
      NGDA = ncol(GDA)
    }
    if(length(SNPaa)>0){
      SNPaa1 = SNP1[SNPaa-NGDA]
      SNPaa2 = SNP2[SNPaa-NGDA]
    }
    myGWAS = matrix(NA,length(index.s),5)
    myGWAS = as.data.frame(myGWAS)
    if(length(SNPa)>0){
      SNPa.c = seqQTN[SNPa]
      myGWAS[1:length(SNPa),1] = as.character(GM[SNPa.c,1])
    }
    if(length(SNPaa)>0){
      myGWAS[(length(SNPa)+1):nrow(myGWAS),1] = as.character(GM[SNPaa1,1])
      myGWAS[(length(SNPa)+1):nrow(myGWAS),2] = as.character(GM[SNPaa2,1])
    }

    myGWAS[,3] = as.numeric(GWAS[,4])
    myGWAS[,4] = as.numeric(GWAS[,5])
    myGWAS[,5] = as.numeric(GWAS[,6])
    colnames(myGWAS) = c("SNP1","SNP2","u","LOD","Rsquare%")

  }else{
    print("No epistatic effects detected!")
    myGWAS = NULL
  }
  ##
  return(GWAS = myGWAS)
}



