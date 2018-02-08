HDGENE_full = function(CV= NULL,GD=NULL,GDD = NULL,GM=NULL,Y=NULL,seqQTN = NULL,p.threshold = NULL,maxIter=3,GWAS=NULL,LD.pruning = F,position.pruning=NULL,LD=0.7,model=NULL){
  # AA
  print("---------------------Welcome to HDGENE, Additive by Additive is working...-----------------")
  # YY = Y
  index.y = which(!is.na(Y[,2]))
  Y = Y[index.y,]
  CV = CV[index.y,]
  GD = GD[,index.y] - 1
  if(model == "AA"){
    GDD = GD
  }else if (model =="AD"){
    GDD = 1 - abs(GD)
  }else if (model == "DD"){
    GD = 1 - abs(GD)
    GDD = GD
  }else{
    stop("incorrect model!")
  }
  m = nrow(GD)
  # force to do LD pruning when marker number larger than 10e5
  if(m>10e4 & is.null(position.pruning)) LD.pruning = TRUE
  if(LD.pruning & is.null(position.pruning)){
    position.pruning = HDGENE.pruning(GD,GM,10000,LD)
  }
  if(is.null(GWAS)){
    print("---------------------Additive by Additive model-----------------")
    GWAS = HDGENE.pairDetection(CV= CV, GD = GD,GDD = GDD,GM = GM,Y = Y,seqQTN = seqQTN,p.threshold = p.threshold, position.pruning = position.pruning,model=model)
    # write.table(GWAS,paste(colnames(Y)[2],"_GWAS_eps.txt",sep=""),col.names=T,row.names = F,quote=F,sep="\t")
  }else{    GWAS = GWAS
  }
  index.aa = !is.na(GWAS[,2])
  index.aa.save = index.aa
  SNP1.save = NULL
  SNP2.save = NULL
  if(sum(index.aa)==0){
    Done = TRUE
  }else{
    Done = FALSE
  }
  iter = 1
  CV1 = CV
  while(!Done & iter <= maxIter){
    iter = iter +1
    SNP1 = append(SNP1.save,match(GWAS[index.aa,1],GM[,1]))
    SNP2 = append(SNP2.save,match(GWAS[index.aa,2],GM[,1]))
    if(sum(index.aa)==0|iter>maxIter) Done =TRUE
    if(length(SNP2)==length(SNP2.save)) Done = TRUE
    SNP1.save = SNP1
    SNP2.save = SNP2
    GD1 = GD[SNP1,]
    GD2 = GDD[SNP2,]

    GDI = GD1*GD2
    if(length(SNP1)>1) GDI = t(GDI)
    CV = cbind(CV1,GDI)

    if(!Done){
      print("---------------------Additive by Additive model-----------------")
      GWAS = HDGENE.pairDetection(CV= CV, GD = GD,GDD = GDD, GM = GM,Y = Y,seqQTN = seqQTN,p.threshold = p.threshold,position.pruning = position.pruning,model=model)
      # write.table(GWAS,paste(colnames(Y)[2],"_GWAS_eps.txt",sep=""),col.names=T,row.names = F,quote=F,sep="\t")
    }else{
      if(is.null(seqQTN)){
        GDA = NULL
      }else if(length(seqQTN)==1){
        GDA = as.matrix(GD[seqQTN,])
      }else{
        GDA = t(as.matrix(GD[seqQTN,]))
      }
      GDAI = cbind(GDA,GDI)
      GMAI = matrix(NA,ncol(GDAI),3)
      GMAI[,1]=paste("rs",seq(1:ncol(GDAI)),sep="")

      print(paste(ncol(GDAI),"markers are fitted in the EM-BLASSO model ..."))
      if(ncol(GDAI)>1){
        myEM = EM_LASSO(CV=CV1,GD = GDAI,y=Y[,2],GM=GMAI)
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
      colnames(myGWAS) = c("SNP1","SNP2","u","LOD","Variance_explained")
      GWAS = myGWAS
    }
    index.aa = !is.na(GWAS[,2])
  }
  write.table(GWAS,paste(colnames(Y)[2],"_GWAS_eps.txt",sep=""),col.names=T,row.names = F,quote=F,sep="\t")
  return(GWAS = GWAS)
}
