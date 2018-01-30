HDGENE_A  = function(CV = NULL,Y = NULL, GD = NULL,GM = NULL,maxLoop = 10, file.out = TRUE,getback=FALSE,LD = 0.7,p.threshold = NA , LD.num = 50,GS.prediction = F){
  print("-----------------------Welcome to HDGENE, additive model is working----------------------")
  myEB = EBGWAS(Y = Y,GD = GD,CV= CV,GM = GM, maxLoop = maxLoop,file.output = F, getback = getback,LD = LD,p.threshold = p.threshold, LD.num = LD.num,GS.prediction = GS.prediction)
  beta = myEB$beta
  trait.name = colnames(Y)[2]
  Done = F
  if(GS.prediction) Done = T
  GWAS = myEB$GWAS
  if(is.null(GWAS)){
    print("No significant QTNs,stop at first iteration!")
  }else{
    GWAS.back = GWAS
    index = match(GWAS[GWAS$LOD>2.99,1],GM[,1])
    index.y = which(!is.na(Y[,2]))
    if(length(index)>1){
      myCV = t(as.matrix(GD[index,]))
    }else if(length(index)==1){
      myCV = as.matrix(GD[index,])
    }else{
      Done = TRUE
    }
    CV.back = CV
    while(!Done){
      CV = cbind(CV.back,myCV)
      ## to avoid the colinearity of CV
      CV1 = cbind(1,CV)
      CV1[1,1] = 2 ## to avoid the SD = 0 for others
      CV1 = t(as.matrix(CV1))
      Psort=Blink.LDRemove(Porder=1:nrow(CV1),GDneo=CV1,bound=FASLE,LD=0.99,model="A",orientation="row",LD.num = nrow(CV1))
      CV = t(CV1[Psort,])[,-1]
      myEB = EBGWAS(CV = CV,Y=Y,GD=GD,GM=GM, maxLoop = maxLoop,file.output = F,getback = getback,LD = LD,p.threshold = p.threshold,LD.num = LD.num,GS.prediction=GS.prediction)
      beta = myEB$beta
      GWAS.new = myEB$GWAS
      if(is.null(GWAS.new)){
        Done = TRUE
      }else{
        GWAS.back = GWAS.new
        index1 = match(GWAS.new[GWAS.new$LOD>3,1],GM[,1])
        index.new = union(index,index1)
        if(length(index)==length(index.new)){
          Done = TRUE
        }else{
          colnames(GWAS.new) = colnames(GWAS)
          GWAS = rbind(GWAS,GWAS.new)
          index = index.new
          GWAS = GWAS[!duplicated(GWAS[,1]),]
          index.CV = match(GWAS[,1],GM[,1])
          cat(length(index.CV),"QTNs added as covariates for the next iteration","\n")
          myCV=  t(as.matrix(GD[index.CV,]))
          if(length(setdiff(GWAS[,1],GWAS.back[,1]))==0 & length(setdiff(GWAS.back[,1],GWAS[,1]))==0){
            Done = TRUE
          }
        }
      }
    }
    if(length(setdiff(GWAS[,1],GWAS.back[,1]))==0 & length(setdiff(GWAS.back[,1],GWAS[,1]))==0){
      if(file.out)   write.table(GWAS,paste(trait.name,"_GWAS.txt",sep=""),sep="\t",col.names = T,row.names = F,quote = F)
    }else{
      index = match(GWAS$rs,GM[,1])
      GD1 = t(as.matrix(GD[index,index.y]))
      if(length(index)==1)  GD1 = as.matrix(GD[index,index.y])
      GM1 = GM[index,]
      y = Y[index.y,2]
      y = as.matrix(y)
      if(ncol(GD1)<2){
        cat("Only SNP",index,"is significant!")
        print("No output for GWAS, General Linear Model recommented for this SNP")
        GWAS = NULL
        beta = NULL
      }else{
        print(paste(ncol(GD1),"SNPs were fitted in EM-BLASSO..."))
        if(GS.prediction){
          maxSNP = nrow(Y)/log(nrow(Y))
          if (maxSNP > 40) maxSNP = 40
                               
          if (ncol(GD1) > maxSNP) {
            GD1 = GD1[,1:maxSNP]
            GM1 = GM1[1:maxSNP,]
          }
          print(paste(ncol(GD1),"SNPs were fitted in EM-BLASSO..."))
          myEM = EM_LASSO(CV = CV.back[index.y,],GD = GD1,y = y,GM = GM1)
        }else{
          myEM = EM_LASSO(CV = CV.back[index.y,],GD = GD1,y = y,GM = GM1)
        }
        
        GWAS = myEM$GWAS
        beta = myEM$beta
        if(file.out)   write.table(GWAS,paste(trait.name,"_GWAS.txt",sep=""),sep="\t",col.names = T,row.names = F,quote = F)
      }
    }
  }
  return(list(GWAS=GWAS,beta = beta))
}
