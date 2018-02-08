`HDGENE.cor2`<-function(Y,GD,w=NULL,ms=10000,model = NULL){
  #Objects: calculate R value with covariates
  #Input: pheontype(nx1), ms is marker size for slicing the genotype, genotype(orientation="row", mxn or orientation="col", nxm,) and covariates(nxp)
  #		n is individual number, m is marker number, p is covariate number
  #Output: abs(r)
  #Author: Yao Zhou
  #Last updated: Jun 28, 2016
  if(!is.matrix(Y)) Y=as.matrix(Y)
  n = nrow(Y)
  m = nrow(GD)
  # Orthogonolize phenotype w.r.t. covariates
  {
    if(!is.null(w)){
      w = cbind(1,w)
    }else{
      w = matrix(1,n,1)
    }
    if(!is.matrix(w)) w = as.matrix(w)
    qw = qr(w)
    if( min(abs(diag(qr.R(qw)))) < .Machine$double.eps * m ) {
      stop("Colinear or zero covariates detected");
    }
    w = qr.Q(qw)
    tw=t(w)
    rm(qw)
  }

  # Orthogonolize phenotype w.r.t. covariates

  {
    Y = Y - w%*%crossprod(w,Y)
    colsq = colSums(Y^2)
    div = sqrt(colsq)
    Y = Y/div
    rm(colsq,div)
  }
  #Orthogonolize genotype w.r.t. covariates
  {
    rabs = matrix(NA,nrow = nrow(GD),ncol = nrow(GD))
    m = nrow(GD)
    for(marker in 1:m){
      ntest = nrow(GD)-marker + 1
      ns = ceiling(ntest/ms)
      for(i in 1:ns){
        bottom=(ms*(i-1)+marker)
        # if(bottom<marker) bottom = marker
        if(i<ns){
          up=ms*i+marker-1
        }else{
          up = m
        }

        GDs = GDD[bottom:up,]
        if(bottom==up){
          GDs = as.matrix(GDs)
        }else{
          GDs = t(GDs)
        }
        if(model == "AA"){
          GDs = GD[marker,] * GDs
        }else if (model == "AD"){
          GDs = GD[marker,] * (1 - abs(GDs))
        }else if (model == "DD"){
          GDs = (1- abs(GD[marker,])) * (1 - abs(GDs))
        }else{
          stop("Model type incorrect!")
        }
        GDs = GDs- crossprod(tw,tw%*%GDs)
        colsq= colSums(GDs^2)
        div = sqrt(colsq)
        GDs=t(GDs)/div
        rabs[bottom:up,marker] = abs(GDs%*%Y)
      }
    }

    rm(GDs,div)
  }
  return(rabs)
}

