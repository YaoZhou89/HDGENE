`EM_LASSO` <- function(CV,GD,y,GM){
  # code was part of mrMLM with minor revison
  Lasso_EM<-function(x,z,y)
  {
    print("EM algorithm is working...")
    n<-nrow(z);k<-ncol(z)
    b<-solve(crossprod(x,x))%*%(crossprod(x,y))
    v0<-as.numeric(crossprod((y-x%*%b),(y-x%*%b))/n)
    u<-matrix(rep(0,k),k,1)
    v<-matrix(rep(0,k),k,1)
    s<-matrix(rep(0,k),k,1)
    for(i in 1:k)
    {
      zz<-z[,i]
      s[i]<-((crossprod(zz,zz)+1e-100)^(-1))*v0
      u[i]<-s[i]*crossprod(zz,(y-x%*%b))/v0
      v[i]<-u[i]^2+s[i]
    }
    vv<-matrix(rep(0,n*n),n,n);
    for(i in 1:k)
    {
      zz<-z[,i]
      vv=vv+tcrossprod(zz,zz)*v[i]
    }
    vv<-vv+diag(n)*v0
    iter<-0;err<-1000;iter_max<-5000;err_max<-1e-10
    tau<-0;omega<-0
    while((iter<iter_max)&&(err>err_max))
    {
      iter<-iter+1
      v01<-v0
      v1<-v
      b1<-b
      vi<-solve(vv)
      xtv<-crossprod(x,vi)
      if(ncol(x)==1)
      {
        b<-((xtv%*%x)^(-1))*(xtv%*%y)
      }else
      {
        b<-solve(xtv%*%x)%*%(xtv%*%y)
      }
      r<-y-x%*%b
      ss<-matrix(rep(0,n),n,1)
      for(i in 1:k)
      {
        zz<-z[,i]
        zztvi<-crossprod(zz,vi)
        u[i]<-v[i]*zztvi%*%r
        s[i]<-v[i]*(1-zztvi%*%zz*v[i])
        v[i]<-(u[i]^2+s[i]+omega)/(tau+3)
        ss<-ss+zz*u[i]
      }
      v0<-as.numeric(crossprod(r,(r-ss))/n)
      vv<-matrix(rep(0,n*n),n,n)
      for(i in 1:k)
      {
        zz<-z[,i]
        vv<-vv+tcrossprod(zz,zz)*v[i]
      }
      vv<-vv+diag(n)*v0
      err<-(crossprod((b1-b),(b1-b))+(v01-v0)^2+crossprod((v1-v),(v1-v)))/(2+k)
      beta<-t(b)
      sigma2<-v0
    }
    print(paste("EM converged at ",iter,"th iteration!"))
    return (list(u=u,sigma2=sigma2,beta=beta))
  }
  multivanormal<-function(y,mean,sigma)
  {
    pdf_value<-(1/sqrt(2*3.14159265358979323846*sigma))*exp(-(y-mean)*(y-mean)/(2*sigma));
    return (pdf_value)
  }
  likelihood<-function(xxn,xxx,yn,bbo)
  {
    nq<-ncol(xxx)
    ns<-nrow(yn)
    at1<-0
    ww1<-as.matrix(which(abs(bbo)>1e-5))
    at1<-dim(ww1)[1]
    lod<-matrix(rep(0,nq),nq,1)
    if(at1>0.5){
      ad<-cbind(xxn,xxx[,ww1])
    }else{
      ad<-xxn
    }
    if(abs(min(eigen(crossprod(ad,ad))$values))<1e-6){
      bb<-solve(crossprod(ad,ad)+diag(ncol(ad))*0.01)%*%crossprod(ad,yn)
    }else{
      bb<-solve(crossprod(ad,ad))%*%crossprod(ad,yn)
    }
    vv1<-as.numeric(crossprod((yn-ad%*%bb),(yn-ad%*%bb))/ns)
    ll1<-sum(log(abs(multivanormal(yn,ad%*%bb,vv1))))
    sub<-1:ncol(ad)
    if(at1>0.5)
    {
      for(i in 1:at1)
      {
        ij<-which(sub!=sub[i+ncol(xxn)])
        ad1<-ad[,ij]
        if(abs(min(eigen(crossprod(ad1,ad1))$values))<1e-6){
          bb1<-solve(crossprod(ad1,ad1)+diag(ncol(ad1))*0.01)%*%crossprod(ad1,yn)
        }else{
          bb1<-solve(crossprod(ad1,ad1))%*%crossprod(ad1,yn)
        }
        vv0<-as.numeric(crossprod((yn-ad1%*%bb1),(yn-ad1%*%bb1))/ns);
        ll0<-sum(log(abs(multivanormal(yn,ad1%*%bb1,vv0))))
        lod[ww1[i]]<--2.0*(ll0-ll1)/(2.0*log(10))
      }
    }
    return (lod)
  }
  if (is.null(CV)==FALSE)
  {
    z<-cbind(matrix(1,nrow(GD),1),CV)
    z = as.matrix(z)
  }else{
    z<-matrix(1,nrow(GD),1)
  }
  GD = as.matrix(GD)
  y = as.matrix(y)
  u1<-Lasso_EM(z,GD,y)
  result1<-u1$u
  Res1 = as.vector(result1)
  le2<-length(na.omit(which(abs(Res1)>1e-5)))
  # print("le2:")
  # print(le2)
  if(le2==0){
    print("There is no QTN significant")
    GWAS = NULL
    beta = NULL
  }else{
    sig1<-which(abs(Res1)>1e-5)
    # print("sig1:")
    # print(sig1)
    bbo<-matrix(0,le2,1)
    for (i in 1:le2){
      bbo[i,]=Res1[sig1[i]]
    }
    her<-vector(length=le2)
    # AA = max(GD[,1:2])
    for (i in 1:le2){
      p1<-length(as.vector(which(GD[,sig1[i]]==1)))/length(GD[,sig1[i]])
      p2<-1-p1
      her[i]=(((p1+p2)-(p1-p2)^2)*(Res1[sig1[i]])^2)/var(y)
    }
    xxxx<-GD[,sig1]
    xxxx = as.matrix(xxxx)
    yn<-as.matrix(y)
    xxn<-z
    lod<-likelihood(xxn,xxxx,yn,bbo)
    beta=u1$beta
    GWAS = cbind(GM[sig1,],Res1[sig1],lod,her)
    colnames(GWAS)[-(1:ncol(GM))]=c("u","LOD","Variance_explained")
  }

  return(list(GWAS=GWAS,beta=beta))
}



