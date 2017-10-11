`FarmCPU.LM` <-
#cmpfun(
function(y,w=NULL,GDP,orientation="col",model="A",ncpus=2,myModel=NULL,seqQTN=NULL,npc=0){
    #Object: 1. To quickly sovel LM with one variable substitute multiple times
    #Object: 2. To fit additive and additive+dominace model
    #intput: y - dependent variable
    #intput: w - independent variable
    #intput: GDP - independent variable of substitution (GDP)
    #intput: model - genetic effects. Options are "A" and "AD"
    #Output: estimate, tvalue, stderr and pvalue ( plus the P value of F test on both A and D)
    #Straitegy: 1. Separate constant covariates (w) and dynamic coveriates (x)
    #Straitegy: 2. Build non-x related only once
    #Straitegy: 3. Use apply to iterate x
    #Straitegy: 4. Derive dominance indicate d from additive indicate (x) mathmaticaly
    #Straitegy: 5. When d is not estimable, continue to test x
    #Authors: Xiaolei Liu and Zhiwu Zhang
    #Start  date: March 1, 2013
    #Last update: March 6, 2013
    ##############################################################################################
    #print("FarmCPU.LM started")
    #print(date())
    #print(paste("No. Obs: ",length(y),sep=""))
    #print("diminsion of covariates and markers")
    if(!is.null(w))#print(dim(w))

    #print("Memory used at begining of LM")
    #print(memory.size())
    gc()
    #Constant section (non individual marker specific)
    #---------------------------------------------------------
    #Configration
    nd=20 #number of markes for checking A and D dependency
    threshold=.99 # not solving d if correlation between a and d is above this
    N=length(y) #Total number of taxa, including missing ones
    direction=2
    if(orientation=="row")direction=1
    #print("direction")
    #print(direction)
    #Handler of non numerical y a and w

    if(!is.null(w)){
        nf=length(w)/N
        w=matrix(as.numeric(as.matrix(w)),N,nf  )
        w=cbind(rep(1,N),w)#add overall mean indicator
        q0=ncol(w) #Number of fixed effect excluding gnetic effects
    }else{
        w=rep(1,N)
        nf=0
        q0=1
    }

    y=matrix(as.numeric(as.matrix(y)),N,1  )

    #print("Adding overall mean")
    #print(date())
    #print("Build the static section")
    #print(date())

    #n=nrow(w) #number of taxa without missing
    n=N
    if(nd>n)nd=n #handler of samples less than nd
    k=1 #number of genetic effect: 1 and 2 for A and AD respectively
    if(model=="AD")k=2

    q1=(q0+1) # vecter index for the posistion of genetic effect (a)
    q2=(q0+1):(q0+2) # vecter index for the posistion of genetic effect (a and d)
    df=n-q0-k #residual df (this should be varied based on validating d)

    iXX=matrix(0,q0+k,q0+k) #Reserve the maximum size of inverse of LHS
    #theNA=c(rep(NA,q0),rep(0,k)) # this should not be useful anymore

    ww=crossprod(w,w)
    wy=crossprod(w,y)
    yy=crossprod(y,y)
    wwi=ginv(ww)

    #print("Prediction")
    #print(date())

    #Statistics on the reduced model without marker
    rhs=wy
    beta <- crossprod(wwi,rhs)
    ve=(yy-crossprod(beta,rhs))/df
    se=sqrt(diag(wwi)*ve)
    tvalue=beta/se
    pvalue <- 2 * pt(abs(tvalue), df,lower.tail = FALSE)
    P0=c(beta[-1],tvalue[-1],se[-1],pvalue[-1])
    yp=w%*%beta

    if(npc!=0){
        betapc = beta[2:(npc+1)]
        betapred = beta[-c(1:(npc+1))]
    }else{
        betapc = NULL
        betapred = beta[-1]
    }
    #print("Detecting genotype coding system")
    #print(date())

    #Finding the middle of genotype coding (1 for 0/1/2 and 0 for -1/0/1)
    s=5 # number of taxa sampled
    t0=which(GDP[1:s,]<0)
    t1=which(GDP[1:s,]>1)
    middle=0
    if(length(t0)<length(t1)) middle=1

    #print("Memory used after setting LM")
    #print(memory.size())
    gc()
    #Dynamic section on individual marker
    #print("Iterating.................")
    #print(date())
    #print("dimension of GD")
    #print(dim(x))
    #sfInit(parallel=ncpus>1, cpus=ncpus)
    ##print(sprintf('%s cpus are used', sfCpus()))

    #---------------------------------------------------------
    #P <- matrix(NA,nrow=nrow(GDP),ncol=4*(nf+1))
    if(orientation=="row"){
        P <- matrix(NA,nrow=nrow(GDP),ncol=nf+1)
        ntest=nrow(GDP)
    }else{
        P <- matrix(NA,nrow=ncol(GDP),ncol=nf+1)
        ntest=ncol(GDP)
    }

    if(orientation=="row"){
        B <- matrix(NA,nrow=nrow(GDP),ncol=1)
    }else{
        B <- matrix(NA,nrow=ncol(GDP),ncol=1)
    }

    for(i in 1:ntest){
        if(orientation=="row"){
            x=GDP[i,]
        }else{
            x=GDP[,i]
        }

        #P <- apply(x,direction,function(x){
        #P <- sfApply(x,direction,function(x){
        r=1 #initial creteria for correlation between a and d
        if(model=="AD"){
            d=1-abs(x-middle)
            r=abs(cor(x[1:nd],d[1:nd]))
            if(is.na(r))r=1
            if(r<=threshold) x=cbind(x,d) # having both a and d as marker effects
        }

        #Process the edge (marker effects)
        xy=crossprod(x,y)
        xx=crossprod(x,x)

        if(model=="AD"&r<=threshold){
            xw=crossprod(x,w)
            wx=crossprod(w,x)
            iXX22 <- solve(xx-xw%*%wwi%*%wx)
            iXX12 <- (-wwi)%*%wx%*%iXX22
            iXX21 <- (-iXX22)%*%xw%*%wwi
            iXX11 <- wwi + wwi%*%wx%*%iXX22%*%xw%*%wwi
        }else{
            xw=crossprod(w,x)
            B21 <- crossprod(xw, wwi)
            t2=B21%*%xw #I have problem of using crossprod and tcrossprod here
            B22 <- xx - t2
            invB22=1/B22
            NeginvB22B21 <- crossprod(-invB22,B21)
            iXX11 <- wwi + as.numeric(invB22)*crossprod(B21,B21)
        }

        #Derive inverse of LHS with partationed matrix
        iXX[1:q0,1:q0]=iXX11

        if(r>threshold){
            iXX[q1,q1]=invB22
            iXX[q1,1:q0]=NeginvB22B21
            iXX[1:q0,q1]=NeginvB22B21
        }else{
            iXX[q2,q2]=iXX22
            iXX[q2,1:q0]=iXX21
            iXX[1:q0,q2]=iXX12
        }

        #statistics
        rhs=c(wy,xy) #the size varied automaticly by A/AD model and validated d

        if(abs(r)>threshold & model=="AD"){
            beta <- crossprod(iXX[-(q0+k),-(q0+k)],rhs) #the last one (d) dose not count
            df=n-q0-1
        }else{
            beta <- crossprod(iXX,rhs)   #both a and d go in
            df=n-q0-2
        }
        if(model=="A") df=n-q0-1 #change it back for model A

        ve=(yy-crossprod(beta,rhs))/df #this is a scaler

        #using iXX in the same as above to derive se
        if(abs(r)>threshold & model=="AD"){
            se=sqrt(diag(iXX[-(q0+k),-(q0+k)])*ve)
        }else{
            se=sqrt(diag(iXX)*ve)
        }

        tvalue=beta/se
        pvalue <- 2 * pt(abs(tvalue), df,lower.tail = FALSE)

        #Handler of dependency between  marker are covariate
        if(!is.na(abs(B22[1,1]))){
            if(abs(B22[1,1])<10e-8)pvalue[]=NA}

        #Calculate P value for A+D effect
        if(model=="AD"){
            #the last bit could be d or a, the second last may be marker effect not even not
            #In either case, calculate F and P value and correct them later
            markerbits=(length(beta)-1):length(beta)
            SSM=crossprod(beta[markerbits],rhs[markerbits])
            F=(SSM/2)/ve
            PF=df(F,2,df)

            #correcting PF with P from t value
            if(r>threshold) PF=pvalue[length(pvalue)]
        }

        #in case AD model and a/d dependent, add NA column at end
        if(r>threshold & model=="AD"){
            beta=c(beta,NA)
            tvalue=c(tvalue,NA)
            se=c(se,NA)
            pvalue=c(pvalue,NA)
        }

        if(model=="AD"){
            result=c(beta[-1],tvalue[-1],se[-1],pvalue[-1],PF)
        }else{
            #result=c(beta[-1],tvalue[-1],se[-1],pvalue[-1])
            #P[i,]=c(beta[-1],tvalue[-1],se[-1],pvalue[-1])
            P[i,c(1:(nf+1))]=pvalue[-1]
            B[i,]=beta[length(beta)]
            #P[i,c(1:(nf+1))]=beta[-1]
            #P[i,c((nf+2):(2*nf+2))]=pvalue[-1]
            #P[i,c((nf+2):(2*nf+2))]=tvalue[-1]
            #P[i,c((2*nf+3):(3*nf+3))]=se[-1]
            #P[i,c((3*nf+4):(4*nf+4))]=pvalue[-1]
        }
    }
    #}
    #}) #end of defyning apply function
    #sfStop()

    #print("iteration accoplished")
    #print(date())
    #print("Memory used after iteration")
    #print(memory.size())
    gc()

    #Final report
    #---------------------------------------------------------
    #P=t(as.matrix(P))
    #P=as.matrix(P)

    PF=P[,ncol(P)]
    if(model=="AD")P=P[,-ncol(P)]

    #print("FarmCPU.LM accoplished")
    #print(date())

    #print(dim(P))
    #print(P[1:5,])
    #print("Memory used at end of LM")
    #print(memory.size())
    gc()
    #print(head(P))
    return(list(P=P,P0=P0,PF=PF,Pred=yp,betapc=betapc,betapred=betapred,B=B))
} #end of function(
#)#end of cmpfun(
