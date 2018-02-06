# HDGENE
HDGENE: High dimensional genetic association analysis package

# Installation
    library(devtools)

    install_github("YaoZhou89/HDGENE")

# Running the tests
    library(HDGENE)
## # read data

    myGD = read.big.matrix("demo.dat",head=F,sep="\t",type="char")

    myGM = read.table("demo.map",head=T)

    myY = read.table("demo.txt",head=T)

    myCV = read.table("demo.cov",head=T)

## # run additive only model

    myHDGENE = HDGENE_A(Y=myY,GD=myGD,GM=myGM,CV=myCV, maxLoop = 10,file.out = T)

## # run epistatical models:
###  1. calcuating the blup residuals: (model: AA: A by A only , full: AA, AD, DA and DD)

    yblup = HDGENE.BLUPresidual(Y = myY[,c(1,2)],GD = myGD,model = "AA")
### 2. run additive and interaction models

    myHDGENE_EPI = HDGENE_EPI(Y= yblup,GD=myGD,GM=myGM,CV=myCV, maxLoop = 10,file.out = T,model = "AA")
    
# Authors
Dr. Yao Zhou (yao.zhou@genetics.ac.cn)
