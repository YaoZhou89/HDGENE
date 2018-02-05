# HDGENE
HDGENE: High dimensional genetic association analysis package

# Installation
    library(devtools)

    install_github("YaoZhou89/HDGENE")

# Running the tests
    library(HDGENE)
##read data
       
    myGD = read.big.matrix("demo.dat",head=F,sep="\t",type="char")

    myGM = read.table("demo.map",head=T)

    myY = read.table("demo.txt",head=T)

    myCV = read.table("demo.cov",head=T)

##run additive model

    myHDGENE = HDGENE_A(Y=myY,GD=myGD,GM=myGM,CV=myCV, maxLoop = 10,file.out = T)

# Authors
Dr. Yao Zhou (yao.zhou@genetics.ac.cn)
