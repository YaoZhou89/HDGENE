# HDGENE
HDGENE: High dimensional genetic association analysis package

# Installation
library(devtools)

install_github("YaoZhou89/HDGENE")

# Running the tests
library(HDGENE)

##read data
setwd("~/data/emaize")
myGD = read.big.matrix("emaize_220K.dat",head=F,sep="\t",type="char")
myGM = read.table("emaize_220K.map",head=T)
myY = read.table("pheno_emaize.txt",head=T)
myCV = read.table("emaize_PCA3.txt",head=T)

##run additive model
myHDGENE = HDGENE_A(Y=myY,GD=myGD,GM=myGM,CV=myCV, maxLoop = 10,file.out = T)

# Authors
Yao Zhou (yao_wsu.zhou@wsu.edu)
