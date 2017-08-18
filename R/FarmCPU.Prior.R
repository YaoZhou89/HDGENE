`FarmCPU.Prior` <-
function(GM,P=NULL,Prior=NULL,kinship.algorithm="FARM-CPU"){
    #Object: Set prior on existing p value
    #Input: GM - m by 3  matrix for SNP name, chromosome and BP
    #Input: Prior - s by 4  matrix for SNP name, chromosome, BP and Pvalue
    #Input: P - m by 1 matrix containing probability
    #Requirement: P and GM are in the same order, Prior is part of GM except P value
    #Output: P - m by 1 matrix containing probability
    #Authors: Zhiwu Zhang
    # Last update: March 10, 2013
    ##############################################################################
    #print("FarmCPU.Prior Started")
    #print("dimension of GM")
    #print(dim(GM))
    
    if(is.null(Prior)& kinship.algorithm!="FARM-CPU")return(P)
    if(is.null(Prior)& is.null(P))return(P)
    
    #get prior position
    if(!is.null(Prior)) index=match(Prior[,1],GM[,1],nomatch = 0)
    
    #if(is.null(P)) P=runif(nrow(GM)) #set random p value if not provided (This is not helpful)
    #print("debug set prior  a")
    
    #Get product with prior if provided
    if(!is.null(Prior) & !is.null(P)  )P[index]=P[index]*Prior[,4]
    
    #print("debug set prior   b")
    return(P)
}#The function FarmCPU.Prior ends here
