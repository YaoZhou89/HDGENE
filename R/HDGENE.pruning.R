`HDGENE.pruning` <- function(GD=NULL,GM=NULL,block = 10000,LD=0.7){
  # Objects: pruning the Genotype based on correlation
  # input: GD: m by n matrix
  #        GM: genotype map information
  #        block: # of markers considering will have LD
  #        LD: correlation, above this threshold will be removed
  CHR = unique(GM[,2])
  SNP.index = NULL
  nchr = 0 
  for(chr in CHR){
    index.chr = which(GM[,2]==chr)
    GD1 = GD[index.chr,]
    SNP.index1 = pruning(GD1,block,LD)
    SNP.index1 = SNP.index1 + nchr
    SNP.index = append(SNP.index,SNP.index1)
    nchr = nchr + length(index.chr)
  }
  return(SNP.index)
}

`pruning` <- function(GD,block,LD){
  m = nrow(GD)
  k = ceiling(m/block)
  GD1.index = NULL
  for( i in 1:k){
    down = (i-1) * block + 1
    up = i * block
    if(up > m) up = m
    GD1 = GD[down:up,]
    GD.index = Blink.LDRemove(GDneo = GD1, LD=LD, Porder = 1:nrow(GD1),block = 10000, LD.num = 10000)
    GD1.index = append(GD1.index,GD.index)
  }
  return(GD1.index)
}


