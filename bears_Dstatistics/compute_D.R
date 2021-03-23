#!/usr/bin/env Rscript
# Theo Tricou

# compute D-statistique for Bears quartet

library('seqinr')
files <- list.files(pattern = ".snp_sites.aln")
filtering<-function(x){
  x2=x
  # x2=x[!x %in% c("n", '-')]
  if (length(unique(x2)) != 1){
    if (length(x2)>3){
      if (length(which(table(x2) >= 2)) > 1){
        return(TRUE)
      }else{FALSE}
    }else{FALSE}
  }else{FALSE}
}

Der<-function(x){
  if (x[1] == x[3]){
    return("baba")
  }else if(x[2] == x[3]){
    return("abba")
  }else if(x[2] == x[1]){
    return("bbaa")
  }
}

DD<-function(x){
  d<-as.data.frame(as.matrix(read.alignment(x, format = "fasta")))
  d2=d[,which(unlist(lapply(d,function(x) filtering(x))))]
  temp=table(unlist(lapply(d2, function(x) Der(x))))
  D=((temp['abba']-temp['baba'])/(temp['abba']+temp['baba']))
  return(c("abba"=as.integer(temp['abba']), "baba"=as.integer(temp['baba']), "D"=as.numeric(D)))
}

# compute block D-stat and Nabba and Nbaba
library(parallel)
# Ds=lapply(files, function(x) DD(x))
Ds= mclapply(files, function(x) DD(x), mc.cores = 6) # parallel run

res=as.data.frame(do.call('rbind',Ds))

abba=sum(res$abba, na.rm=T)
baba=sum(res$baba, na.rm=T)

# D-stat
D=(abba-baba)/(abba+baba)
D

# Z-score
n_blocks <- nrow(res)
D_sd <- sd(res$D, na.rm=T)
D_err <- D_sd/sqrt(n_blocks)
D_Z <- D / D_err
D_Z

prin(D, D_Z)
