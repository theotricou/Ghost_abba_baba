#!/usr/bin/env Rscript
# Théo

cat("Starting ms simulation \n")
args = commandArgs(trailingOnly=TRUE)

cat("step 1 : initializing \n")

output=args[1]
# output="."
source(paste(output,"ms_command.R", sep ="/"))
source(paste(output,"sim_parameters", sep ="/"))
library(parallel)

is_abba <- function(seg_site){
  if(sum(seg_site) == 2 & seg_site[1] == seg_site[4]){
    return(1)
  }else{
    return(0)
  }
}

is_baba <- function(seg_site){
  if(sum(seg_site) == 2 & seg_site[1] == seg_site[3]){
    return(1)
  }else{
    return(0)
  }
}

validateandreorderD<-function(arr, dist) {
  # created by Damien
  submat<-dist[arr,arr]
  if (sum(submat==max(submat))==6) {
    diag(submat)<-1
    return(names(sort(apply(submat,1,prod))))
  }
}

getquatuors<-function(tr) {
  # created by Damien
  dist<-cophenetic(tr)
  allquat<-combn(tr$tip.label,4)
  RES<-do.call(rbind,apply(allquat, 2, function(x) validateandreorderD(x, dist)))
  return(RES)
}


leaf_descendants <- function(leaf, comp_leaf) {
  all_anc_leaf <- c(leaf, Ancestors(tree, leaf, type = "all"))
  first_son_leaf <- all_anc_leaf[match(mrca.phylo(tree, c(leaf, comp_leaf)), all_anc_leaf) - 1]
  des_leaf <- c(first_son_leaf,Descendants(tree, as.integer(first_son_leaf), type = "all"))
  return(des_leaf)
}

leaf_ancestors<-function(p, p_ext){
  all_anc <- c(p, Ancestors(tree, p, type = "all"))
  anc_only <- all_anc[1:match(mrca.phylo(tree, c(p, p_ext)), all_anc) - 1]
  return(anc_only)
}

where_is_Raldo<-function(p1,p2,p3,p4){
  if (true_recip %in% leaf_ancestors(p1,p2)) {
    return("P1")
  }else if (true_recip %in% leaf_ancestors(p2,p1)) {
    return("P2")
  }else if (true_recip %in% leaf_ancestors(p3,p1)) {
    return("P3")
  }else if (true_recip %in% leaf_ancestors(p4,p1)) {
    return("P4")
  }else {
    if (!true_recip %in% Descendants(tree,mrca.phylo(tree, c(p1, p4)), type="all")) {
      return("O")
    }else {
      if (true_recip %in% leaf_ancestors(p1,p3)) {
        return("N1")
      }else {
        if (true_recip %in% leaf_ancestors(p1,p4)) {
          return("N2")
        }else{
          return("O")
        }
      }
    }
  }
}

where_is_Daldo<-function(p1,p2,p3,p4){
  if (true_donor %in% leaf_descendants(p1,p2)) {
    return("P1")
   }else if (true_donor %in% leaf_descendants(p2,p1)) {
     return("P2")
   }else if (true_donor %in% leaf_descendants(p3,p1)) {
     return("P3")
   }else if (true_donor %in% leaf_descendants(p4,p1)) {
    return("P4")
  }else {
    if (!true_donor %in% Descendants(tree,mrca.phylo(tree, c(p1, p4)), type="all")) {
      return("O")
    }else {
      if (!true_donor %in% Descendants(tree,mrca.phylo(tree, c(p1, p3)), type="all")) {
        return("N2")
      }else{
        return("N1")
      }
    }
  }
}

N2_diversity<-function(q){
  a<-mrca.phylo(tree, spnd[unlist(q)])
  p1p4Desc<-Descendants(tree, mrca.phylo(tree, spnd[q]),type="tip")[[1]]
  p1p3Desc<-Descendants(tree, mrca.phylo(tree, spnd[q][c(1:3)]),type="tip")[[1]]
  p0<-c(q[4],Ancestors(tree, q[4], type = "all"))
  p4Desc<-Descendants(tree,p0[which(p0 == a)-1],type="tip")[[1]]
  all<-c(p1p3Desc,p4Desc)
  `%notin%` <- Negate(`%in%`)
  N2_tips<-p1p4Desc[p1p4Desc %notin% all]
  n_hidden<-length(N2_tips)
  n_ext<-length(grep("ext>",tree$tip.label[N2_tips]))
  n_uns<-length(grep("uns>",tree$tip.label[N2_tips]))
  return(c("H"=n_hidden, "E"=n_ext, "U"=n_uns, "I"=length(p1p3Desc)))
}

D_stat <- function(stat_simulation, q){
  p1 = which(da$nu == q[1])
  p2 = which(da$nu == q[2])
  p3 = which(da$nu == q[3])
  p4 = which(da$nu == q[4])
  abba <- sum(apply(stat_simulation[c(p1, p2, p3, p4),], 2, function(x) is_abba(x)))
  baba <- sum(apply(stat_simulation[c(p1, p2, p3, p4),], 2, function(x) is_baba(x)))
  if ((abba + baba) != 0) {D = (abba - baba) / (abba + baba)} else {D = 0}
  donor = where_is_Daldo(q[1], q[2], q[3], q[4])
  recip = where_is_Raldo(q[1], q[2], q[3], q[4])
  N2_tips<-N2_diversity(c(q[1], q[2], q[3], q[4]))
  data <- c(
    "P1" = q[1],
    "P2" = q[2],
    "P3" = q[3],
    "P4" = q[4],
    "abba" = abba,
    "baba" = baba,
    "D" = D,
    "Pvalue" = binom.test(c(abba, baba), p = 0.5)$p.value,
    "Donor" = donor,
    "Recip" = recip,
    "dP1P2" = as.integer(strsplit(spnd[mrca.phylo(tree, c(q[1], q[2]))], "@")[[1]][2]),
    "dP1P3" = as.integer(strsplit(spnd[mrca.phylo(tree, c(q[1], q[3]))], "@")[[1]][2]),
    "dP1P4" = as.integer(strsplit(spnd[mrca.phylo(tree, c(q[1], q[4]))], "@")[[1]][2]),
    "N2_M"=N2_tips[[1]],
    "N2_E"=N2_tips[[2]],
    "N2_U"=N2_tips[[3]],
    "N2_I"=N2_tips[[4]]
  )
  return(data)
}

leaf_descendants <- function(leaf, comp_leaf) {
  all_anc_leaf <- c(leaf, Ancestors(tree, leaf, type = "all"))
  first_son_leaf <- all_anc_leaf[match(mrca.phylo(tree, c(leaf, comp_leaf)), all_anc_leaf) - 1]
  des_leaf <- c(first_son_leaf,Descendants(tree, as.integer(first_son_leaf), type = "all"))
  return(des_leaf)
}

leaf_ancestors<-function(p, p_ext){
  all_anc <- c(p, Ancestors(tree, p, type = "all"))
  anc_only <- all_anc[1:match(mrca.phylo(tree, c(p, p_ext)), all_anc) - 1]
  return(anc_only)
}
uniquer<-function(x){
  if (ncol(x$seg_sites[[1]]) == 1) { # ajouter filtre sur site sum(seg) == 1
    return(x$seg_sites[[1]][[1]])
  }
}

is_sister<-function(x){
  donor<-as.character(unlist(x[9]))
  if (donor %in% c("O",'N2','N1')){
    return('Nope')
  }else{
    col<-as.numeric(strsplit(donor,'P')[[1]][2])
    node=as.numeric(unlist(x[col]))
    if (true_donor %in% c(node, Ancestors(tree, node))){
      return('True')
    }else{
      return('Sister')
    }
  }
}

######## function for D3
D3_v2<-function(matrice, t){
  P1<-da[da$nu == t[1],3]
  P2<-da[da$nu == t[2],3]
  P3<-da[da$nu == t[3],3]
  D3<-(matrice[P2,P3]-matrice[P1,P3])/(matrice[P2,P3]+matrice[P1,P3])
  D3_donor<-position_D_in_tree(t[1],t[2],t[3])
  D3_recip<-position_R_in_tree(t[1],t[2],t[3])
  data <- c(
    "P1" = t[1],
    "P2" = t[2],
    "P3" = t[3],
    "dp23" = matrice[P2,P3],
    "dp13" = matrice[P1,P3],
    "D3" = as.numeric(D3),
    "Donor" = D3_donor,
    "Recip" = D3_recip,
    "dP1P2" = as.integer(strsplit(spnd[mrca.phylo(tree, c(t[1], t[2]))], "@")[[1]][2]),
    "dP1P3" = as.integer(strsplit(spnd[mrca.phylo(tree, c(t[1], t[3]))], "@")[[1]][2])
  )
  return(data)
}

validateandreorderD3<-function(arr, dist) {
  # created by Damien
  submat<-dist[arr,arr]
  return(names(sort(apply(submat,1,sum))))
}

gettrios<-function(tr) {
  # created by Damien
  dist<-cophenetic(tr)
  alltrios<-combn(tr$tip.label,3)
  RES<-t(apply(alltrios, 2, function(x) validateandreorderD3(x, dist)))
  return(RES)
}

devil_mambojambo<-function(range){
  return(Reduce('+', lapply(single_trees[as.integer(fromto[range,1]):as.integer(fromto[range,2])], function(x) cophenetic(x)[ord,ord])))
}


position_D_in_tree<-function(p1,p2,p3){
  if (true_donor %in% leaf_descendants(p3,p2)){
    return("P3")
  }else if (true_donor %in% leaf_descendants(p2,p1)){
    return("P2")
  }else if (true_donor %in% leaf_descendants(p1,p2)){
    return("P1")
  }else if (true_donor %in% Descendants(tree,mrca.phylo(tree, c(p1, p3)), type="all")){
    return("N1")
  }else{
    return("O")
  }
}

position_R_in_tree<-function(p1,p2,p3){
  if (true_recip %in% leaf_ancestors(p3,p2)){
    return("P3")
  }else if (true_recip %in% leaf_ancestors(p2,p1)){
    return("P2")
  }else if (true_recip %in% leaf_ancestors(p1,p2)){
    return("P1")
  }else if (true_recip %in% Descendants(tree,mrca.phylo(tree, c(p1, p3)), type="all")){
    return("N1")
  }else{
    return("O")
  }
}


############# function for Dfoil
where_is_Raldoil<-function(p1,p2,p3,p4,p5){
  if (true_recip %in% leaf_ancestors(p1,p2)){
    return("P1")}
  else if (true_recip %in% leaf_ancestors(p2,p1)){
    return("P2")}
  else if (true_recip %in% leaf_ancestors(p3,p4)){
    return("P3")}
  else if (true_recip %in% leaf_ancestors(p4,p3)){
    return("P4")}
  else if (true_recip %in% leaf_ancestors(p5,p1)){
    return("P5")}
  else {
    if (!true_recip %in% Descendants(tree,mrca.phylo(tree, c(p1, p5)), type="all")){
      return("O")
    }else {
      if (true_recip %in% Descendants(tree,mrca.phylo(tree, c(p1, p4)), type="all")){
        if (true_recip %in% leaf_ancestors(p1,p3)){
          return("P1P2")
        }else if (true_recip %in% leaf_ancestors(p3,p1)){
          return("P3P4")
        }else{
          return('O')
        }
      }else{
        if (true_recip %in% leaf_ancestors(p3,p5)){
          return("N2")
        }else{
          return('O')
        }
      }
    }
  }
}

where_is_Daldoil<-function(p1,p2,p3,p4,p5){
  if (true_donor %in% leaf_descendants(p1,p2)){
    return("P1")
  }else if (true_donor %in% leaf_descendants(p2,p1)){
    return("P2")
  }else if (true_donor %in% leaf_descendants(p3,p4)){
    return("P3")
  }else if (true_donor %in% leaf_descendants(p4,p3)){
    return("P4")
  }else if (true_donor %in% leaf_descendants(p5,p3)){
    return("P5")
  }else{
    if (!true_donor %in% Descendants(tree,mrca.phylo(tree, c(p1, p5)), type="all")){
      return("O")
    }else {
      if (!true_donor %in% Descendants(tree,mrca.phylo(tree, c(p1, p4)), type="all")){
        return("N2")
      }else{
        if (true_donor %in% leaf_descendants(p1,p3)){
          return("P1P2")
        }else if (true_donor %in% leaf_descendants(p3,p1)){
          return("P3P4")
        }
      }
    }
  }
}

eval_Dfoil<-function(a,pa,b,pb,c,pc,d,pd){
  one="0"
  two="0"
  tre="0"
  fou="0"
  pvalue=0.005
  if (pa<=pvalue){
    if (a>0){
      one='+'
    }else{
      one='-'
    }
  }
  if (pb<=pvalue){
    if (b>0){
      two='+'
    }else{
      two='-'
    }
  }
  if (pc<=pvalue){
    if (c>0){
      tre='+'
    }else{
      tre='-'
    }
  }
  if (pd<=pvalue){
    if (d>0){
      fou='+'
    }else{
      fou='-'
    }
  }
  res=c('one'=one,
  "two"=two,
  'tre'=tre,
  'fou'=fou)
  return(res)
}

EVfoil<-c("+++0","+0++","--0+","-0++","++-0","0+--","--0-","0---","++00","++00","--00","--00",'0000','0000','0000','0000')
patoil<-c("P1-P3","P3-P1","P1-P4","P4-P1","P2-P3","P3-P2","P2-P4","P4-P2","P1P2-P3","P3-P1P2","P4-P1P2","P1P2-P4","P1-P2","P2-P1","P3-P4","P4-P3")
dfoil_dict<-as.data.frame(cbind(EVfoil,patoil))

is_true_Dfoil<-function(pattern, EV){
  if (pattern %in% dfoil_dict[which(dfoil_dict == EV),1]){
    return("TP")
  }else if (pattern %in% dfoil_dict[,1]){
    return("FP")
  }else if (!pattern %in% dfoil_dict[,1]){
    return('New_pattern')
  }else if (EV %in% dfoil_dict[,1]){
    return("FN")
  }else{
    return("none")
  }
}

filter<-function(x){
  if (sum(x)!=5){
    if(sum(x)!=0){
      if(x[5]==1){
        return(+(!x))
      }else{
        return(x)
      }
    }
  }
}

sitesquintet<-function(site){
  if (sum(site) == 1){
    if (site[1]==1) {
      return("baaaa")
    }else if(site[2]==1) {
      return("abaaa")
    }else if(site[3]==1) {
      return("aabaa")
    }else if(site[4]==1) {
      return("aaaba")
    }
  }else if (sum(site) == 2){
    if (site[1]+site[3] == 2){
      return("babaa")
    }else if (site[1]+site[4] == 2){
      return("baaba")
    }else if (site[2]+site[3] == 2){
      return("abbaa")
    }else if (site[2]+site[4] == 2){
      return("ababa")
    }
  }else if (sum(site) == 3){
    if (site[1]==0) {
      return("abbba")
    }else if(site[2]==0) {
      return("babba")
    }else if(site[3]==0) {
      return("bbaba")
    }else if(site[4]==0) {
      return("bbbaa")
    }
  }else if (sum(site) == 4){
    if (site[1]==0) {
      return("baaaa")
    }else if(site[2]==0) {
      return("abaaa")
    }else if(site[3]==0) {
      return("aabaa")
    }else if(site[4]==0) {
      return("baaba")
    }
  }
}

Dfoiler<-function(site, top){
  # top=nn
  nn=unlist(lapply(as.character(top),function(x) da[which(da$na == x),2]))
  top=unlist(lapply(as.character(top),function(x) which(da$na == x)))
  # top=match(top, da$nu)
  names<-c("aaaba","aabaa","abaaa","baaaa","babaa","baaba","abbaa","ababa","bbbaa","bbaba","babba","abbba")
  pat<-as.data.frame(names)
  colnames(pat)<-"Var1"
  if (nrow(site) == 5){
    data<-t(site[top,])
  }else{
    data<-as.data.frame(do.call('rbind',apply(site[top,],2,filter)))
  }

  p<-as.data.frame(table(unlist(apply(data,1, function(x) sitesquintet(x)))))
  d<-t(merge(p,pat,all=TRUE))
  d[is.na(d)]<-0
  colnames(d)<-as.matrix(d[1, ])
  d<-as.data.frame(t(d[-1, ]))
  d<-lapply(d, function(x) as.numeric(as.character(x)))

  Dfo<- ((d$babaa + d$bbbaa + d$ababa + d$aaaba) - (d$baaba + d$bbaba + d$abbaa + d$aabaa)) /
        ((d$babaa + d$bbbaa + d$ababa + d$aaaba) + (d$baaba + d$bbaba + d$abbaa + d$aabaa))
  PDfo<-binom.test(c((d$babaa + d$bbbaa + d$ababa + d$aaaba), (d$baaba + d$bbaba + d$abbaa + d$aabaa)), p = 0.5)$p.value

  Dil<- ((d$abbaa + d$bbbaa + d$baaba + d$aaaba) - (d$ababa + d$bbaba + d$babaa + d$aabaa)) /
        ((d$abbaa + d$bbbaa + d$baaba + d$aaaba) + (d$ababa + d$bbaba + d$babaa + d$aabaa))
  PDil<-binom.test(c((d$abbaa + d$bbbaa + d$baaba + d$aaaba), (d$ababa + d$bbaba + d$babaa + d$aabaa)), p = 0.5)$p.value

  Dfi<- ((d$babaa + d$babba + d$ababa + d$abaaa) - (d$abbaa + d$abbba + d$baaba + d$baaaa)) /
        ((d$babaa + d$babba + d$ababa + d$abaaa) + (d$abbaa + d$abbba + d$baaba + d$baaaa))
  PDfi<-binom.test(c((d$babaa + d$babba + d$ababa + d$abaaa), (d$abbaa + d$abbba + d$baaba + d$baaaa)), p = 0.5)$p.value

  Dol<- ((d$baaba + d$babba + d$abbaa + d$abaaa) - (d$ababa + d$abbba + d$babaa + d$baaaa)) /
        ((d$baaba + d$babba + d$abbaa + d$abaaa) + (d$ababa + d$abbba + d$babaa + d$baaaa))
  PDol<-binom.test(c((d$baaba + d$babba + d$abbaa + d$abaaa), (d$ababa + d$abbba + d$babaa + d$baaaa)), p = 0.5)$p.value

  Don<-where_is_Daldoil(nn[1],nn[2],nn[3],nn[4],nn[5])
  Rec<-where_is_Raldoil(nn[1],nn[2],nn[3],nn[4],nn[5])

  tEv=eval_Dfoil(Dfo,PDfo,Dil,PDil,Dfi,PDfi,Dol,PDol)

  res<-c("P1"=nn[1],
        "P2"=nn[2],
        "P3"=nn[3],
        "P4"=nn[4],
        "P5"=nn[5],
        "Dfo"=Dfo, "Pvalue_Dfo"=PDfo,
        "Dil"=Dil, "Pvalue_Dil"=PDil,
        "Dfi"=Dfi, "Pvalue_Dfi"=PDfi,
        "Dol"=Dol, "Pvalue_Dol"=PDol,
        "Donor"=Don,
        "Recip"=Rec,
        "dP1P2"=mdist[top[1],top[2]],
        "dP3P4"=mdist[top[3],top[4]],
        "dP1P4"=mdist[top[1],top[4]],
        "dP1P5"=mdist[top[1],top[5]],
        "tEv"=tEv,
        "EV"=paste(Don, Rec,sep="-")
      )
  return(res)
}

getquintet<-function(x){
  submat<-mdist[x,x]
  if (identical(as.numeric(table(submat)),c(5,2,2,8,8))){
    return(names(sort(apply(submat, 1, sum))))
  }
}

is_true_Dfoil<-function(pattern, EV){
  if (pattern %in% dfoil_dict[which(dfoil_dict$patoil == EV),1]){
    return("TP")
  }else if (!EV %in% dfoil_dict[,2] & pattern=="0000"){
    return("TN")
  }else if (pattern %in% dfoil_dict[,1]){
    return("FP")
  }else if (EV %in% dfoil_dict[,2]){
    return("FN")
  }else if (!pattern %in% dfoil_dict[,1]){
    return('New_pattern')
  }else{
    return("none")
  }
}

################################################################################

### D stat

tree <- read.tree(file.path(output, "spe_tree"))
extant<-grep("ali>",tree$tip.label)
tr<-keep.tip(tree, tree$tip.label[extant])
da<-data.frame(na=tree$tip.label[extant], nu=extant, ge=paste("s",1:length(extant),sep=""))
# max_d_to_extant_outgroup <- max(cophenetic(tr))/2

cat("step 2 : running simulation \n")
# N_SIMULATION=50000
# rep = simulate(model, nsim = N_SIMULATION, cores = )
if (SEED == 0) {rep = simulate(model, nsim = N_SIMULATION, cores = N_CORE)
} else {
  library(parallel)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(SEED)
  M <- N_CORE
  s <- .Random.seed
  for (i in 1:M){s <- nextRNGStream(s)}
  rep <- mclapply(X = 1:N_SIMULATION,
               FUN = function(x) simulate(model),
               mc.cores = M,
               mc.set.seed = TRUE
  )
}

uniq<-mclapply(rep, function(x) uniquer(x),mc.cores = N_CORE)

single_trees<-unlist(lapply(which(uniq != "NULL"), function(x) rep[[x]]$trees[[1]]))

sites <- do.call("cbind", uniq)
# sites2 <- sites[, colSums(sites != 0) > 1]
pretopologiesD <- getquatuors(tr)
topologiesD <- apply(pretopologiesD,c(1,2), function(x) da[as.character(da$na) == as.character(x),2])
spnd<-c(tree$tip.label, tree$node.label)
true_donor = which(spnd == name_donor) # from variable in sourced file
true_recip = which(spnd == name_recip) # from variable in sourced file

cat("step 3 : computing D stat \n")

if (N_CORE > nrow(topologiesD)){
  N_c<-nrow(topologiesD)
  results<-as.data.frame(do.call(rbind,mclapply(as.data.frame(t(topologiesD)), function(x) D_stat(sites, x),mc.cores = N_c)))
}else{
  results<-as.data.frame(do.call(rbind,mclapply(as.data.frame(t(topologiesD)), function(x) D_stat(sites, x),mc.cores = N_CORE)))
}
results$is_sister<-apply(results,1, function(x) is_sister(x))

cat("step 4 : output D\n")

outfile_d <- paste(output, "data.txt", sep = "/")
write.table(results, outfile_d, sep = "\t", row.names = F, append = F, quote=F)

# cat("step 5 : computing Dfoid \n")
# mdist<-cophenetic(tr)
# comb_D<-combn(tr$tip.label,5)
# if (ncol(comb_D) > 1){
#   pretopologiesDfoil<-do.call('rbind',apply(t(comb_D),1, function(x) getquintet(x)))
#   topologiesDfoil<-apply(pretopologiesDfoil,c(1,2), function(x) da[as.character(da$na) == as.character(x),2])
#   results2<-as.data.frame(do.call("rbind",mclapply(as.data.frame(t(pretopologiesDfoil)),function(x) Dfoiler(sites,x),mc.cores = N_CORE)))
#   # results2$TF<-apply(results2,1,function(x) is_true_Dfoil(do.call('paste',c(as.list(x[c(20:23)]),sep="")),as.character(x[24])))
# }else if (ncol(comb_D) == 1){
#   pretopologiesDfoil<-getquintet(t(comb_D))
#   topologiesDfoil<-unlist(lapply(pretopologiesDfoil, function(x) da[as.character(da$na) == as.character(x),2]))
#   results2<- as.data.frame(t(Dfoiler(sites,c(pretopologiesDfoil))))
# }
#
# # cat("step 6 : output Dfoil\n")
#
# outfile <- paste(output, "data_Dfoil.txt", sep = "/")
# write.table(results2, outfile, sep = "\t", row.names = F, append = F, quote=F)
#
cat("End \n")

#GNU Terry Pratchett
