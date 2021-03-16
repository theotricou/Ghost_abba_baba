## Dec 11, 2020
## -> Generate species tree
## -> Create all quartets
## -> For each quartet 
##     - compute all possible transfers that would make this quartet have signification Dstat 
##     - get the relative distance between outgroup and ingroup
## -> Separate possible transfers in two lists: ghost transfers (coming from MIDgroup)
##    and nonghost transfers (coming from INgroup))
## -> sample an introgression (change effect of distance on proba of introg)
## -> compute the proportion of quartets that would have significant Dstat
## -> compute the proportion of erroneous interpretation

################
# DEPENDENCIES #
################
require(ape)
require(phangorn)


#############
# FUNCTIONS #
#############
validateandreorder<-function(arr, dist) {
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
  RES<-do.call(rbind,apply(allquat, 2, function(x) validateandreorder(x, dist)))
  return(RES)
}
getstartstopofedges<-function(mat,xx) {
    if (!is.null(mat)) {
        return(t(apply(mat,1,function(x,xx) c(x,xx[x]),xx=xx)))
    }
    else return(NULL)
}
#TEST IF TWO INTERVALS OVERLAP
overlap<-function(ab,cd) {
    #true if the two intervals overlap
    return(!(ab[2]<cd[1] | cd[2]<ab[1]))
}
events.donor.receiver<-function(quat, tr, xy) { #xy are x and y coordinates in the plot
    p1<-quat[1]
    p2<-quat[2]
    p3<-quat[3]
    p4<-quat[4]
    #################
    ### P1 AND P2 ###
    #################
    p1p2<-mrca.phylo(tr, c(p1,p2))
    p1nb<-which(tr$tip.label==p1)
    p2nb<-which(tr$tip.label==p2)
    desc<-sapply(Children(tr, p1p2), function(x,tree) unique(c(x,Descendants(tree, x, "all"))), tree=tr, simplify=FALSE)
    if (is.element(p2nb, desc[[1]])) desc<-rev(desc) ##so that p1 remains in first position. 
    ##WE STORE ALL POSSIBLE BRANCHES
    P1DONOR<-tr$edge[match(desc[[1]],tr$edge[,2]),, drop=FALSE]
    P2DONOR<-tr$edge[match(desc[[2]],tr$edge[,2]),,drop=FALSE]

    ## NOW WE LOOK AT THE SUM OF BL AS RECEIVERS. fOR AN INTROGRESSION TO BE DETECTED, 
    ## THE RECEIVING BRANCH MUST BE AN ANCESTOR OF P1 (RESP. P2) BUT NO OLDER THAN THE MRCA OF P1P2
    p1.up<-c(p1nb, Ancestors(tr, p1nb))
    p1.up.small<-p1.up[1:(which(p1.up==p1p2)-1)]
    p2.up<-c(p2nb, Ancestors(tr, which(tr$tip.label==p2)))
    p2.up.small<-p2.up[1:(which(p2.up==p1p2)-1)]
    ##WE STORE ALL POSSIBLE BRANCHES
    P1RECEIVER<-tr$edge[match(p1.up.small,tr$edge[,2]),, drop=FALSE]
    P2RECEIVER<-tr$edge[match(p2.up.small,tr$edge[,2]),, drop=FALSE]

    #################
    ###     P3    ###
    #################
    ## FOR P3 and P4 IT IS A BIT MORE TRICKY BECAUSE WE MUST EXCLUDE BRANCHES THAT END BEFORE MRCA(P1,P2)
    ## AND INTEGRATE PART OF THOSE THAT ARE CROSSING THIS LINE
    ## NOTE: P4 WILL NEVER BE THE RECEIVER BECAUSE IT IS GHOST.
    p1p3<-mrca.phylo(tr, c(p1,p3)) #same than p2p3
    p3nb<-which(tr$tip.label==p3)
    desc1.3<-sapply(Children(tr, p1p3), function(x,tree) unique(c(x,Descendants(tree, x, "all"))), tree=tr, simplify=FALSE)
    if (is.element(p1nb, desc1.3[[1]])) desc1.3<-rev(desc1.3)
    desc1.3<-desc1.3[[1]]
    subedge<-tr$edge[match(desc1.3,tr$edge[,2]),, drop=FALSE]
    ##
    xp1p2<-xy$xx[p1p2] #time of p1p2 speciation. p3 donor or receiver must be more recent.
    xp1<-xy$xx[p1nb]
    ##
    subedgetime<-cbind(subedge, t(apply(subedge,1,function(x,xx) xx[x], xx=xy$xx)))
    subedgetime<-subedgetime[(subedgetime[,3]>xp1p2)|(subedgetime[,4]>xp1p2),, drop=FALSE] ##select those that have one side at least in the range encounterd by p1p2
    subedgetime[subedgetime[,3]<xp1p2,3]<-xp1p2 #change start time of the edges traversing the line to only count possible bl.
    ##WE STORE ALL POSSIBLE BRANCHES    
    P3DONOR<-subedge

    p3.up<-c(p3nb, Ancestors(tr, p3nb))
    p3.up.subedgetime<-match(p3.up, subedgetime[,2])
    subedgetime.up<-subedgetime[p3.up.subedgetime[!is.na(p3.up.subedgetime)],, drop=FALSE]
    ##WE STORE ALL POSSIBLE BRANCHES    
    P3RECEIVER<-subedgetime.up[,1:2, drop=FALSE]
    #################
    ###     P4    ###
    #################
    ## P4 represents a possible ghost lineage that branches between the ancestor of (P1P2P3) and the ancestor of (P1P2P3P4)
    p1p4<-mrca.phylo(tr, c(p1,p4))
    p4nb<-which(tr$tip.label==p4)
    possible.p4.nodes<-intersect(Descendants(tr, p1p4, "all"), Ancestors(tr, p1p3))
    if (length(possible.p4.nodes)==0) {
        P4DONOR<-NULL
    }
    else {
        p4.to.check<-setdiff(unique(unlist(Descendants(tr, possible.p4.nodes,"all"))), Descendants(tr, p1p3,"all"))
        subedge.p4<-tr$edge[match(p4.to.check,tr$edge[,2]),]
        subedgetime.p4<-cbind(subedge.p4, t(apply(subedge.p4,1,function(x,xx) xx[x], xx=xy$xx)))
        subedgetime.p4<-subedgetime.p4[(subedgetime.p4[,3]>xp1p2)|(subedgetime.p4[,4]>xp1p2),,drop=FALSE] ##select those that have one side at least in the range encounterd by p1p2
        subedgetime.p4[subedgetime.p4[,3]<xp1p2,3]<-xp1p2 #change start time of the edges traversing the line to only count possible bl.
        P4DONOR<-subedgetime.p4[,1:2,drop=FALSE]
        if (nrow(P4DONOR)==0) P4DONOR<-NULL
    }

    ###############################################################################
    ### STORE DISTANCE BETWEEN OUTGROUP AND P1P2P3 (in proportion of full size) ###
    ###############################################################################
    relativedist2outgroup<-(xy$xx[p1p3]-xy$xx[p1p4])/(xy$xx[p1nb]-xy$xx[p1p4])
    relativedistp1top3<-(xy$xx[p1nb]-xy$xx[p1p3])/(xy$xx[p1nb]-xy$xx[p1p4])
    relativedistp1top2<-(xy$xx[p1nb]-xy$xx[p1p2])/(xy$xx[p1nb]-xy$xx[p1p4])
    ############################################
    ### COMPUTE EDGE-EDGE POSSIBLE TRANSFERS ###
    ############################################
    DONOR.RECEIVER<-list(P1DONOR=P1DONOR, P2DONOR=P2DONOR,P3DONOR=P3DONOR,P4DONOR=P4DONOR,P1RECEIVER=P1RECEIVER,P2RECEIVER=P2RECEIVER,P3RECEIVER=P3RECEIVER)
    DONOR.RECEIVER<-lapply(DONOR.RECEIVER, getstartstopofedges, xx=xy$xx)
    ##THEN FOR EACH POSSIBLE TRABSFERSTORY (P1P3, P3P1, etc.) we only keep those possible in terms of temporal co-occurence of edges.
    getpossibletransfers<-function(mat1,mat2) {
        if (!is.null(mat1)&!is.null(mat2)) {
            TFmat<-apply(mat1[,3:4, drop=FALSE],1,function(x,m2) apply(m2[,3:4, drop=FALSE], 1, overlap,cd=x), m2=mat2)
            if (is.null(dim(TFmat))) TFmat<-t(t(TFmat))
            else TFmat<-t(TFmat)
            ##in the matrix, rows are rows of mat1 and cols are cols of 
            matchingedges<-which(TFmat, arr.ind=T)
            ##and for each match, we recover its condensed name.
            res<-apply(matchingedges,1,function(x,m1,m2) paste(paste(m1[x[1],1:2],collapse="-"),paste(m2[x[2],1:2],collapse="-"),sep="->"), m1=mat1, m2=mat2)
            return(res)
        }
    }
    P1P3<-getpossibletransfers(DONOR.RECEIVER$P1DONOR, DONOR.RECEIVER$P3RECEIVER)
    P2P3<-getpossibletransfers(DONOR.RECEIVER$P2DONOR, DONOR.RECEIVER$P3RECEIVER)
    P3P1<-getpossibletransfers(DONOR.RECEIVER$P3DONOR, DONOR.RECEIVER$P1RECEIVER)
    P3P2<-getpossibletransfers(DONOR.RECEIVER$P3DONOR, DONOR.RECEIVER$P2RECEIVER)
    P4P1<-getpossibletransfers(DONOR.RECEIVER$P4DONOR, DONOR.RECEIVER$P1RECEIVER)
    P4P2<-getpossibletransfers(DONOR.RECEIVER$P4DONOR, DONOR.RECEIVER$P2RECEIVER)
    TRANSFERSINVOLVED<-list(NOGHOST=c(P1P3,P2P3,P3P1,P3P2),GHOST=c(P4P1,P4P2))

    return(list(relativedist2outgroup=relativedist2outgroup, relativedistp1top3=relativedistp1top3,relativedistp1top2=relativedistp1top2, transfersinvolved=TRANSFERSINVOLVED))
}

SelectIntrogressionDistanceEffect<-function(tree, xy, alpha, lwd, ForbiddenNodes=NULL) {
    #for a given tree, sample an introgression, with proba proportional to nb of branches at this time.
    #sampling: we put all branches each after each other and sample in this new array uniformly. Should work. 
    #alpha controls the effect of phylog distance when choosing the recipient. high alpha mean that the effect is big. If alpha=0, there is no effect of distance 
    ##  We can forbid some nodes so that branches connected to these nodes
    ### will never be allowed
    testforbid<-0
    while(testforbid==0) {
        sample<-runif(1,0,max(cumsum(tree$edge.length)))
        if (sample<cumsum(tree$edge.length)[1]) {
            quellebranche<-1
            oudanslabranche<-sample
        }
        else {
            quellebranche<-max(which(sample>cumsum(tree$edge.length)))+1
            oudanslabranche<-sample-cumsum(tree$edge.length)[quellebranche-1] #disatnce au début de la branche.
        }
        ##we decide that this is the donor branch. The receiver branch will be another one present at exactly the same moment. 
        edgeandtime<-getstartstopofedges(tree$edge, xy$xx)
        donorbranch<-edgeandtime[quellebranche,]
        edgeandtime2<-edgeandtime[-quellebranche,] #matrix without donot to not pick it.
        momentintrogression<-donorbranch[3]+oudanslabranche
        possiblerecipient<-apply(edgeandtime2[,3:4],1, function(x,m) (m>=x[1])&(m<x[2]), m=momentintrogression)

        #Distance To Possible Recipients = 
        MRCAofpossiblerecipients<-sapply(edgeandtime2[possiblerecipient,2], function(x, tr,tip) getMRCA(tr, c(x,tip)),tr=tree, tip=donorbranch[2])
        ageofMRCA<-edgeandtime2[match(MRCAofpossiblerecipients, edgeandtime2[,1]),3]
        DtoRecip<-(xy$xx[donorbranch[1]]+oudanslabranche)-ageofMRCA
        if (length(DtoRecip)>1) SampledOne<-sample(1:length(DtoRecip),1,prob=exp(-alpha * DtoRecip/sum(DtoRecip)))
        else SampledOne<-1
        recipientbranch<-edgeandtime2[possiblerecipient,,drop=FALSE][SampledOne,]
        if (sum(is.element(c(donorbranch,recipientbranch),ForbiddenNodes))==0) testforbid<-1
        else {
            cat("TRY AGAIN\n")
        }
    }

    #uncomment after if you want to see the transfers on the tree.
	arrows(xy$xx[donorbranch[1]]+oudanslabranche, xy$yy[donorbranch[2]],xy$xx[donorbranch[1]]+oudanslabranche, xy$yy[recipientbranch[2]], angle=12, length=0.1, col="red", lwd=lwd)
    return(paste(paste(donorbranch[1:2],collapse="-"),paste(recipientbranch[1:2],collapse="-"),sep="->"))
}


#for a given transfer event this function returns for each 
GetDistribDistanceToOutgroup<-function(event, allscenario) {
	#Number of quartets that will show significant Dstats because of ingroup introgressions: 
	noghost<-unlist(lapply(allscenario[unlist(lapply(allscenario, function(x, event) sum(is.element(x$transfersinvolved$NOGHOST, event)), event=event))==1], function(x) x$relativedist2outgroup))
	#Number of quartets that will show significant Dstats because of midgroup introgressions: 
	ghost<-unlist(lapply(allscenario[unlist(lapply(allscenario, function(x, event) sum(is.element(x$transfersinvolved$GHOST, event)), event=event))==1], function(x) x$relativedist2outgroup))
	return(list(noghost=noghost, ghost=ghost))
}

SummarizeResults<-function(res) {
	#proportion d'introgressions qui vont être détectables (avec un ou plusieurs quartets)
	ProportionDetectable<-sum(apply(res,1,sum)!=0)/nrow(res)
	PropErroneousInterpretALL<-apply(res,1,function(x) x[2]/sum(x))
	PropErroneousInterpretMEAN<-mean(PropErroneousInterpretALL, na.rm=TRUE)
	PropErroneousInterpretSD<-sd(PropErroneousInterpretALL, na.rm=TRUE)
	PropErroneousInterpretMEDIAN<-median(PropErroneousInterpretALL, na.rm=TRUE)
	return(unlist(list(PropIntrogDetectable=ProportionDetectable, meanPropErroneousInterp=PropErroneousInterpretMEAN, sdPropErroneousInterp=PropErroneousInterpretSD, medianPropErroneousInterp=PropErroneousInterpretMEDIAN)))
}
colless<-function(tr) {
    leftright4all<-t(sapply((Ntip(tr)+1):(Ntip(tr)+Nnode(tr)),function(x,tree) Children(tree,x), tree=tr))
    leftright4all.size<-apply(leftright4all,2,function(x,tree) unlist(lapply(Descendants(tree, x, "tips"), length)),tree=tr)
    index<-sum(abs(apply(leftright4all.size,1,diff)))
    return(index)
}
ComputeAll<-function(N=20, pextinction=runif(1), Nreal=sample(seq(20,100,by=20),1), alpha=sample(c(0,1,10,100,1000),1), Nintrog=100, fileout="res_output") {
	VALINFOS<-c(pextinction, alpha, Nintrog)
	tr<-rphylo(Nreal,1,pextinction, fossils=TRUE)
	extant<-paste("t",1:N, sep="") #easy because of the way rphylo is naming tips.
	extant.tr<-keep.tip(tr, extant)
	coll<-colless(extant.tr) ##how much is the tree balanced
    ##
    Ntotalold<-Ntip(tr)
    tr<-extract.clade(tr, node=getMRCA(tr, extant))
    Ntotal<-Ntip(tr)
	TREEINFOS<-c(N, Nreal, Ntotal,coll)
    cat(paste("Ntotalold",Ntotalold, "\n",sep=" - "))
    cat(paste("Ntotal",Ntotal, "\n", sep=" - "))

	#get quartets
	quatuors <- getquatuors(extant.tr)
	QUATINFOS <-nrow(quatuors)
	#get coordinates in plot (easier for bl computations)
	plot(tr, show.tip.label=FALSE)
	xy<-get("last_plot.phylo", envir = .PlotPhyloEnv)
	#get all transfers that would make each quartet having a significant Dstat
	ScenrarioPerQuartet<-apply(quatuors, 1, events.donor.receiver, tr=tr, xy=xy)
	#Simulate Nintrog introgression events
	events<-replicate(Nintrog, SelectIntrogressionDistanceEffect(tr, xy, alpha=alpha, lwd=1))
	result<-sapply(events, GetDistribDistanceToOutgroup, allscenario=ScenrarioPerQuartet, simplify=FALSE)
	all.threshold<-sapply(seq(0,0.9,by=0.1), function(tt) do.call(rbind, lapply(result, function(x) unlist(lapply(x, function(y,th) sum(y>=th), th=tt)))), simplify=F)
	ALLRESULTS<-array(do.call(rbind, lapply(all.threshold, SummarizeResults)))
	##WE CAN WRITE EVERYTHING DOWN IN ONE LINE: 
	#TAMBOUILLE POUR GARDER TOUS LES SCORES : 
	NOGHOST<-lapply(result, function(x) x$noghost)
	GHOST<-lapply(result, function(x) x$ghost)
	NOGHOST.save<-paste(paste(rep(names(NOGHOST),times=unlist(lapply(NOGHOST, length))),unlist(NOGHOST), sep="|"),collapse=";")
	GHOST.save<-paste(paste(rep(names(GHOST),times=unlist(lapply(GHOST, length))),unlist(GHOST), sep="|"),collapse=";")

	Res<-paste(c(as.character(c(VALINFOS, TREEINFOS, QUATINFOS, ALLRESULTS)), GHOST.save, NOGHOST.save), collapse=" ")
	cat(Res, sep="\n", file=fileout, append=TRUE)
	return(c(VALINFOS, TREEINFOS))
}
