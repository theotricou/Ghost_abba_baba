# simulate alpha one phylogenetic tree
require(ape)
require(phangorn)


source(dstat-topology.R)

getstartstopofedges<-function(mat,xx) {
    if (!is.null(mat)) {
        return(t(apply(mat,1,function(x,xx) c(x,xx[x]),xx=xx)))
    }
    else return(NULL)
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
	# arrows(xy$xx[donorbranch[1]]+oudanslabranche, xy$yy[donorbranch[2]],xy$xx[donorbranch[1]]+oudanslabranche, xy$yy[recipientbranch[2]], angle=12, length=0.1, col="red", lwd=lwd)
    return(paste(paste(donorbranch[1:2],collapse="-"),paste(recipientbranch[1:2],collapse="-"),sep="->"))
}

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


alpha_sim<-function(alpha, n_tr){
	#Simulate Nintrog introgression events
	events<-replicate(n_tr, SelectIntrogressionDistanceEffect(tr, xy, alpha=alpha, lwd=1))
	dist_events<-function(x){
		fr_to<-strsplit(x, "-")[[1]]
		dist= length(nodepath(tr, as.numeric(fr_to[2]), as.numeric(fr_to[4]))) -2
		if (dist > 1){
			return(dist)
		}
	}
	return(mean(unlist(lapply(events,function(x) dist_events(x))), na.rm=T))
}

force.ultrametric<-function(tree,method=c("nnls","extend")){
    method<-method[1]
    if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
        rooted=TRUE,trace=0)
    else if(method=="extend"){
        h<-diag(vcv(tree))
        d<-max(h)-h
        ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
            y=tree$edge[,2])
        tree$edge.length[ii]<-tree$edge.length[ii]+d
    } else
        cat("method not recognized: returning input tree\n\n")
    tree
}


alphas=c(0,5,10,50,100,500,1000)



tree<-read.tree('bear_Hailer_2012')
tr=force.ultrametric(tree)
quatuors <- getquatuors(tr)
plot(tr, show.tip.label=T, plot=F)
xy<-get("last_plot.phylo", envir = .PlotPhyloEnv)
ScenrarioPerQuartet<-apply(quatuors, 1, events.donor.receiver, tr=tr, xy=xy)
data = lapply(alphas, function(x) unlist(replicate(100,alpha_sim(x, 1))))
transfer_score=3

png(filename = "bear_score_no_sister_num.png", width = 600, height = 600, units = "px", pointsize = 12, bg = "white")
boxplot(data, names=alphas,
	main="Mean divergence between donors and recipients\nin bear phylogeny (Hailer et al. 2012)",
	xlab="Alpha", ylab="Simulated mean divergence between donors and recipients")

legend("topright", legend = c(paste("observed = ", transfer_score, sep="")),
	col = c("red", "blue"),
	pch = c(19,19), pt.cex = 2, text.col = "black")
abline(h=transfer_score, col="red")
dev.off()




tree<-read.tree('spider_Leduc-Robert_2018')
tr=force.ultrametric(tree)
quatuors <- getquatuors(tr)
plot(tr, show.tip.label=T, plot=F)
xy<-get("last_plot.phylo", envir = .PlotPhyloEnv)
ScenrarioPerQuartet<-apply(quatuors, 1, events.donor.receiver, tr=tr, xy=xy)
data = lapply(alphas, function(x) unlist(replicate(100,alpha_sim(x, 4))))
transfer_score=mean(c(2,3,4,4))

png(filename = "spider_score_no_sister_num.png", width = 600, height = 600, units = "px", pointsize = 12, bg = "white")
boxplot(data, names=alphas,
	main="Mean divergence between donors and recipients\nin spider phylogeny (Leduc-Robert et al. 2018)",
	xlab="Alpha", ylab="Simulated mean divergence between donors and recipients")

legend("topright", legend = c(paste("observed = ", transfer_score, sep="")),
	col = c("red", "blue"),
	pch = c(19,19), pt.cex = 2, text.col = "black")
abline(h=transfer_score, col="red")
dev.off()



tr<-read.tree('bos_Wu_2018')
quatuors <- getquatuors(tr)
plot(tr, show.tip.label=T, plot=F)
xy<-get("last_plot.phylo", envir = .PlotPhyloEnv)
ScenrarioPerQuartet<-apply(quatuors, 1, events.donor.receiver, tr=tr, xy=xy)
data = lapply(alphas, function(x) unlist(replicate(100,alpha_sim(x, 5))))
transfer_score=mean(c(5,5,4,4,4))

png(filename = "bos_score_no_sister_num.png", width = 600, height = 600, units = "px", pointsize = 12, bg = "white")
boxplot(data, names=alphas,
	main="Mean divergence between donors and recipients\nin bos complexe phylogeny (Wu et al. 2018)",
	xlab="Alpha", ylab="Simulated mean divergence between donors and recipients ")

legend("topright", legend = c(paste("observed = ", transfer_score, sep="")),
	col = c("red", "blue"),
	pch = c(19,19), pt.cex = 2, text.col = "black")
abline(h=transfer_score, col="red")
dev.off()


tr<-read.tree('mosquitoes_Fontaine_2015')
quatuors <- getquatuors(tr)
plot(tr, show.tip.label=T, plot=F)
xy<-get("last_plot.phylo", envir = .PlotPhyloEnv)
ScenrarioPerQuartet<-apply(quatuors, 1, events.donor.receiver, tr=tr, xy=xy)
data = lapply(alphas, function(x) unlist(replicate(100,alpha_sim(x, 3))))
transfer_score=mean(c(3,4,4))

png(filename = "mosquitoes_score_no_sister_num.png", width = 600, height = 600, units = "px", pointsize = 12, bg = "white")
boxplot(data, names=alphas,
	main="Mean divergence between donors and recipients\nin Gambiae complexe phylogeny (Fontaine et al. 2015)",
	xlab="Alpha", ylab="Simulated mean divergence between donors and recipients")

legend("topright", legend = c(paste("observed = ", transfer_score, sep="")),
	col = c("red", "blue"),
	pch = c(19,19), pt.cex = 2, text.col = "black")
abline(h=transfer_score, col="red")
dev.off()


tr<-read.tree('woodcreepers_pulido_2019')
quatuors <- getquatuors(tr)
plot(tr, show.tip.label=T, plot=F)
xy<-get("last_plot.phylo", envir = .PlotPhyloEnv)
ScenrarioPerQuartet<-apply(quatuors, 1, events.donor.receiver, tr=tr, xy=xy)
data = lapply(alphas, function(x) unlist(replicate(100,alpha_sim(x, 5))))
transfer_score=mean(c(2,2,3,2,4))

png(filename = "woodcreepers_score_no_sister_num.png", width = 600, height = 600, units = "px", pointsize = 12, bg = "white")
boxplot(data, names=alphas,
	main="Mean divergence between donors and recipients\nin woodcreepers phylogeny (Pulido-Santacruz et al. 2020)",
	xlab="Alpha", ylab="Simulated mean divergence between donors and recipients")

legend("topright", legend = c(paste("observed = ", transfer_score, sep="")),
	col = c("red", "blue"),
	pch = c(19,19), pt.cex = 2, text.col = "black")
abline(h=transfer_score, col="red")
dev.off()
