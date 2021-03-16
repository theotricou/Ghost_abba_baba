## march, 16 2021
## data used here is produced by function ComputeAll() in dstat-topology.R 

##################
## DEPENDENCIES ##
##################
require(ggplot2)
require(scales)


#########
# PLOTS #
#########

############
## READ data

COLNAMES<-c("pext", "alpha", "Nintrog","Nextant","Nreal", "Ntotal","Colless","Nbquartets", paste(rep(seq(0,0.9,by=0.1),4),rep(c("PropIntrogDetected","PropErroneousInterp.mean", "PropErroneousInterp.sd", "PropErroneousInterp.median"), each=10),sep="-"), "ghost.all","noghost.all")
X<-read.table("all_res", sep=" ", col.names=COLNAMES, colClasses=c(rep("numeric",48),rep("character",2)))

############
## FIGURE - effect of distance to outgroup

X2<-X
start<-which(colnames(X2)=="X0.PropErroneousInterp.mean")
stop<-start+9
all<-start:stop
sequence<-seq(0,0.9,by=0.1)
DF<-NULL
for (i in 1:length(all)) {
	print(i)
	DF<-rbind(DF, data.frame(alpha=X2$alpha, thres=rep(sequence[i],nrow(X2)), Nreal=X2$Nreal, Ntotal=X2$Ntotal, PropErroneousInterp.mean=X2[,all[i]]))
}

ggplot(DF, aes(x=thres, y=PropErroneousInterp.mean, group=thres)) + geom_boxplot(aes(fill=factor(alpha)), varwidth=TRUE) + facet_wrap(~alpha, nrow=5) + xlab("Threshold of the relative distance to outgroup considered") + ylab("Proportion of erroneous interpretations of the D-statistics") + labs(fill="Phylogenetic\ndistance\neffect")

############
## FIGURE - same but only for alpha = 0

ggplot(DF[DF$alpha==0,], aes(x=thres, y=PropErroneousInterp.mean, group=thres)) + geom_boxplot(varwidth=TRUE, fill="#F8766D") + xlab(expression(paste("Threshold of the relative distance to outgroup (",italic(R),") considered"))) + ylab("Proportion of erroneous interpretations of the D-statistics")


##############
#### FIGURE  - Complete cases

ggplot(DF, aes(x=factor(Nreal), y=PropErroneousInterp.mean, group=Nreal, fill=factor(alpha)))+geom_boxplot() + facet_grid(alpha~thres) + theme(legend.position = "none") + xlab("Number of extant species before sampling 20 of them") + ylab("Proportion of erroneous interpretation of the D-statistics")

###############
##### FIGURE - zoom

ggplot(DF[DF$alpha==100 & DF$thres==0.5,], aes(x=factor(Nreal), y=PropErroneousInterp.mean, group=Nreal)) +geom_boxplot(fill="#619dff") + theme(legend.position = "none") + xlab("Percentage of extant species considered") + ylab("Proportion of erroneous interpretation of the D-statistics")

###############
##### FIGURE - heatmaps

DF2<-NULL
for (i in unique(DF$alpha)) {
	for (j in unique(DF$thres)) {
		for (k in unique(DF$Nreal)) {
			meanprop<-mean(DF[DF$alpha==i & DF$thres==j & DF$Nreal==k,]$PropErroneousInterp.mean, na.rm=TRUE)
			DF2<-rbind(DF2, data.frame(alpha=i,thres=j,Nreal=k,erroneous=meanprop))
		}
	}
}

ggplot(DF2, aes(x=factor(alpha), y=thres)) + geom_tile(aes(fill=erroneous)) + facet_wrap(~Nreal) +  scale_fill_gradient2(low="blue", high="red",midpoint=0.5, limits=c(0,1), mid="pink") + labs(fill="Mean proportion\nof erroneous interpretation") + xlab("Phylogenetic distance effect on introgressions") + ylab("Threshold of the relative distance to outgroup considered")




