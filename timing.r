#Adapted from Elizabeth Purdom
#Timing for ICGC- PCAWG11
#Pavana Anur

#Load packages
library(cancerTiming)
library(GenomicRanges)
library(gplots)

#set working directory
setwd("./pcawg_timing") 

#########################################################################################
#Load functions
#########################################################################################
#First 5 columns of segData must have column names 'chr','end','start','strand',and 'width'
#"names of mutData must include 'chr' and 'position'

mut2Seg<-function(mutData,segData,verbose=TRUE){
	require(GenomicRanges)
	require(cancerTiming)
	if(any(sort(colnames(segData)[1:5])!=c('chr','end','start','strand','width'))) stop("First 5 columns of segData m
ust have column names 'chr','end','start','strand',and 'width'")
	if(!all(c("chr","position")%in%colnames(mutData))) stop("names of mutData must include 'chr' and 'position'")
	
	#check same chromosome names:
	segChr<-unique(segData$chr)
	mutChr<-unique(mutData$chr)
	newChr<-union(segChr,mutChr)
	if(!all(newChr%in%intersect(segChr,mutChr))){
		#warning("Not all of the chromosomes shared in mutData and segData, will allow the union of them")
		if(any(!newChr%in%segChr)){
			if(verbose){
				cat("Chromosomes not in segData:\n")
				print(newChr[!newChr%in%segChr])
			}
		}
		if(any(!newChr%in%mutChr)){
			if(verbose){
				cat("Chromosomes not in mutData:\n")
				print(newChr[!newChr%in%mutChr])
			}
		}
	}	
	segGr<-GRanges(seqnames =factor(segData$chr,levels=newChr),
	          ranges = IRanges(segData$start, end=segData$end),strand="*",
	          segData[,6:ncol(segData)])
	mutGr<-GRanges(seqnames=factor(mutData$chr, levels=newChr), ranges=IRanges(mutData$position,end=mutData$position)
)
	ov<-findOverlaps(segGr,mutGr)
	combDf<-data.frame(mutData[subjectHits(ov),],
		seg_start=start(segGr)[queryHits(ov)],seg_end=end(segGr)[queryHits(ov)],
			as.data.frame(values(segGr)[queryHits(ov),]))

	#check that none match more than 1
	if(length(unique(subjectHits(ov)))!=length(subjectHits(ov))){
		isDup<-duplicated(subjectHits(ov))
		whDup<-which(subjectHits(ov)%in%subjectHits(ov)[which(isDup)])
		ndups<-table(subjectHits(ov)[whDup])
		if(verbose){
			cat("Overlapping Segments with Mutations matching:\n")
			print(segData[unique(queryHits(ov)[whDup]),])
		}
#		stop(length(ndups)," mutations matched more than 1 segment")
		
	}
	#add those with no hits!
	nmissing<-nrow(mutData)-length(subjectHits(ov))
	if(nmissing>0){
		dummyData<-data.frame(seg_start=NA,seg_end=NA,matrix(NA,ncol=ncol(values(segGr)),nrow=nmissing))
		names(dummyData)<-c("seg_start","seg_end",colnames(values(segGr)))
		if(length(subjectHits(ov))>0) missDat<-cbind(mutData[-subjectHits(ov),],dummyData)
		else{
			missDat<-cbind(mutData,dummyData) #means there were no mutations in these segments!
			if(verbose) cat("there was no overlap between the mutations and the segments\n")
		}
		colnames(missDat)<-make.names(colnames(missDat))
		colnames(combDf)<-make.names(colnames(combDf))
		combDf<-rbind(combDf,missDat)
	}
	combDf<-combDf[order(numChromosome(combDf$chr),combDf$position),]
	return(combDf)
}

getpi<-function(estList){
	pi0<-sapply(unlist(estList,recursive=FALSE),function(x){x$pi["Stage0"]})
	if(length(pi0)>0){
		nam<-sapply(strsplit(names(unlist(estList,recursive=FALSE)),":"),.subset2,1)	
		chr<-sapply(strsplit(nam,"[.]"),.subset2,2)
		type<-sapply(strsplit(nam,"[.]"),.subset2,1)
		nam<-sapply(strsplit(names(unlist(estList,recursive=FALSE)),"[.]"),.subset2,2)	
		N<-sapply(unlist(estList,recursive=FALSE),function(x){x$summaryTable[2,]})
		ui<-sapply(unlist(estList,recursive=FALSE),function(x){x$piCI["Stage0",2]})
		li<-sapply(unlist(estList,recursive=FALSE),function(x){x$piCI["Stage0",1]})
		out<-data.frame(pi0=pi0,lCI=li,uCI=ui,N=N,type=type,chr=chr,ascatId=nam)
		row.names(out)<-NULL
		return(out)
	}
	else return(data.frame(pi0=NA,lCI=NA,uCI=NA,N=NA,type=NA,chr=NA,id=NA,locProb=NA))
}

runEvent<-function(eventArgs){
	mapply(whCall,normCont,allData,FUN=function(x,nc,dat){
		ACNLOH<-makeEventHistory(totalCopy=2,type="LOH")[[1]]
		eventCNLOH<-lapply(x[["CNLOH"]],function(ascatId){
			subdat<-dat[which(dat$ascatId==ascatId),]
			#print(subdat$position)
	#		print(dim(subdat))
			out<-do.call(eventTiming,c(list(x=subdat$t_alt_count, m=subdat$t_depth, 
				history=ACNLOH,totalCopy=2,returnAssignments=T,mutationId=paste(subdat$ascatId,subdat$pos
ition,sep=":"),type="CNLOH",
				normCont=nc),eventArgs))
		})
		names(eventCNLOH)<-x[["CNLOH"]]
		
		AGain<-makeEventHistory(totalCopy=3,type="gain")[[1]]
		eventGain<-lapply(x[["SingleGain"]],function(ascatId){
			subdat<-dat[which(dat$ascatId==ascatId),]
			print(ascatId)
			# print(dim(subdat))
			do.call(eventTiming,c(list(x=subdat$t_alt_count, m=subdat$t_depth, 
				history=AGain,totalCopy=3,returnAssignments=T,mutationId=paste(subdat$ascatId,subdat$posi
tion,sep=":"),type="gain",
				normCont=nc),eventArgs))
		})
		names(eventGain)<-x[["SingleGain"]]
		
		DGain<-makeEventHistory(totalCopy=4,type="gain")[[1]]
		eventDoubleGain<-lapply(x[["DoubleGain"]],function(ascatId){
			subdat<-dat[which(dat$ascatId==ascatId),]
			print(ascatId)
			# print(dim(subdat))
			do.call(eventTiming,c(list(x=subdat$t_alt_count, m=subdat$t_depth, 
				history=DGain,totalCopy=4,returnAssignments=T,mutationId=paste(subdat$ascatId,subdat$posi
tion,sep=":"),type="gain",
				normCont=nc),eventArgs))
		})
		names(eventDoubleGain)<-x[["DoubleGain"]]
		
		return(list(gains=eventGain,CNLOH=eventCNLOH,doubleGains=eventDoubleGain))
		
	},SIMPLIFY=FALSE)
}

getPerLocation<-function(estList){
ll<-unlist(estList,recursive=FALSE)
nms<-unlist(sapply(estList,names))
whSuccess<-which(sapply(ll,function(x){x$success})) #find out which ones were timed successfully
probs<-mapply(ll[whSuccess],nms[whSuccess],FUN=function(x,nam){
whAF<-sapply(c("mutationId","x","^m$"),grep,x=colnames(x$perLocationProb),perl=TRUE)
if(x$success){
pp<-x$perLocationProb[,-whAF]
assignVal<-apply(pp,1,function(rr){if(any(rr>.8)) return(colnames(pp)[rr>.8]) else return("no assignment")})
return(data.frame(segId=nam,x$perLocationProb[,whAF],"assignedAF"=assignVal)) 
}
else{
return(data.frame(segId=nam,x$perLocationProb[,whAF],"assignedAF"=NA))
}},SIMPLIFY=FALSE)
#return(probs)
out<- do.call("rbind",probs)
row.names(out)<-NULL
return(out)
}

##########################################################################################
#code for running timing on a sample
##########################################################################################
# trailingOnly=TRUE means that only your arguments are returned

args <- commandArgs(trailingOnly = TRUE)
print(args)

uuid<-args[1]

#bring in the mutation data -- must fix the names and locations
sampleNames<-uuid
mutDataFiles<-sampleNames
allMutData<-mapply(mutDataFiles,sampleNames,FUN=function(f,n){
	cat("-------",n,"-------\n")
	dat<-read.table(paste("/home/groups/atlas/anurpa/pcawg_timing/MutationData/",f,"_mutations.txt",sep=""),header=T,
sep="\t",stringsAsFactors=FALSE)
	return(dat)
},SIMPLIFY=FALSE)
names(allMutData)<-sampleNames

#bring in the purity/ploidy information from the ascat run
sampleQuality<-read.table("/home/groups/atlas/anurpa/pcawg_timing/purity/purity_ploidy.txt",header=T,sep="\t")

#bring in the segmentations
ascatSegs<-lapply(sampleNames,function(x){read.table(sprintf("Segmentations/%s_segmentation.txt",x),sep="\t",header=TRUE)
})
names(ascatSegs)<-sampleNames

#merge with the merged segmentation data so know what segment each mutation falls into
mData<-mapply(allMutData,ascatSegs,FUN=mut2Seg,SIMPLIFY=FALSE)

#write combination to file
mapply(mData,names(mData),FUN=function(x,nam){write.table(x,file=sprintf("MutationData/%s_mutDataWithSegmentation.txt",na
m),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)})

##########################################################################################
#Run the timing algorithm for those with at least 20 mutations for CNLOH,single Gain and double gain:
##########################################################################################
mutDataFiles<-sampleNames
allData<-mapply(mutDataFiles,sampleNames,FUN=function(f,n){
	cat("-------",n,"-------\n")
	dat<-read.table(paste("MutationData/",f,"_mutDataWithSegmentation.txt",sep=""),header=T,sep="",stringsAsFactors=F
ALSE)
	return(dat)
},SIMPLIFY=FALSE)
names(allData)<-sampleNames

for(samp in sampleNames)
{
allData[[samp]]$call<-"Other"
allData[[samp]]$call<-ifelse(allData[[samp]]$c1N==0 & allData[[samp]]$c2N==2,"CNLOH",allData[[samp]]$call)
allData[[samp]]$call<-ifelse(allData[[samp]]$c1N==2 & allData[[samp]]$c2N==0,"CNLOH",allData[[samp]]$call)
allData[[samp]]$call<-ifelse(allData[[samp]]$c1N==1 & allData[[samp]]$c2N==2,"SingleGain",allData[[samp]]$call)
allData[[samp]]$call<-ifelse(allData[[samp]]$c1N==2 & allData[[samp]]$c2N==1,"SingleGain",allData[[samp]]$call)
allData[[samp]]$call<-ifelse(allData[[samp]]$c1N==1 & allData[[samp]]$c2N==1,"Diploid",allData[[samp]]$call)
allData[[samp]]$call<-ifelse(allData[[samp]]$c1N==1 & allData[[samp]]$c2N==3,"DoubleGain",allData[[samp]]$call)
allData[[samp]]$call<-ifelse(allData[[samp]]$c1N==3 & allData[[samp]]$c2N==1,"DoubleGain",allData[[samp]]$call)
}

whCall<-lapply(allData,function(x){tapply(x$ascatId,factor(x$call,levels=c("Other","CNLOH","SingleGain","Diploid","Double
Gain")),unique)})

for(samp in sampleNames){
	normCont<-1-sampleQuality[sampleQuality$sample==samp,"purity"]
	}
	
#Make segmentation out of ids, need to fix this

segs<-do.call("rbind",mapply(allData,names(allData),FUN=function(x,nam){
	chr<-tapply(x$chr[!is.na(x$ascatId)],x$ascatId[!is.na(x$ascatId)],unique)
	id<-names(chr)
	start<-as.numeric(sapply(strsplit(sapply(strsplit(id,":"),.subset2,2),"-"),.subset2,1))
	end<-as.numeric(sapply(strsplit(sapply(strsplit(id,":"),.subset2,2),"-"),.subset2,2))
	lab<-labelSeg(chr,start,end)
	call<-tapply(x$call,x$ascatId,unique)
	nMut<-tapply(x$ascatId,x$ascatId,length)
	aveReadDepth<-tapply(x$t_depth,x$ascatId,mean)
	medReadDepth<-tapply(x$t_depth,x$ascatId,median)
	return(data.frame(Sample=nam,chr=chr,start=start,end=end,id=id,label=lab,nMut=nMut,call=call,aveRD=aveReadDepth,m
edRD=medReadDepth))
},SIMPLIFY=FALSE))
row.names(segs)<-NULL
segs$width<-segs$end-segs$start+1
segs$mutRateMb<-segs$nMut/segs$width*1e6

#write.table(segs,file=paste("Timings/segsForTiming/",uuid,"_segsForTiming.txt",sep=""),sep="\t",col.names=TRUE,row.names
=FALSE,quote=FALSE)

#Standard MLE

eventArgs<-list(method="fullMLE",seqError=0,bootstrapCI="parametric",B=500,minMutations=20,coverageCutoff=10)
mleEventEst<-runEvent(eventArgs)
xmle<-lapply(mleEventEst,getpi)

xmle<-data.frame(Sample=rep(names(xmle),times=sapply(xmle,nrow)),do.call(rbind,xmle))
write.table(xmle,file=paste("Timings/mleEst/",uuid,"_timingsMle.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE)

save(mleEventEst,file=paste("Timings/rdata/",uuid,".RData"))

##########################################################################################
#Plot
##########################################################################################

seg_mle<-merge(segs,xmle,by.x=c("id"),by.y=c("ascatId"))
seg_mle<-seg_mle[,c(-13,-19)]
#write.table(seg_mle,file=paste("Timings/segMle/",uuid,"_seg_mle.txt",spe=""),row.names=FALSE,col.names=TRUE,quote=F,sep=
"\t")

options(stringsAsFactors=FALSE)
xmle<-na.omit(xmle)
xord<-order(xmle$pi0)
n<-length(xord)

pdf(paste("Timings/plots/",uuid,".pdf",sep=""),width=10,height=6)
plotCI(xmle$pi0[xord],ui=xmle$uCI[xord],li=xmle$lCI[xord],xlab="Chromosome",
	ylab=expression(pi[0]),xaxt="n",sfrac=0.001,cex=.9,labels=xmle$N[xord],pch=c(" "),
	main="CI for segments")
sm<-match(xmle$ascatId[xord],segs$id)
axis(1,1:n,paste(unname(segs$chr[sm]),unname(segs$label[sm]),sep=""),las=2,cex.axis=0.7)	
dev.off()
##########################################################################################
#get mle est per mutation
##########################################################################################

mutInfo<-lapply(mleEventEst,getPerLocation)
mutInfo<-data.frame(Sample=rep(names(mleEventEst),times=sapply(mutInfo,nrow)),do.call(rbind,mutInfo))
write.table(mutInfo,file=paste("Timings/mutInfo/",uuid,"_mutInfo.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE,col.na
mes=TRUE)

#Multiple entries running into error
mutInfo_mle<-merge(mutInfo,xmle,by.x=c("Sample","segId"),by.y=c("Sample","ascatId"))
mutInfo_mle<-mutInfo_mle[,c("Sample","segId","x","m","assignedAF","pi0","type")]
write.table(mutInfo,file=paste("Timings/perLocation/",uuid,"_mut_mle.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE,co
l.names=TRUE)

##########################################################################################







