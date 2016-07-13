#'Converts full text data to stm style document term matrix. Redundant with stm function textProcessor()
#' @param source.dir A path to a folder containing a plain text file for each document.
#' @param string A character string where each element is a full document.
#' @param names A vector of document names. Coerced to character. If unspecified no names will be given. If unspecified and source given, files names will be assigned instead.
ftxt2stmbow.f<-function(source.dir=NULL,string=NULL,names=NULL,lt=1,save.to.disk=F,check.for.saved.output=F,out.dir){

	sfn<-"ftxt2stmbow.RData"
	if(check.for.saved.output) if(any(grepl(sfn,dir(recursive=T,full.names=T,ignore.case=T)))) {
		warning(paste('Loading and returning first saved',sfn),call.=F)
		l<-ls()
		load(dir(pattern=sfn,full.names=T,recursive=T,ignore.case=F)[1])
		return(get(setdiff(ls(),c(l,'l'))))
	}

	if(all(sapply(list(source.dir,string),is.null))) stop('Specify either source.dir or string.')
	if(all(!sapply(list(source.dir,string),is.null))) stop('Specify either source.dir or string but not both.')
	require(stm,quietly = T)
	# 	require(tm,quietly = T)
	# 	require(SnowballC,quietly = T)
	require(data.table,quietly = T)

	# 1. Preprocessing functions in base R

	if(!is.null(source.dir)) {
		docs<-list() # container for docs
		files<-list.files(source.dir,full.names=T)
		if(!is.null(sample.docs)) files<-sample(files,sample.docs)
		files<-files[order(sub(paste('.*','(.+)',sep=.Platform$file.sep),'',files))] # helps later to have files in alpha order by document name
		for(i in files) docs[[i]]<-readLines(i,warn=F)
		docs<-lapply(docs,FUN=paste,collapse=' ') # each doc is one long character string
		string<-sapply(docs,FUN=strsplit,split='\\s') # split docs into words "\\s+" is a regex. "[[:space:]]+" also works but is R specific
		names(string)<-names
	}

	if(!is.null(string)) {
		docs<-string
		string<-sapply(docs,FUN=strsplit,split='\\s') # split docs into words "\\s+" is a regex. "[[:space:]]+" also works but is R specific
		names(string)<-names
	}

	docs<-lapply(docs,FUN=tolower) # transform to lower case
	docs<-lapply(docs,FUN=removePunctuation) # ...
	docs<-lapply(docs,FUN=removeNumbers)
	docs<-sapply(docs,FUN=strsplit,split='\\s') # split docs into words "\\s+" is a regex. "[[:space:]]+" also works but is R specific
	docs<-lapply(docs,FUN=removeWords,stopwords('english'))
	docs<-lapply(docs,FUN=stemDocument,language='english')

	## here we end with a list of untabulated tokenized character vectors in original document order, with blanks preserved.

	# 2. Match tokens to originals for applying stm results for qualitative cross validation.
	ftxt2stmbow<-list()
	ftxt2stmbow$map<-mapply(function(o,t) {t<-c(t,rep('',length(o)-length(t)));data.table(o,t,ord=1:length(o))},o=string,t=docs,SIMPLIFY = F) #original and tokens

	docs<-lapply(docs,FUN=function(x) x[!!nchar(x)]) #remove blanks

	# 3. Sparse matrix format expected by stm.

	ftxt2stmbow$vocab<-sort(unique(unlist(docs))) # stm expects as input a list of matrices where each element of the list is a document and where the first row of each matrix is the index of a word and the second row is the count of the word. This is a more memory efficient form since zeros are not stored.
	for(i in 1:length(docs)){
		t<-table(docs[[i]])
		ftxt2stmbow$documents[[i]]<-rbind(vocab.index=which(ftxt2stmbow$vocab%in%names(t)),frequency=t)
	}
	names(ftxt2stmbow$documents)<-names

	p<-prepDocuments(ftxt2stmbow$documents,ftxt2stmbow$vocab,meta = names,lower.thresh=lt)
	ftxt2stmbow$documents<-p$documents
	ftxt2stmbow$vocab<-p$vocab
	names(ftxt2stmbow$documents)<-p$meta

	if(save.to.disk) try(save(ftxt2stmbow,file=paste(grep(paste(out.dir,'$',sep=''),dir(include.dirs = T,full.names=T,recursive=F,ignore.case=F),value = T)[1],sfn,sep=.Platform$file.sep)))

	ftxt2stmbow
}

#'Imports from ICPSR Congressional Record files saved as a CTAWG list, and performs conventional preprocessing to yield an stm formatted term-document "sparse matrix".
#' @param icpsr.cong A CTAWG list of 3 data frames named according to ICPSR plain text tables and containing the "GPOspeech", "SpeakerID", and "GPOdescr" substrings.
icpsr2stmbow.f<-function(icpsr,sample.size=100,wrdmin=0,wrdmax=0,lt=1,save.to.disk=F,check.for.saved.output=F,out.dir=NULL){

	sfn<-paste("icpsr2stmbow-samp",ifelse(sample.size,sample.size,'all'),".RData",sep="")
	if(check.for.saved.output) if(any(grepl(sfn,dir(recursive=T,full.names=T,ignore.case=T)))) {
		warning(paste('Loading and returning first saved',sfn),call.=F)
		l<-ls()
		load(dir(pattern=sfn,full.names=T,recursive=T,ignore.case=F)[1])
		return(get(setdiff(ls(),c(l,'l'))))
	}

	require(stm,quietly = TRUE)
	s<-1:nrow(icpsr[[grep('GPOspeech',names(icpsr))]])
	if(wrdmin) {wn<-icpsr[[grep('GPOdescr',names(icpsr))]]$word.count>=wrdmin} else {wn<-rep(T,length(s))}
	if(wrdmax) {wx<-icpsr[[grep('GPOdescr',names(icpsr))]]$word.count<=wrdmax} else {wx<-rep(T,length(s))}
	s<-s[wn&wx]
	if(sample.size) s<-sample(s,sample.size)
	GPOspeech<-grep('GPOspeech',names(icpsr))
	icpsr2stmbow<-textProcessor(
		documents=icpsr[[GPOspeech]]$speech[s]
		,lowercase=TRUE
		,removestopwords=TRUE
		,removenumbers=TRUE
		,removepunctuation=TRUE
		,stem=TRUE
		,wordLengths=c(3,Inf)
		,sparselevel=1
		,language="en"
		,verbose=TRUE
		,onlycharacter= FALSE
		,striphtml=FALSE
		,customstopwords=NULL)
	nm<-icpsr[[GPOspeech]]$speechID[s]
	if(length(icpsr2stmbow$docs.removed)) nm<-nm[-icpsr2stmbow$docs.removed]
	p<-prepDocuments(icpsr2stmbow$documents,icpsr2stmbow$vocab,meta = nm,lower.thresh = lt)
	icpsr2stmbow$documents<-p$documents
	icpsr2stmbow$vocab<-p$vocab
	names(icpsr2stmbow$documents)<-p$meta
	if(save.to.disk) try(save(icpsr2stmbow,file=paste(grep(paste(out.dir,'$',sep=''),dir(include.dirs = T,full.names=T,recursive=F,ignore.case=F),value = T)[1],sfn,sep=.Platform$file.sep)))
	icpsr2stmbow
}

#' Browse ICPSR records using either speech or speaker IDs.
#' @param icpsr A CTAWG list of 3 data frames named according to ICPSR plain text tables and containing the "GPOspeech", "SpeakerID", and "GPOdescr" substrings.
#' @param speechID The name of the field containing the speech ID. First three digits indicate the nth Congress, next 4 indicate the year, following 7 are the numerical index in order of appearance.
#' @param speakerID Speaker ID, which the function links via the "GPOdesc" table to merge speaker and speech data.
browse.icpsr<-function(icpsr,speechID=NULL,speakerID=NULL,print.max=50){
	require(data.table,quietly = T)
	if(is.null(speechID)&is.null(speakerID)) stop('Must specify either speechID or speakerID.')

	icpsr<-icpsr.cong107
	for(i in grep('GPOspeech',names(icpsr))) icpsr[[i]]$speechID<-as.character(icpsr[[i]]$speechID)
	GPOspeech<-rbindlist(icpsr[grep('GPOspeech',names(icpsr))])
	icpsr[grep('GPOspeech',names(icpsr))]<-NULL
	setkey(GPOspeech,speechID)

	for(i in grep('SpeakerID',names(icpsr))) icpsr[[i]]$id<-as.character(icpsr[[i]]$id)
	SpeakerID<-rbindlist(icpsr[grep('SpeakerID',names(icpsr))])
	icpsr[grep('SpeakerID',names(icpsr))]<-NULL
	setnames(SpeakerID,'id','speakerID')
	setkey(SpeakerID,speakerID)

	for(i in grep('GPOdesc',names(icpsr))) {
		icpsr[[i]]$speechID<-as.character(icpsr[[i]]$speechID)
		icpsr[[i]]$speakerID<-as.character(icpsr[[i]]$speakerID)
	}
	GPOdesc<-rbindlist(icpsr[grep('GPOdesc',names(icpsr))])
	icpsr[grep('GPOdesc',names(icpsr))]<-NULL

	if(!is.null(speechID)){
		speechID<-as.character(speechID)
		setkey(GPOdesc,speechID)
		ret<-merge(GPOdesc,GPOspeech[speechID])
		setkey(ret,speakerID)
		ret<-merge(SpeakerID,ret)
		if(nrow(ret)<print.max) {cat(paste('Printing first',floor(print.max),'documents.\n'));for(i in 1:ifelse(nrow(ret)>print.max,print.max,nrow(ret))) {print(t(ret[i]));cat('\n')}}
		return(ret)
	}
	if(!is.null(speakerID)){
		speakerID<-as.character(speakerID)
		setkey(GPOdesc,speakerID)
		ret<-merge(GPOdesc,SpeakerID[speakerID])
		setkey(GPOdesc,speechID)
		setkey(ret,speechID)
		ret<-merge(ret,GPOspeech)
		cat(paste('Printing first',floor(print.max),'documents.'))
		if(nrow(ret)<print.max) {cat(paste('Printing first',floor(print.max),'documents.\n'));for(i in 1:ifelse(nrow(ret)>print.max,print.max,nrow(ret))) {print(t(ret[i]));cat('\n')}}
		return(ret)
	}
}

#' Takes stm style bag of words as input and gives LDA.
#' @param stmbow stm formatted bag of words, as given by textProcessor command or stmbow importer.
#' @param out.dir: Character path to a file folder where you want your output files to be stored once text processing is completed.
#' @param k Number of topics.
#' @param alpha The alpha parameter must be greater than 0. Alpha < 1 assumes that each document is constructed from relatively few topics. Alpha > 1 assumes that each is constructed from many topics. If you aren't sure choose the convention: 50/k, which will also be the default if you specify nothing.
#' @param sample.docs If not NULL, an integer specifying the size of a sample taken from the original documents. Helpful for debugging.
#' @param visualize.results Logical, if TRUE provides visualization tools to help interpret LDA output.
#' @param verbose Print iteratio progress (very verbose!).
#' @param save.to.disk Save an .RData file of the object; helpful for lengthy estimations. May be buggy across environments.
stmbow2lda.f<-function(stmbow,out.dir,k,alpha=NULL,visualize.results=F,verbose=F,save.to.disk=F,check.for.saved.output=F)
{
	if(is.null(alpha)) alpha<-50/k
	if(check.for.saved.output) if(any(grepl(paste("stm-model-k",k,"-alpha",round(alpha,3),".RData",sep=""),dir(recursive=T,full.names=T,ignore.case=T)))) {
		warning(paste('Loading and returning first saved',paste("stm-model-k",k,"-alpha",round(alpha,3),".RData.",sep="")),call.=F)
		load(dir(pattern=paste("stm-model-k",k,"-alpha",round(alpha,3),".RData",sep=""),full.names=T,recursive=T,ignore.case=F)[1])
		return(stmbow2lda)
	}
	# Check package requirements and arguments
	require(stm,quietly = T)


	stmbow2lda<-list()
	stmbow2lda$model<-stm(documents=stmbow$documents,vocab=stmbow$vocab,K=k,control=list(alpha=alpha),verbose=verbose)
	stmbow2lda$top.word.phi.beta<-sapply(data.frame(stmbow2lda$model$beta$logbeta),function(x) sapply(x,function(y) ifelse(is.infinite(y),.Machine$double.eps,exp(y)))) # called beta by stm, epsilon closest thing to zero the machine can represent, necessary to prevent error
	colnames(stmbow2lda$top.word.phi.beta)<-stmbow2lda$model$vocab
	stmbow2lda$doc.top.theta<-stmbow2lda$model$theta
	rownames(stmbow2lda$doc.top.theta)<-names(stmbow$documents)
	stmbow2lda$doc.length<-sapply(stmbow$documents,ncol)
	stmbow2lda$vocab<-stmbow2lda$model$vocab
	tn<-data.table(do.call(rbind,sapply(stmbow$documents, t)))
	setnames(tn,c('ix','freq'))
	setkey(tn,ix)
	tn<-tn[,list('freq'=sum(freq)),by=ix]
	stmbow2lda$term.frequency<-tn$freq
	names(stmbow2lda$term.frequency)<-stmbow$vocab[tn$ix]

	if(save.to.disk){
		pfs<-.Platform$file.sep
		f<-paste(out.dir,paste("stm-model-k",k,"-alpha",round(alpha,3),".RData",sep=""),sep=pfs)
		save(stmbow2lda,file=f)
	}

	stmbow2lda
}

#' Different approaches to reduce and visualize the standard topic-document and topic-word word matrices output by LDA and other topic modelling estimations.
lda2viz.f<-function(stmbow2lda){
	# from http://cpsievert.github.io/LDAvis/newsgroup/newsgroup.html
	# http://nlp.stanford.edu/events/illvi2014/papers/sievert-illvi2014.pdf
	# http://glimmer.rstudio.com/cpsievert/xkcd/
	require(LDAvis,quietly = T)
	require(servr,quietly = T)

	json <- createJSON(
		phi = stmbow2lda$top.word.phi.beta
		,theta = stmbow2lda$doc.top.theta
		,vocab = stmbow2lda$vocab
		,doc.length = stmbow2lda$doc.length
		,term.frequency = stmbow2lda$term.frequency
		,reorder.topics = F
	)
	save(json,file='viz.RData')
	serVis(json, out.dir = "vis", open.browser = F)
}

#'Network visualizations and clustering.
lda2netviz.f<-function(stmbow2lda,thresh="choose"){

	require(network,quietly = T)

	bam<-stmbow2lda$doc.top.theta
	b<-quantile(bam,seq(0,1,.05))
	h<-hist(bam,breaks=b,col="black")
	abline(v=b,col=gray(0,.5))
	text(
		x=rev(h$breaks)
		,y=seq(0,max(h$density),length.out=21)
		,labels=rev(paste('|',names(b),' <= ',round(b,3),sep=''))
		,pos=4
		,offset=0
		,cex=1
	)
	if(thresh=="choose"){
		cat("\nPlease choose an edge weight threshold by clicking on the histogram near the x-axis where you would like to cut the distribution of probabilities that a document draws words from a particular topic (i.e. theta or the document-topic probability matrix). Relationships between documents and topics that fall below this threshold will be ignored.\n")
		thresh<-locator(n=1,type="p",col="red")
		abline(v=thresh$x,col="red")
		text(
			x=thresh$x
			,y=thresh$y
			,labels=rev(paste('|',round(mean(bam<thresh$x)*100,2),'% <= ',round(thresh$x,3),sep=''))
			,pos=4
			,offset=0
			,cex=1
			,col="red"
		)
	}
	bam[bam<thresh$x]<-0
	m1am<-bam%*%t(bam)
	m1net<-network(m1am,directed=F,loops=F)
	network.vertex.names(m1net)<-sub(paste('.*','(.+)',sep=.Platform$file.sep),'\\1',rownames(m1am))
	#	pdf('doc-by-top-net.pdf')
	plot(m1net
			 ,displaylabels=T
			 ,label=paste('T',1:nrow(m1am),sep='')
			 ,label.pos=5
			 ,label.col="white"
			 ,label.cex=.25
			 ,vertex.col="black"
			 ,vertex.cex=2
	)
	#	dev.off()

	m2am<-t(bam)%*%bam
	m2net<-network(m2am,directed=F,loops=F)

	#	pdf('top-by-doc-net.pdf')
	plot(m2net
			 ,displaylabels=T
			 ,label=paste('D',1:nrow(m2am),sep='')
			 ,label.pos=5
			 ,label.col="black"
			 ,label.cex=.75
			 ,vertex.col="white"
			 ,vertex.cex=2
	)
	#	dev.off()

	# 	b<-quantile(m1am,seq(0,1,.1))
	# 	h<-hist(m1am,breaks=b)
	# 	abline(v=b,col'pink')
	# 	h<-hist(m2am,breaks=quantile(m2am,seq(0,1,.1)))

	bel<-which(bam>0,arr.ind=T)
	w<-bam[bel]
	bel<-cbind(sub(paste(".+","(.+)",sep=.Platform$file.sep),"\\1",rownames(bel)),paste("t",bel[,2],sep=""))
	o<-order(bel[,1],bel[,2])
	bel<-data.frame(bel[o,])
	w<-w[o]
	colnames(bel)<-c("document","topic")
	nm1<-length(unique(bel$document))
	nm2<-length(unique(bel$topic))
	bnet<-network(bel,bipartite=nm1,matrix.type="edgelist")
	bnet%e%"w"<-w
	pdf('bimodal-net.pdf')
	plot(
		bnet
		,displaylabels=T
		,label=c(
			paste(1:nm1)
			,network.vertex.names(bnet)[-(1:nm1)]
		)
		,label.pos=5
		,label.col=c(rep("white",nm1),rep("black",nm2))
		,label.cex=.75
		,vertex.col=c(rep("black",nm1),rep("white",nm2))
		,vertex.cex=2
		#,vertex.sides=c(rep(3,nm1),rep(20,nm2))
	)
	dev.off()


}

#'Color coding original documents in topic highlighting.
lda2ftxt.f<-function(map,doc.top,top.word,lda2rel,intensify=T,intensity=.3,sample=10,out.dir=NULL,index=0,spacing=.4,ptsize=10,axes=F,pdf=F,mxdoc.word.prop=0)
{
	d<-dir(include.dirs = T,recursive=T,full.names = T)
	out.dir<-grep(paste(out.dir,'$',sep=''),d,value=T)
	if(!length(out.dir)) {warning('out.dir not found, saving to current working directory.');out.dir<-getwd()}

	require(colorspace,quietly = T)

	mx<-max(sapply(1:nrow(doc.top), function(x) max(doc.top[x,]*top.word)))
	s<-1:length(map)
	if(sample) s<-sort(sample(s,sample))
	if(index[1]) s<-index
	for(i in s){

		cat(c('\nRendering document',i))
		w<-order(doc.top[i,],decreasing=T)[1:2]
		tn<-doc.top[i,w]
		names(tn)<-paste('T',w,sep='')
		p<-round(sum(tn)*100,1)

		cat('\nCalculating document\'s topic by word probability matrix =\n(document by topic prob vector) * (global topic by word prob matrix)')
		m<-t(doc.top[i,w]*top.word[w,])

		cat('\nOriginal range of predicted document\'s topic by word probabilities:')
		print(range(m))

		if(intensify) {
			m<-lintran(m,c(0,mx),c(intensity,1))
			cat('\nRange of document\'s topic by word probabilities after linear amplification:')
			print(range(m))
		}

		colnames(m)<-names(tn)
		m<-data.table(t=rownames(m),m,r=apply(m,1,function(x) x[2]/sum(x)),w=apply(m,1,function(x) sum(x)))
		setkey(m,t)
		setkey(map[[i]],t)
		m<-merge(map[[i]],m,all.x=T,all.y=F)
		col<-hex(HSV(H=240+120*m$r,S=1,V=m$w))
		col[is.na(m$w)]<-gray(.3)
		m[,'col':=col]
		setkey(m,ord)
		mr<-.6

		par(mar=rep(0,4))

		plot.new()
		pw<-8
		ph<-1
		par(family='Times',ps=ptsize,fin=c(pw,ph))
		plot.window(xlim=c(0,pw),ylim=c(0,1))

		o<-as.vector(rbind(m$o,' '))
		t<-as.vector(rbind(m$t,' '))
		c<-as.vector(rbind(m$col,'#000000'))
		shu<-strheight(LETTERS,units = 'inches')
		swu<-strwidth(o,units = 'inches')*1.25

		fl<-list()
		x<-list()
		sww<-swu
		ct<-0
		f<-cumsum(sww)%/%pw
		while(max(f)>0){
			ft<-f==0
			x[[length(x)+1]]<-cumsum(sww[ft])
			sww<-sww[!ft]
			fl[[length(fl)+1]]<-rep(ct,sum(ft))
			f<-cumsum(sww)%/%pw
			ct<-ct+1
		}
		ft<-f==0
		x[[length(x)+1]]<-cumsum(sww[f==0])
		fl[[length(fl)+1]]<-rep(ct,sum(ft))
		f<-unlist(fl)
		x<-unlist(sapply(x,function(x) c(0,x[-length(x)])))
		l<-max(f)*2+1.5

		lda2relc<-copy(lda2rel)
		setkey(lda2relc,Category)
		lda2relc<-lda2relc[names(tn)]
		lda2relc[,inc:=Term%in%t]

		leg<-sum(c(
			head=4
			,rows=nrow(lda2relc)/length(unique(lda2relc$Category))/length(unique(lda2relc$lambda))
			,foot=0
		))

		ph<-sum(c(lines=l,legend=leg,divider=1))

		par(mar=c(b=2.25,l=2.25,t=3.25,r=2.25),fin=c(pw,ph*spacing),yaxs='i')
		plot.window(xlim=c(0,pw),ylim=c(ph,1))
		box()
		if(axes){
			axis(1)
			axis(2)
		}
		abline(h=leg)

		#Topic Super Column headings
		text(x=seq(0,pw,length.out = 9)[c(3,7)],y=2,labels=names(tn),col=c('Blue','Red'),font=2,pos=1,offset=0)
		#Rank listed
		text(x=seq(0,pw,length.out = 9)[-c(1,5,9)],y=3,labels=c('Anchor','~','Common'),font=3,pos=1,offset=0)

		loc<-data.table(expand.grid(Category=names(tn),lambda=unique(lda2relc$lambda)),col=c(2,6,3,7,4,8))
		setkey(lda2relc,Category,lambda)
		setkey(loc,Category,lambda)
		lda2relc<-merge(lda2relc,loc)
		lda2relc[,clr:=ifelse(inc,'black','gray')]

		# List of terms
		text(x=seq(0,pw,length.out = 9)[lda2relc$col],y=lda2relc$ord+3,labels=lda2relc$Term,col = lda2relc$clr,pos=1,offset=0)
		#Original text
		text(x=x,y=f*2+leg+1,labels=o,pos=4,offset=0,cex=1,col=c,ps=11)
		#Tokenized text
		text(x=x,y=f*2+leg+1.5,labels=t,pos=4,offset=0,cex=1,font=3,col='gray5')
		#Line numbers
		text(x=0,y=unique(f)*2+leg+1.25,labels=unique(f)+1,pos=4,offset=-0.5,cex=.75,col='lightgray',ps=11,srt=90)

		title(main=paste(c('Document',i,':',paste(c('Blue','Red'),names(tn),round(tn*100,2),'%'),paste('=',round(sum(tn)*100,2),'%')),collapse = ' '))

		if(pdf) try(dev.copy2pdf(file=paste(out.dir,paste('doc',i,'.pdf',sep=''),sep=.Platform$file.sep), width=pw, height=ph*spacing),silent = T)

		cat('\nRange of ',ifelse(intensify,'amplified ',''),'probabilites for ',names(m)[[4]],sep='')
		print(range(na.omit(m[[4]])))
		cat('\nRange of ',ifelse(intensify,'amplified ',''),'probabilites for ',names(m)[[5]],sep='')
		print(range(na.omit(m[[5]])))
	}

	m
}

#' Original authors code to compute relevance. 'Forked' from https://github.com/cpsievert/LDAvis/blob/6f93aa85499b705c9ae6c56e5985df637f9f5132/R/createJSON.R
lda2rel.f<- function(stmbow2lda,R = 10,lambda.step = 0.5,reorder.topics = TRUE,save.to.disk=F,check.for.saved.output=F,out.dir,...) {

	sfn<-'lda2rel.RData'
	if(check.for.saved.output) if(any(grepl(sfn,dir(recursive=T,full.names=T,ignore.case=T)))) {
		warning(paste('Loading and returning first saved',sfn),call.=F)
		l<-ls()
		load(dir(pattern=sfn,full.names=T,recursive=T,ignore.case=F)[1])
		return(get(setdiff(ls(),c(l,'l'))))
	}

	require(data.table,quietly = T)
	phi = stmbow2lda$top.word.phi.beta
	theta = stmbow2lda$doc.top.theta
	vocab = stmbow2lda$vocab
	doc.length = stmbow2lda$doc.length
	term.frequency = stmbow2lda$term.frequency
	#rm(stmbow2lda)
	# Set the values of a few summary statistics of the corpus and model:
	dp <- dim(phi)  # should be K x W
	dt <- dim(theta)  # should be D x K

	N <- sum(doc.length)  # number of tokens in the data
	W <- length(vocab)  # number of terms in the vocab
	D <- length(doc.length)  # number of documents in the data
	K <- dt[2]  # number of topics in the model

	# check that certain input dimensions match
	if (dp[1] != K) stop("Number of rows of phi does not match
											 number of columns of theta; both should be equal to the number of topics
											 in the model.")
	if (D != dt[1]) stop("Length of doc.length not equal
											 to the number of rows in theta; both should be equal to the number of
											 documents in the data.")
	if (dp[2] != W) stop("Number of terms in vocabulary does
											 not match the number of columns of phi (where each row of phi is a
											 probability distribution of terms for a given topic).")
	if (length(term.frequency) != W) stop("Length of term.frequency
																				not equal to the number of terms in the vocabulary.")
	if (any(nchar(vocab) == 0)) stop("One or more terms in the vocabulary
																	 has zero characters -- all terms must have at least one character.")

	# check that conditional distributions are normalized:
	phi.test <- all.equal(rowSums(phi), rep(1, K), check.attributes = FALSE)
	theta.test <- all.equal(rowSums(theta), rep(1, dt[1]),
													check.attributes = FALSE)
	if (!isTRUE(phi.test)) stop("Rows of phi don't all sum to 1.")
	if (!isTRUE(theta.test)) stop("Rows of theta don't all sum to 1.")

	# compute counts of tokens across K topics (length-K vector):
	# (this determines the areas of the default topic circles when no term is
	# highlighted)
	topic.frequency <- colSums(theta * doc.length)
	topic.proportion <- topic.frequency/sum(topic.frequency)

	# re-order the K topics in order of decreasing proportion:
	if(reorder.topics) {o <- order(topic.proportion, decreasing = TRUE)} else {o <- seq_along(topic.proportion)}

	phi <- phi[o, ]
	theta <- theta[, o]
	topic.frequency <- topic.frequency[o]
	topic.proportion <- topic.proportion[o]

	# compute intertopic distances using the specified multidimensional
	# scaling method:
	# 	mds.res <- mds.method(phi)
	# 	if (is.matrix(mds.res)) {
	# 		colnames(mds.res) <- c("x", "y")
	# 	} else if (is.data.frame(mds.res)) {
	# 		names(mds.res) <- c("x", "y")
	# 	} else {
	# 		warning("Result of mds.method should be a matrix or data.frame.")
	# 	}
	# 	mds.df <- data.frame(mds.res, topics = seq_len(K), Freq = topic.proportion*100,
	# 											 cluster = 1, stringsAsFactors = FALSE)
	# note: cluster (should?) be deprecated soon.

	# token counts for each term-topic combination (widths of red bars)
	term.topic.frequency <- phi * topic.frequency

	# compute term frequencies as column sums of term.topic.frequency
	# we actually won't use the user-supplied term.frequency vector.
	# the term frequencies won't match the user-supplied frequencies exactly
	# this is a work-around to solve the bug described in Issue #32 on github:
	# https://github.com/cpsievert/LDAvis/issues/32
	term.frequency <- colSums(term.topic.frequency)
	stopifnot(all(term.frequency > 0))

	# marginal distribution over terms (width of blue bars)
	term.proportion <- term.frequency/sum(term.frequency)

	# Old code to adjust term frequencies. Deprecated for now
	# adjust to match term frequencies exactly (get rid of rounding error)
	#err <- as.numeric(term.frequency/colSums(term.topic.frequency))
	# http://stackoverflow.com/questions/3643555/multiply-rows-of-matrix-by-vector
	#term.topic.frequency <- sweep(term.topic.frequency, MARGIN=2, err, `*`)

	# Most operations on phi after this point are across topics
	# R has better facilities for column-wise operations
	phi <- t(phi)

	# compute the distinctiveness and saliency of the terms:
	# this determines the R terms that are displayed when no topic is selected
	topic.given.term <- phi/rowSums(phi)  # (W x K)
	kernel <- topic.given.term * log(sweep(topic.given.term, MARGIN=2,
																				 topic.proportion, `/`))
	distinctiveness <- rowSums(kernel)
	saliency <- term.proportion * distinctiveness

	# Order the terms for the "default" view by decreasing saliency:
	default.terms <- vocab[order(saliency, decreasing = TRUE)][1:R]
	counts <- as.integer(term.frequency[match(default.terms, vocab)])
	Rs <- rev(seq_len(R))
	default <- data.frame(Term = default.terms, logprob = Rs, loglift = Rs,
												Freq = counts, Total = counts, Category = "Default",
												stringsAsFactors = FALSE,lambda=NA)
	topic_seq <- rep(seq_len(K), each = R)
	category <- paste0("Topic", topic_seq)
	lift <- phi/term.proportion

	# Collect R most relevant terms for each topic/lambda combination
	# Note that relevance is re-computed in the browser, so we only need
	# to send each possible term/topic combination to the browser
	find_relevance <- function(i) {
		relevance <- i*log(phi) + (1 - i)*log(lift)
		idx <- apply(relevance, 2,
								 function(x) order(x, decreasing = TRUE)[seq_len(R)])
		# for matrices, we pick out elements by their row/column index
		indices <- cbind(c(idx), topic_seq)
		data.frame(Term = vocab[idx], Category = category,
							 logprob = round(log(phi[indices]), 4),
							 loglift = round(log(lift[indices]), 4),
							 stringsAsFactors = FALSE)
	}
	lambda.seq <- seq(0, 1, by=lambda.step)
	#if (missing(cluster)) {
	tinfo <- lapply(as.list(lambda.seq), function(x) {x<-data.frame(find_relevance(x),lambda=as.character(x));data.frame(x,ord=1:R)})
	#} else {
	#	tinfo <- parallel::parLapply(cluster, as.list(lambda.seq), find_relevance)
	#}
	tinfo <- unique(do.call("rbind", tinfo))
	tinfo$Total <- term.frequency[match(tinfo$Term, vocab)]
	rownames(term.topic.frequency) <- paste0("Topic", seq_len(K))
	colnames(term.topic.frequency) <- vocab
	tinfo$Freq <- term.topic.frequency[as.matrix(tinfo[c("Category", "Term")])]
	#tinfo <- rbind(default, tinfo)
	tinfo$Category<-sub('opic','',tinfo$Category)

	lda2rel<-data.table(tinfo)

	if(save.to.disk) try(save(lda2rel,file=paste(grep(paste(out.dir,'$',sep=''),dir(include.dirs = T,full.names=T,recursive=F,ignore.case=F),value = T)[1],sfn,sep=.Platform$file.sep)))

	lda2rel
}

lintran<-function(x,s1=c(0,1),s2=c(0,1)) {a=diff(s2)/diff(s1);b=s2[1]-a*s1[1];return(a*x+b)}
