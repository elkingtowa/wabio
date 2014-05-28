smallvariants.by.ancestry<-function(vartype,varT,varN,ancestry,makeplots,colnum,pdir) {

	ancestry.id<-function(ancestry,colnum){
		anc.label<-unique(ancestry[,colnum])
		mydata<-NULL

		for(a in 1:length(anc.label)) {
			mydata[a]<-list(grep(anc.label[a],ancestry[,colnum]))
		}

		return(list(anc.label=anc.label,mydata=mydata))
	}


	mysmallvar<-function(ancestry,varT,varN,colnum){
		ancestry.data<-ancestry.id(ancestry,colnum)
		varT.anc<-lapply(seq(1,length(ancestry.data$mydata)),function(i){varT[ancestry.data$mydata[[i]]]})
		varN.frac.anc<-lapply(seq(1,length(ancestry.data$mydata)),function(i){varN[ancestry.data$mydata[[i]]]})
		out<-lapply(seq(1,length(varN.frac.anc)),function(i){return(list(medval=median(varN.frac.anc[[i]]),seval=sd(varN.frac.anc[[i]])/sqrt(length(varN.frac.anc[[i]]))))})
		return(list(ancestry.data=ancestry.data,varT.anc=varT.anc,varN.frac.anc=varN.frac.anc,var.stats=out))	

	}

	vardata<-mysmallvar(ancestry,varT,varN,colnum)

	if (makeplots) {
		##generate a plot showing both distribution of total counts for a given smallvariant type
		##overlay novelty rate on this

		ancestry.data<-vardata$ancestry.data
		varT.anc<-vardata$varT.anc
		varN.frac.anc<-vardata$varN.frac.anc
		varN.frac.anc.med<-unlist(vardata$var.stats)[grep('medval',names(unlist(vardata$var.stats)))]
		varN.frac.anc.se<-unlist(vardata$var.stats)[grep('seval',names(unlist(vardata$var.stats)))]

		my.len<-length(ancestry.data$anc.label)
		anc.col<-rep(1,my.len)	
		if (my.len==12) 
			anc.col<-c(2,2,2,4,5,4,4,4,4,6,6,4)
	
		print(file.path(pdir,paste(vartype,'count_novrate.png',sep='_')))
		png(file.path(pdir,paste(vartype,'count_novrate.png',sep='_')),width=800)
		par(mar=c(5,4,4,4)+0.3)
		##1: first the boxplot of counts
		bx.p<-boxplot(varT.anc,plot=F)
		bxp(bx.p,boxlwd=0.2,xaxt='n',main='Diversity Panel: 69 samples',ylab=paste('Total',vartype,'count',sep=' '),boxfill=anc.col,xlab='Distribution of total count(boxplot) and novelty rate(green)',main='Diversity Panel: small variants')
		mtext(side=1,ancestry.data$anc.label,at=seq(1,length(ancestry.data$anc.label)))
		legend('topright',legend=c('African','American-admixture','European','Asian'),text.col=c(2,4,5,6),bty='n')
		##2: novelty rate w/ error bars
		par(new=T)
		x<-c(1,seq(2.01,(my.len+0.01)))
		plot(x,varN.frac.anc.med, type="b", axes=F, bty="n", xlab="", ylab="",col='green',xlim=c(0.5,(my.len+0.5)))
		par(new=T)
		arrows(x,varN.frac.anc.med+varN.frac.anc.se,x,varN.frac.anc.med-varN.frac.anc.se,length=.05,angle=90,code=3,col='green')
		axis(4, at=pretty(range(varN.frac.anc.med)))
		mtext('Novelty rate',side=4,line=2)
		graphics.off()
	}

	return(vardata)
}