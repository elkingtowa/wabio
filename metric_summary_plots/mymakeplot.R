mymakeplot<-function(pname,pvar,pvar2,ptype='density',xlabel,ylabel=NULL,axistxt=NULL,axistxt2=NULL,mlabel='Diversity Panel:69 samples',height=480,width=480,pointsize=12,logy=T,pdir){
	mycols<-c('black','blue','red','magenta','brown','orange')
	print(file.path(pdir,pname))
	png(file.path(plotdir,pname),height=height,width=width,pointsize=pointsize)

	if (ptype=='density')
		plot(density(pvar),xlab=xlabel,main=mlabel)

	if (ptype=='multidensity') {
		yrange<-xrange<-NULL
		for (i in 1:length(pvar)) {
			yrange<-c(yrange,density(pvar[[i]])$y)
			xrange<-c(xrange,density(pvar[[i]])$x)
		}
		plot(density(pvar[[1]]),ylim=c(min(yrange),max(yrange)),xlim=c(min(xrange),max(xrange)),xlab=xlabel,main=mlabel)
		for (i in 2:length(pvar)) {
			lines(density(pvar[[i]]),col=mycols[i])
		}
		legend('top',legend=axistxt,text.col=mycols[1:length(pvar)],bty='n')
	}

	if (ptype=='boxplot' || ptype=='boxplot_errbar') {
		if (ptype=='boxplot_errbar')
			par(mar=c(5,4,4,4)+0.3)
		if (logy) {
			boxplot(pvar,col=mycols[1:length(pvar)],log='y',yaxt='n',xaxt='n',xlab=xlabel, ylab=ylabel,main=mlabel)
			axis(2,c(1,10,100,1000,10000,100000,1000000,10000000))
		} else {
			boxplot(pvar,col=mycols[1:length(pvar)],xaxt='n',xlab=xlabel, ylab=ylabel,main=mlabel)
		}
		mtext(axistxt,side=1,line=1,at=seq(1,length(pvar)))
	}

	if (ptype=='boxplot_errbar') {
		nov.frac.med<-nov.frac.se<-NULL
		for (i in 1:length(pvar2)) {
			nov.frac.med<-c(nov.frac.med,median(pvar2[[i]]))
			nov.frac.se<-c(nov.frac.se, sd(pvar2[[i]])/sqrt(length(pvar2[[i]])))
		}
		par(new=T)
		x<-c(1,seq(2.01,(length(pvar2)+0.01)))
		plot(x,nov.frac.med, type="b", axes=F, bty="n", xlab="", ylab="",col='green',xlim=c(0.5,4.5))
		par(new=T)
		arrows(x,nov.frac.med+nov.frac.se,x,nov.frac.med-nov.frac.se,length=.05,angle=90,code=3,col='green')
		axis(4, at=pretty(range(nov.frac.med)))
		mtext(axistxt2,side=4,line=2)
	}


	graphics.off() 
}
