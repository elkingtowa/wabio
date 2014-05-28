analysis.summary<-function(makeplots=1,mydir=getwd()) {

	setwd(mydir)

#################setup some relevant directories###################
	codedir<-file.path(mydir,'CODE')
	plotdir<-file.path(mydir,'PLOTS')
	system(paste('mkdir -p',plotdir,sep=' ')) 


#################source some modules###################
	source.modules<-c('ancestry.LUT.R',
			'getfiles.R',
			'make.summarymatrix.R',
			'smallvariants.by.ancestry.R',
			'mymakeplot.R')

	for (s in 1:length(source.modules))
		source(file.path(codedir,source.modules[s]))

#################some simple functions#################
	myfunc<-function(myvar){
		id<-grep(myvar,rownames(summary.mat),ignore.case=T)
		if (!is.null(id))
			return(as.numeric(as.character(unlist(summary.mat[id,]))))
	}



#################load ancestry information#################
	ancestry<-ancestry.LUT()


#################read in the filelist for summary file locations#################
	files<-getfiles(ancestry,37)


#################generate the summary matrix#################
##combine summary data across all files in 'files'(above)
	summary.mat<-make.summarymatrix(files,ancestry)



#################Now generate distributions from summary level metrics across all files #################
########A: Extract/Plot gross mapping yield metric ########
	gross.mapping.yield<-myfunc('gross mapping yield')/2.86
	if (makeplots) mymakeplot(pname='gross_mapping_yield.png',pvar=gross.mapping.yield,xlabel='Gross coverage: unique+non-unique mapped bases',pdir=plotdir)

########B: Extract/Plot metrics for mapping of both mates of the DNB ########
##Note: This information is only available for v1.11+
	gen.ws.ge5x<-gen.ws.ge10x<-gen.ws.ge20x<-gen.ws.ge30x<-gen.ws.ge40x<-NULL
	tmp<-myfunc('genome fraction where weightSumSequenceCoverage>=5x')
	if (length(tmp)>0) gen.ws.ge5x<-tmp
	tmp<-myfunc('genome fraction where weightSumSequenceCoverage>=10x')
	if (length(tmp)>0) gen.ws.ge10x<-tmp
	tmp<-myfunc('genome fraction where weightSumSequenceCoverage>=20x')
	if (length(tmp)>0) gen.ws.ge20x<-tmp
	tmp<-myfunc('genome fraction where weightSumSequenceCoverage>=30x')
	if (length(tmp)>0) gen.ws.ge30x<-tmp
	tmp<-myfunc('genome fraction where weightSumSequenceCoverage>=40x')
	if (length(tmp)>0) gen.ws.ge40x<-tmp

#################
########C: Extract/Plot coverage variability metric ########
	cov.var<-myfunc('coverage variability')
	if (makeplots) mymakeplot(pname='coverage_variability.png',pvar=cov.var,xlabel='100k normalized coverage variability across samples',pdir=plotdir)

#################
########D: Extract/Plot genome fraction called metrics ########
##fc: fully-called; pc: partially-called; nc:no-called
	fc.genome.frac<-myfunc('fully called genome fraction')
	pc.genome.frac<-myfunc('partially called genome fraction')
	nc.genome.frac<-myfunc('no-called genome fraction')
	if (makeplots) mymakeplot(pname='fully_called_genome_frac.png',pvar=fc.genome.frac,ptype='density',xlabel='Distribution of fully-called genome fraction',pdir=plotdir)
##not run
	if (0) { 
		if (makeplots) mymakeplot(pname='partially_called_genome_frac.png',pvar=pc.genome.frac,xlabel='Distribution of partially-called genome fraction',pdir=plotdir)
		if (makeplots) mymakeplot(pname='no_called_genome_frac.png',pvar=nc.genome.frac,xlabel='Distribution of no-called genome fraction',pdir=plotdir)
	}

########Global small variant analysis ########
##E: Total Count for small variants
	snpT=myfunc('snp total count')
	smallvarT<-list(snpT=snpT,insT=myfunc('ins total count'),delT=myfunc('del total count'),subT=myfunc('sub total count'))
##F: Novelty rate for small variants
	snpN.frac<-myfunc('snp novel')
	smallvarN.frac<-list(snpN=snpN.frac,insN=myfunc('ins novel'),delN=myfunc('del novel'),subN=myfunc('sub novel'))
	if (makeplots) 
		mymakeplot(pname='smallvariant_count_novrate.png',pvar=smallvarT,pvar2=smallvarN.frac,ptype='boxplot_errbar',xlabel='Distribution of total count(boxplot) and novelty rate(green)',axistxt=c('snp','ins','del','sub'),axistxt2='Novelty rate',pdir=plotdir)


########G: Small variant analysis by poulation types ########
	snpdata<-smallvariants.by.ancestry(vartype='SNP',varT=smallvarT$snpT,varN=smallvarN.frac$snpN,ancestry,makeplots,colnum=2,pdir=plotdir)
	insdata<-smallvariants.by.ancestry(vartype='INS',varT=smallvarT$insT,varN=smallvarN.frac$insN,ancestry,makeplots,colnum=2,pdir=plotdir)
	deldata<-smallvariants.by.ancestry(vartype='DEL',varT=smallvarT$delT,varN=smallvarN.frac$delN,ancestry,makeplots,colnum=2,pdir=plotdir)
	subdata<-smallvariants.by.ancestry(vartype='SUB',varT=smallvarT$subT,varN=smallvarN.frac$subN,ancestry,makeplots,colnum=2,pdir=plotdir)


########H: Small variant het-hom ratio analysis########
	snp.het.hom=myfunc('snp het/hom')
	smallvar.het.hom<-list(snp.het.hom=snp.het.hom,ins.het.hom=myfunc('ins het/hom'),del.het.hom=myfunc('del het/hom'),sub.het.hom=myfunc('sub het/hom'))
	if (makeplots) mymakeplot(pname='smallvariant_het_hom_ratio.png',pvar=smallvar.het.hom,ptype='multidensity',xlabel='Distribution of het-hom ratio',axistxt=c('snp','ins','del','sub'),pdir=plotdir)


#######I: het-hom ratio for all small variants, specifically for ancestries########
	my.bxp<-list(snp.het.hom[which(ancestry[,3]=='Af')],snp.het.hom[which(ancestry[,3]=='Eu')],snp.het.hom[which(ancestry[,3]=='As')],snp.het.hom[which(ancestry[,3]=='Am-admix')])
	if (makeplots) mymakeplot(pname='snp_het_hom_anc.png',pvar=my.bxp,ptype='boxplot',xlabel='Distribution by ancestry',ylabel='SNP het-hom ratio',axistxt=c('African','European','Asian','American'),height=960,width=1000,pointsize=22,logy=F,pdir=plotdir)


########J: Comparison of SNP counts: total vs. het vs hom vs novel########
##homozygous SNP count
##heterozygous SNP count
##novel SNP count
	snpH<-(snp.het.hom*snpT)/(1+snp.het.hom)
	snph<-snpT-snpH
	snpN<-(snpN.frac/100)*snpT
	mybxp<-list(snpT,snpH,snph,snpN)
	if (makeplots) mymakeplot(pname='smallvariant_count.png',pvar=mybxp,ptype='boxplot',xlabel='Distribution of SNP counts',ylabel='Total Count(log scale)',axistxt=c('Total snp','Heterozygous','Homozygous','Novel'),pdir=plotdir)






#######K: genome-wide transition/transversion ratio########
	snp.transition.transversion<-myfunc('transitions/transversions')
	if (makeplots) mymakeplot(pname='snp_transition_transversion.png',pvar=snp.transition.transversion,xlabel='Distribution of SNP transition/transversion ratio',pdir=plotdir)


#######L: Nonsyn/syn SNP ratio########
	snp.nonsyn.syn<-myfunc('Nonsyn/syn')
	if (makeplots) mymakeplot(pname='snp_nonsyn_syn.png',pvar=snp.nonsyn.syn,xlabel='Distribution of SNP NonSynonymous/Synonymous ratio',pdir=plotdir)


#######M: genome-wide/coding insertion/deletions ratio########
	ins.del.all<-list(ins.del=myfunc('insertion/deletions'),ins.del.coding=myfunc('Coding insertion/deletions ratio'))
	if (makeplots) mymakeplot(pname='ins_del_ratio.png',pvar=ins.del.all,ptype='multidensity',xlabel='Distribution of insertions/deletions ratio',axistxt=c('ins.del_genomewide','ins.del_coding'),pdir=plotdir)


#######N: genome-wide frame shifting to preserving ratio########
	shift.preserve.ratio<-myfunc('Frame-shifting/preserving ratio')
	if (makeplots) mymakeplot(pname='shift_preserve.png',pvar=shift.preserve.ratio,xlabel='Distribution of Frame shifting to preserving ratio',pdir=plotdir)



#######O: A bit about libraries########
##A: Mate distribution mean
	mate.dist.mean<-myfunc('Mate distribution mean')
	if (makeplots) mymakeplot(pname='mate_distribution_mean.png',pvar=mate.dist.mean,xlabel='Distribution of mean mate gaps(basepairs)',pdir=plotdir)


	return(list(ancestry=ancestry,summary.mat=summary.mat))
}