do.qc <-
function(ff,cols=NULL,type=NULL,out.dir='.') {
	require(stringr)
	cat('Processing',description(ff)$GUID,',',nrow(ff),'events found\n')
	# channels
	if(!is.null(cols)) {
		if(!all(cols %in% parameters(ff)$name)) {
			warning("Some column names are not found")
		}
		cls <- which(parameters(ff)$name %in% cols)
	} else if(!is.null(type)) {
		if(!all(type %in% parameters(ff)$type)) {
			warning("Some types are not found")
		}
		if('CyTOF' %in% type) {
			warning("CyTOF columns are included by default")
			type <- type[type!='CyTOF']
		}
		cls <- which(parameters(ff)$type %in% type)
	} else {
		stop("You must provide either a list of columns or a list of types to transform")
	}
	
	# qc
	qc.df <- data.frame(exprs(ff)[,cls])
	# replace DNA with isotopes / tags
	if(any(str_detect(na.omit(pData(parameters(ff))$isotope),'Rh103'))) {
		names(qc.df)[names(qc.df)==with(pData(parameters(ff)),name[isotope=='Rh103' & !is.na(isotope)])] <- 'Dead'
	} else {
		names(qc.df)[names(qc.df)==with(pData(parameters(ff)),name[isotope=='Pt196' & !is.na(isotope)])] <- 'Dead'
	}
	names(qc.df)[names(qc.df)==with(pData(parameters(ff)),name[isotope=='Ir191' & !is.na(isotope)])] <- 'Ir191'
	names(qc.df)[names(qc.df)==with(pData(parameters(ff)),name[isotope=='Ir193' & !is.na(isotope)])] <- 'Ir193'
	qc.df$Live <- apply(qc.df[,c('Ir191','Ir193')],1,max)
	# extract & reformat time
	if(any(str_detect(pData(parameters(ff))$name,'Leading_Push'))) {
		sec.push <- exprs(ff)[,'Leading_Push']
	} else {
		sec.push <- exprs(ff)[,'Time']*76.8
	}
	sec.push <- sec.push*13e-6
	max.push <- ceiling(max(sec.push))
	seq.push <- seq(0,max.push,by=1)
	bin.push <- cut(sec.push,breaks=seq.push,labels=F)
	pushes <- table(bin.push)
	# plot
	outfile <- paste(substr(description(ff)$GUID,1,nchar(description(ff)$GUID)-4),'png',sep='.')
	png(file.path(out.dir,outfile),width=960,height=1440)
	par(mfrow=c(3,2))
	# event_length distribution
	hist(exprs(ff)[,'Cell_length'],breaks=seq(10,75,by=1),xlab='Number of Pushes',main='Event Length Distribution')
	# number of cells detected by time
	plot(pushes,type='p',axes=F,xlab="Time (s)",ylab="Number of Events")
	axis(1,at=seq(0,max.push,by=30),las=2)
	axis(2,at=10*round(seq(0,max(pushes),length.out=5)/10,0))
	# DNA channels
	plot(apply(qc.df[,c('Ir191','Ir193')],2,asinh),main='Live DNA Channels') # best
	# 
	hist(asinh(qc.df$Dead),
			xlab='Dead DNA Marker',
			main='Dead Channel',
			sub=paste(sum(asinh(qc.df$Dead)<=1),
					'cells out of',
					nrow(qc.df),
					'are considered live (',
					round(100*sum(asinh(qc.df$Dead)<=1)/nrow(qc.df),0),
					'% )')
	)
	abline(v=c(1),col=2,lty=2)
	plot(apply(qc.df[,c('Live','Dead')],2,asinh),main='Live vs Dead DNA Channels') # best
	boxplot(qc.df,las=2)
	dev.off()
	return(invisible(NULL))
}