do.deconvolution <-
function(x,barcodes=NULL,live=NULL,dest='.',diagnostics=T,add=T) {
	do.bin <- function(x) {
		sum(x*2^seq(length(x)-1,0))
	}
	
	if(!require(stringr)) {
		stop('the stringr library must be installed\n')
	}
	
	if(diagnostics) {
		if(!file.exists(dest)) {
			dir.create(dest,recursive=T)
		}
	}
	
	if(is.null(live)) {
		live <- TRUE
	}
	
	if(!all(is.logical(live))) {
		stop('live must be a logical vector with TRUE and FALSE as live and dead cells\n')
	}
	
	# handle both data.frames and flowFrames
	if(class(x)=='data.frame') {
		df <- x
	} else if(class(x)=='flowFrame') {
		df <- as.data.frame(exprs(x))
		names(df) <- parameters(x)$name
	} else {
		stop('X must be a data.frame or a flow frame (was ',class(x),')\n')
	}
	
	# try to find barcodes by name if they are not provided
	if(is.null(barcodes)) {
		if(class(x)=='flowFrame') {
			if('type' %in% varLabels(parameters(x))) {
				if('barcode' %in% tolower(parameters(x)$type)) {
					barcodes <- sort(parameters(x)$name[parameters(x)$type=='barcode'])
				}
			}
		} else if(any(str_detect(names(df),'BC'))) {
			barcodes <- sort(names(df)[str_detect(names(df),'BC')])
			cat(length(barcodes),'barcode columns found\n')
			cat(paste('\t',head(barcodes),collapse='\n'),'\n')
		} else {
			stop('barcode columns not defined (either set them explicitly or prefix with BC)\n')
		}
	}
	
	nBC <- ceiling(length(barcodes)^0.5)
	
	if(diagnostics) {
		png(file.path(dest,'histogram.png'),960,960)
		par("mfrow"=c(nBC,nBC))
		hist.live.df <- lapply(barcodes,function(bc) {
					bcdf <- subset(df,live,select=bc,drop=T)
					hist(bcdf,main=bc)
				}
		)
		names(hist.live.df) <- barcodes
		dev.off()
	}
	
	# find thresholds
	dens.live.df <- lapply(barcodes,function(bc) {
				val <- subset(df,live,select=bc,drop=T)
				val <- val[val!=0]
				density(val,bw=0.3)
			}
	)
	names(dens.live.df) <- barcodes
	
	if(diagnostics) {
		png(file.path(dest,'density.png'),960,960)
		par("mfrow"=c(nBC,nBC))
		ksink <- lapply(dens.live.df,plot )
		dev.off()
	}
	
	if(diagnostics) {
		png(file.path(dest,'density_threshold.png'),960,960)
		par("mfrow"=c(nBC,nBC))
	}
	
	inf.live.df <- lapply(dens.live.df,function(ds) {
				inflection <- str_locate(with(ds,paste(1+sign(diff(y)),collapse='')),'02')[,'start']+1
				if(diagnostics) {
					plot(ds)
					with(ds,abline(v=x[inflection],col=2))
				}
				return(ds$x[inflection])
			}
	)
	
	if(diagnostics) {
		dev.off()
	}
	
	if(diagnostics) {
		png(file.path(dest,'histogram_threshold.png'),960,960)
		par("mfrow"=c(nBC,nBC))
		ksink <- lapply(barcodes,function(x) {
					plot(hist.live.df[[x]])
					ref <- hist.live.df[[1]]$breaks[sum(hist.live.df[[1]]$breaks<=inf.live.df[[1]])]
					abline(v=ref,col=2)
				}
		)
		dev.off()
	}
	
	bc.df <- lapply(names(inf.live.df),function(x) { as.numeric(df[,x]>=inf.live.df[[x]]) })
	names(bc.df) <- names(inf.live.df)
	bc.df <- do.call('cbind',bc.df)
	bc.df <- as.data.frame(bc.df)
	names(bc.df) <- str_replace(names(bc.df),'Di','bc')
	
	int2bin <- expand.grid(lapply(barcodes,function(x) return(c(0,1))))
	int2bin <- data.frame(
			barcode=apply(int2bin,1,do.bin),
			code=apply(int2bin,1,paste,collapse='')
	)
	int2bin <- int2bin[order(int2bin$barcode),]
	
	bc.df$barcode <- apply(bc.df,1,do.bin)
	bc.df$barcode <- factor(bc.df$barcode,levels=int2bin$barcode,labels=int2bin$code)
	
	df <- cbind(df,bc.df)
	
	if(diagnostics) {
		png(file.path(dest,'cell_length_barcode.png'),960,960)
		par('mfrow'=c(length(barcodes),length(barcodes)))
		ksink <- lapply(levels(df$barcode),function(x) {
					df <- subset(df,barcode==x & live)
					hist(df$Cell_length,
							main=paste(x,'(',nrow(df),'events - ',round(100*nrow(df)/sum(df$live),1),'% )'),
							breaks=20
					)
				}
		)
		dev.off()
	}
	
	if(add & class(x)!='flowFrame') {
		return(df)
	} else {
		return(bc.df)
	}
}