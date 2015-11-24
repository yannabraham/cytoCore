do.transform <-
function(ff,cols=NULL,type=NULL,fun=arcsinhTransform()) {
	if(!is.null(cols)) {
		if(!all(cols %in% parameters(ff)$name)) {
			warning("Some column names are not found")
		}
		cls <- which(parameters(ff)$name %in% cols)
	} else if(!is.null(type)) {
		if(!all(type %in% parameters(ff)$type)) {
			warning("Some types are not found")
		}
		cls <- which(parameters(ff)$type %in% type)
	} else {
		stop("You must provide either a list of columns or a list of types to transform")
	}
	
	exprs(ff)[,cls] <- apply(exprs(ff)[,cls,drop=F],2,fun)
	ff <- correct.range(ff)
	return(ff)
}
