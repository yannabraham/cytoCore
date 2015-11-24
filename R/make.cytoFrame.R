make.cytoFrame <-
function(ff,channels) {
	params <- parameters(ff)
	# extract parameters & values
	newParams <- params[which(params$name %in% unlist(channels) | params$desc %in% unlist(channels)),]
	newExprs <- exprs(ff[,params$name %in% unlist(channels) | params$desc %in% unlist(channels)])
	# compute & rename meaningful parameters
	newParams$isotope <- str_extract(newParams$name,'\\([A-Z]{1}[a-z]{0,1}[0-9]+\\)')
	newParams$isotope <- str_extract(newParams$name,'\\([A-Z]{1}[a-z]{0,1}[0-9]+\\)')
	newParams$isotope <- str_replace_all(newParams$isotope,'[()]','')
	
	# replace names not in channels by corresponding desc provided they are in desc (name has priority)
	newParams$name[newParams$desc %in% unlist(channels) & !newParams$name %in% unlist(channels)] <- 
			newParams$desc[newParams$desc %in% unlist(channels) & !newParams$name %in% unlist(channels)]
	
	# add description
	cdf <- lapply(names(channels),function(x) data.frame(type=x,name=channels[[x]]) )
	cdf <- do.call('rbind',cdf)
	rownames(cdf) <- cdf$name
	cdf$type <- factor(cdf$type,levels=names(channels))
	newParams$type <- cdf[newParams$name,'type']
	# sanitize names
	newParams$name <- make.names(newParams$name)
	dimnames(newExprs)[[2]] <- newParams$name
	# create & return cytoFrame
	newff <- flowFrame(parameters=newParams,exprs=newExprs)
	return(newff)
}