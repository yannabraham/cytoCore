correct.range <-
function(ff) {
	parameters(ff)$minRange <- apply(exprs(ff),2,min)
	parameters(ff)$maxRange <- apply(exprs(ff),2,max)
	return(ff)
}
