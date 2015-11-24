\name{correct.range}
\alias{correct.range}

\title{
Updates the range of the \link[flowCore]{flowFrame} object after a transformation or a filtering 
}

\description{
Following transformation or filtering the range columns in parameters may not reflect the current range anymore. Correct.range() will update the range to the current set of values in each of the columns of the exprs slot of the flowFrame 
}

\usage{
correct.range(ff)
}

\arguments{
  \item{ff}{
a \link[flowCore]{flowFrame} object
}

}

\value{
a \link[flowCore]{flowFrame} object with updated parameters
}

\author{
Yann Abraham
}

\examples{
require(flowCore)

ff <- read.FCS(system.file("extdata","Staurosporine_cd4+_H05.fcs",package="cytoCore"))

channels <- list(CyTOF=c('Time','Cell_length'),
     DNA=c('DNA-1','DNA-2'),
     phenotypic=c('CD3','CD45','CD4','CD20','CD33','CD123','CD14','HLA-DR','IgM','CD7'),
     functional=c('pNFkB','pp38','pStat5','pAkt','pStat1','pSHP2','pZap70','pStat3','pSlp76','pBtk','pPlcg2','pErk','pLat','pS6'),
     barcode=c('BC1','BC2','BC3','BC4','BC5','BC6','BC7')
)

cf <- make.cytoFrame(ff,channels)
cf <- correct.range(cf)
}

\keyword{ transform }
\keyword{ filtering }
