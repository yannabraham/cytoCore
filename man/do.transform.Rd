\name{do.transform}
\alias{do.transform}

\title{
Applies transformation to selected columns
}

\description{
Applies a specific transformation to all columns specified either by name or type
}
\usage{
do.transform(ff, cols=NULL, type=NULL, fun=arcsinhTransform())
}

\arguments{
  \item{ff}{
a \link[flowCore]{flowFrame} object
}
  \item{cols}{

}
  \item{type}{
a column type to transform
}
  \item{fun}{
the transformation function to apply (defaults to \link[flowCore]{arcsinhTransform})
}
}

\details{
Type can be a single character string or a vector of column types
}
\value{
a \link[flowCore]{flowFrame} object with specified transformation applied to specified columns or types 
}

\author{
Yann Abraham
}

\examples{
ff <- read.FCS(system.file("extdata","Staurosporine_cd4+_H05.fcs",package="cytoCore"))

channels <- list(CyTOF=c('Time','Cell_length'),
     DNA=c('DNA-1','DNA-2'),
     phenotypic=c('CD3','CD45','CD4','CD20','CD33','CD123','CD14','HLA-DR','IgM','CD7'),
     functional=c('pNFkB','pp38','pStat5','pAkt','pStat1','pSHP2','pZap70','pStat3','pSlp76','pBtk','pPlcg2','pErk','pLat','pS6'),
     barcode=c('BC1','BC2','BC3','BC4','BC5','BC6','BC7')
)

cf <- make.cytoFrame(ff,channels)

# apply default transformation by column name
cf <- do.transform(cf,cols=c('DNA.1','DNA.2'))
# apply default transformation by type
cf <- do.transform(cf,type=c('phenotypic','functional'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
