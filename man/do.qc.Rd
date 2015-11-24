\name{do.qc}
\alias{do.qc}

\title{
Performs standard Quality Control on a flowFrame object
}

\description{
The function generates a set of quality control plots, checking basic features of the flowFrame object  
}
\usage{
do.qc(ff, cols=NULL, type=NULL, out.dir='.')
}

\arguments{
  \item{ff}{
a \link[flowCore]{flowFrame} object
}
  \item{cols}{
a list of columns to transform
}
  \item{type}{
a column type to transform
}
  \item{out.dir}{
the directory where the image will be saved
}
}

\details{
Cols and type can be a single character string or a vector of column types.
CyTOF-specific columns (Leading_push, Time, Cell_length) are included by default. 
}
\value{
Images corresponding to standard quality control plots are saved to out.dir. Function returns an invisible NULL.
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

# turn the flowFrame into a cytoFrame
cf <- make.cytoFrame(ff,channels)

# apply default transformation by type
cf <- do.transform(cf,type=c('DNA','phenotypic','functional'))
}

# perform QC
do.qc(cf,type=c('DNA','phenotypic','functional'),out.dir='.')
