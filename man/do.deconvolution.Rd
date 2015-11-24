\name{do.deconvolution}
\alias{do.deconvolution}

\title{
Automatic deconvolution of barcodes in cytoFrame
}

\description{
Provided a data.frame or a \link[flowCore]{flowFrame} where some columns can be identified as barcode channels, the function will perform an automatic deconvolution using density distribution of each individual channels to identify inflection points and split the barcode into positive and negative. 
}

\usage{
do.deconvolution(x, barcodes = NULL, live = NULL, dest = ".", diagnostics = T, add = T)
}

\arguments{
  \item{x}{
either a data.frame or a \link[flowCore]{flowFrame} object containing barcode channels
}
  \item{barcodes}{
a list of barcode channel names; if NULL, the function will try to guess them from the channel names (see Details) 
}
  \item{live}{
either NULL or a logical vector specifying which cell should be considered live
}
  \item{dest}{
if diagnostics is TRUE, use dest to specify where the diagnostic plots will be saved
}
  \item{diagnostics}{
if diagnostics is TRUE, the function will return plots of the barcode distributions and deconvolution thresholds  
}
  \item{add}{
if add is TRUE, the barcodes will be added to the data rather than returned as a separate file (see value)
}
}

\details{
If barcodes is NULL and x is a \link[flowCore]{flowFrame}, the function will look for a type column: if type is present and contains a barcode category the corresponding channels will be used. Otherwise the function will look for channel names containig the string 'BC'
}

\value{
A data.frame object containing either the individual barcode channels plus an automated identifier; if add is TRUE and x in NOT a \link[flowCore]{flowFrame}, a combination of input data plus barcodes and identifiers
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

# deconvolute the cytoFrame
bc.cf <- do.deconvolution(cf,diagnostics=FALSE,add=FALSE)
table(droplevels(bc.cf$barcode)) # some noise in the deconvolution

# split the flowFrame using the barcode
cSet <- split(cf,bc.cf$barcode,flowSet=TRUE)
}

\keyword{ barcode }
\keyword{ deconvolution }
