\name{cytoCore-package}
\alias{cytoCore-package}
\alias{cytoCore}
\docType{package}
\title{
Simple helper functions to deal with CyTOF data
}
\description{
Simplifies the usage of flowCore with CyTOF data by providing simple function to load, annotate and deconvolute barcodes for CyTOF experiments
}
\details{
\tabular{ll}{
Package: \tab cytoCore\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2013-11-13\cr
License: \tab CC BY-NC-SA\cr
}

}
\author{
Yann Abraham
Maintainer: Yann Abraham <yann.abraham@novartis.com>
}

\keyword{ package }
\keyword{ CyTOF }
\keyword{ mass cytometry }

\seealso{
\code{\link[flowCore:flowCore-package]{flowCore}}
}

\examples{
# load example data
ff <- read.FCS(system.file("extdata","Staurosporine_cd4+_H05.fcs",package="cytoCore"))

# specify channels
channels <- list(CyTOF=c('Time','Cell_length'),
     DNA=c('DNA-1','DNA-2'),
     phenotypic=c('CD3','CD45','CD4','CD20','CD33','CD123','CD14','HLA-DR','IgM','CD7'),
     functional=c('pNFkB','pp38','pStat5','pAkt','pStat1','pSHP2','pZap70','pStat3','pSlp76','pBtk','pPlcg2','pErk','pLat','pS6'),
     barcode=c('BC1','BC2','BC3','BC4','BC5','BC6','BC7')
)

# annotate flowFrame
cf <- make.cytoFrame(ff,channels)

# normalize channels
cf <- do.transform(cf,type=c('phenotypic','functional'))

# deconvolute barcodes
bc.cf <- do.deconvolution(cf,diagnostics=FALSE,add=FALSE)
table(droplevels(bc.cf$barcode)) # some noise in the deconvolution

# split the flowFrame using the barcode
cSet <- split(cf,bc.cf$barcode,flowSet=TRUE)

# correct the range for individual files
cSet <- fsApply(cSet,correct.range)
}