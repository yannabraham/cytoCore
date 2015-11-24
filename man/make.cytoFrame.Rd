\name{make.cytoFrame}
\alias{make.cytoFrame}

\title{
Turns a flowFrame into a cytoFrame
}

\description{
After loading CyTOF data as FCS, updates the parameters slot of the flowFrame so that isotopes & channel name are clearly separated; adds a type column to explicitly set parameter type 
}

\usage{
make.cytoFrame(ff, channels)
}

\arguments{
  \item{ff}{
a \link[flowCore]{flowFrame} object
}
  \item{channels}{
a list of channels sorted by type (see details)
}
}

\details{
Channels 
}

\value{
a \link[flowCore]{flowFrame} object with extra type and isotope columns
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
}

\keyword{ CyTOF }
\keyword{ mass cytometry }
