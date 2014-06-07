\name{get.kegg.byId}
\alias{get.kegg.byId}
\title{
Download given KEGG REACTION entry using KEGGREST package
}
\description{
By implementing KEGGREST package in Bioconductor (http://www.bioconductor.org/packages/release/bioc/html/KEGGREST.html), retrieve given KEGG REACTION in data.frame() object.
}
\usage{
get.kegg.byId(keggId)
}
\arguments{
  \item{keggId}{ single KEGG REACTION ID (R number) or multiple IDs in vector}
}
\value{
data.frame() object where data fields include,
  \item{ENTRY }{KEGG ID}
  \item{NAME }{Enzyme name}
  \item{DEFINITION }{Reaction equation with standard compound name}
  \item{EQUATION }{Reaction equation with KEGG ID}
  \item{ENZYME }{E.C number}
  \item{COMMENT }{Comment}
  \item{RPAIR }{KEGG RPAIR is a collection of substrate-product pairs (reactant pairs) defined for each reaction in KEGG REACTION}
  \item{PATHWAY }{KEGG PATHWAY that this reaction participates}
  \item{ORTHOLOGY }{Cross-link to KEGG ORTHOLOGY database}
  \item{REMARK }{Usually this column indicates glycan reaction-chemical reaction pair}
  \item{REFERENCE}{Publication with PMID, title, and journal name}
}
\references{
KEGGREST package: http://www.bioconductor.org/packages/release/bioc/html/KEGGREST.html \cr
KEGG API: http://www.kegg.jp/kegg/docs/keggapi.html
}
\author{
Byoungnam Min <mbnmbn00@gmail.com>, Byeonghyuk Park, Kyoung Heon Kim and In-Geol Choi 
}
\note{
Because "KEGGREST" package has not been deposited in CRAN, you need to manually install this package by typing: \cr
\cr
> source("http://bioconductor.org/biocLite.R") \cr
> biocLite("KEGGREST") \cr
\cr
For details, \cr
http://www.bioconductor.org/packages/release/bioc/html/KEGGREST.html
}
\seealso{
\code{\link{get.kegg.all}}
}
\examples{
# Not run:

##  get.kegg.byId("R00001")                                  # single KEGG ID
##  get.kegg.byId(c("R00001", "R00002", "R00005", "R00006")) # mutiple KEGG IDs
}
\keyword{ KEGG REACTION }
\keyword{ KEGG REST }
\keyword{ KEGG API }
