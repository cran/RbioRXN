\name{get.unipathway.all}
\alias{get.unipathway.all}
\title{
Download current Unipathway
}
\description{
Download Unipathway in OBO, then parse and generate data.frame() object
}
\usage{
get.unipathway.all(oboPath = "")
}
\arguments{
  \item{oboPath}{ Unipathway in OBO from http://www.grenoble.prabi.fr/obiwarehouse/unipathway. If you don't provide the file, the function tries to download from web
}
}
\value{
R list object. The list() contains two data.frame() objects (reaction and compound) where data fields of each data.frame() include, \cr
reaction includes
  \item{id }{Unipathway ID}
  \item{enzName }{It usually indicates substrate and product in the reaction}
  \item{equation }{Reaction equation with compound standard name}
  \item{KEGG }{Cross-link to KEGG}
  \item{RHEA }{Cross-link to Rhea}
  \item{METACYC }{Cross-link to MetaCyc}
  \item{EC }{EC number each reaction assigned }
  \item{GO }{GO term}
  \item{UNIPROT }{Cross-link to UniProt}
  \item{compoundId }{Compound ID of each reaction participant}
  \item{compoundName }{Compound standard name of each reaction participant}
\cr
compound includes
  \item{id }{Unipathway ID}
  \item{type }{Compound type}
  \item{name }{Compound name}
  \item{KEGG }{Cross-link to KEGG}
  \item{CHEBI }{Cross-link to ChEBI}
  \item{formula }{Molecular formula}
  \item{inchi }{InChI code}
  \item{synonym.kegg }{Synonyms to KEGG}
}
\references{
UniPathway: http://www.grenoble.prabi.fr/obiwarehouse/unipathway
}
\author{
Byoungnam Min <mbnmbn00@gmail.com>, Byeonghyuk Park, Kyoung Heon Kim and In-Geol Choi 
}
\examples{
# Not run:

##  unipathway = get.unipathway.all()
##  dim(unipathway)                   # The number of downloaded reaction entries
}
\keyword{ UniPathway }
