\name{get.kegg.all}
\alias{get.kegg.all}
\title{
Download current KEGG REACTION and KEGG COMPOUND database using KEGGREST package
}
\description{
By implementing KEGGREST package in Bioconductor (http://www.bioconductor.org/packages/release/bioc/html/KEGGREST.html), gather all KEGG REACTION and KEGG COMPOUND datain data.frame() object. 
}
\usage{
get.kegg.all()
}
\value{
R list object. The list contains two data.frame() objects (reaction and comopund) where data fields include,
  \item{ENTRY }{KEGG ID (R number)}
  \item{NAME }{Enzyme name}
  \item{DEFINITION }{Reaction equation with standard compound name}
  \item{EQUATION }{Reaction equation with KEGG ID}
  \item{ENZYME }{E.C number}
  \item{COMMENT }{Comment}
  \item{RPAIR }{KEGG RPAIR is a collection of substrate-product pairs (reactant pairs) defined for each reaction in KEGG REACTION}
  \item{PATHWAY }{KEGG PATHWAY that this reaction participates}
  \item{ORTHOLOGY }{Cross-link to KEGG ORTHOLOGY database}
  \item{REMARK }{Usually this column indicates glycan reaction-chemical reaction pair}
  \item{REFERENCE}{Publication info with PMID, title, and journal name}
\cr
KEGG COMPOUND,
  \item{ENTRY }{KEGG ID (C number)}
  \item{NAME }{Compound name}
  \item{FORMULA }{Molecular formula}
  \item{EXACT_MASS }{Molecular mass}
  \item{MOL_WEIGHT }{Molecular weight}
  \item{REMARK }{}
  \item{REACTION }{KEGG REACTION entries where the compound particiaptes}
  \item{PATHWAY }{KEGG PATHWAY entries where the compound participates}
  \item{ENZYME }{KEGG ENZYME entries where the compound is related}
  \item{ATOM }{}
  \item{BOND }{}
  \item{CAS }{Cross-link to CAS database}
  \item{PubChem }{Cross-link to PubChem database}
  \item{ChEBI }{Cross-link to ChEBI database}
  \item{PDB.CCD }{Cross-link to Chemical Component Dictionary}
  \item{X3DMET }{}
  \item{NIKKAJI }{}
  \item{BRITE }{}
  \item{KNApSAcK }{}
  \item{COMMENT }{}
  \item{BRACKET }{}
  \item{LIPIDMAPS }{}
  \item{LipidBank }{}
  \item{SEQUENCE }{}
  \item{REFERENCE }{}
  
}
\references{
KEGGREST package: http://www.bioconductor.org/packages/release/bioc/html/KEGGREST.html \cr
KEGG API: http://www.kegg.jp/kegg/docs/keggapi.html
}
\author{
Byoungnam Min <mbnmbn00@gmail.com>, Byeonghyuk Park, Kyoung Heon Kim and In-Geol Choi
}
\note{
Because 'KEGGREST' package has not been deposited in CRAN, you need to manually install this package by typing: \cr
\cr
> source("http://bioconductor.org/biocLite.R") \cr
> biocLite("KEGGREST") \cr
\cr
For details, \cr
http://www.bioconductor.org/packages/release/bioc/html/KEGGREST.html \cr
}

\seealso{
\code{\link{get.kegg.byId}}
}

\examples{
# Not run: (it will take 30-40 minutes)

##  keggAll = get.kegg.all()
##  head(keggAll$reaction)
##  head(keggAll$comound)

}
\keyword{ KEGG REACTION } 
\keyword{ KEGG COMPOUND }
\keyword{ KEGG REST }
\keyword{ KEGG API }
