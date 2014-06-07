\name{get.chebi.all}
\alias{get.chebi.all}
\title{
Download and parse current ChEBI
}
\description{
Download all ChEBI in BioPAX level 2 (ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl) and parse to generate data.frame() object
}
\usage{
get.chebi.all()
}
\value{
data.frame() object where data fields include
    \item{ChEBI }{Unique ChEBI identifier}
    \item{name }{Chemical standard name}
    \item{synonyms }{If entity has multiple synonyms, they are separated by "///"}
    \item{smiles }{SMILES representation of chemical structure}
    \item{inchi }{InChI representation of chemical structure}
    \item{kegg }{Cross-link to KEGG COMPOUND database}
    \item{parent }{Parent compound if available}
    \item{formula }{Chemical formula}
}
\references{
ChEBI web address: http://www.ebi.ac.uk/chebi/ \\
BioPAX.org: http://www.biopax.org
}
\author{
Byoungnam Min, Byeonghyuk Park, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\examples{
# Not run: it will take several minutes (5-7 mins)

##  chebi = get.chebi.all()
##  head(chebi)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ChEBI }
\keyword{ BioPAX }% __ONLY ONE__ keyword per line
