\name{get.ChEBI}
\alias{get.ChEBI}
\title{
Get ChEBI data in BioPAX level 2 format
}
\description{
Get ChEBI data in BioPAX level 2 format (ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl)
}
\usage{
get.ChEBI()
}
\value{
Character vector represented in RDF/XML syntax
}
\references{
http://www.ebi.ac.uk/chebi/downloadsForward.do
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\seealso{
\code{\link{parse.ChEBI}}
}
\examples{
# Not run:

## ChEBI = get.ChEBI()

## length(ChEBI)
## head(ChEBI)
}
