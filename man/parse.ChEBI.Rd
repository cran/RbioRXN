\name{parse.ChEBI}
\alias{parse.ChEBI}
\title{
Parse ChEBI raw data (chebi.owl in BioPAX level 2 format))
}
\description{
Parse chemical properties from ChEBI data in BioPAX level 2 format such as compound name, synonyms, SMILES, InChi, KEGG ID, compound class
}
\usage{
parse.ChEBI(owl)
}
\arguments{
  \item{owl}{
get.ChEBI function should be precded and provided as input
}
}
\value{
Data frame with chemical properties
}
\references{
http://www.ebi.ac.uk/chebi/
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\seealso{
\code{\link{get.ChEBI}}
}
\examples{
# Not run:

## ChEBI = get.ChEBI()
## parsed_ChEBI = parse.ChEBI(ChEBI)

## head(parsed_ChEBI)
## dim(parsed_ChEBI) # get dimension of data frame
}
