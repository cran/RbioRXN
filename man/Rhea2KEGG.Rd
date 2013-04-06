\name{Rhea2KEGG}
\alias{Rhea2KEGG}
\title{
Convert ChEBI ID into KEGG ID in reaction equation
}
\description{
Rhea reaction equation is composed of reaction participants expressed with ChEBI ID. This function converts ChEBI ID into counterpart KEGG ID. Note that some reactions connot be converted because of imcomplete crosslink between the databases.
}
\usage{
Rhea2KEGG(parsed_ChEBI, eq)
}
\arguments{
  \item{parsed_ChEBI}{
parse_ChEBI function should be preceded and provided as input
}
  \item{eq}{
Reaction equation to be converted
}
}
\value{
Reaction equation expressed with KEGG ID.
}
\references{
http://www.ebi.ac.uk/rhea/ \cr
http://www.ebi.ac.uk/chebi/
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\seealso{
\code{\link{get.ChEBI}}\cr
\code{\link{parse.ChEBI}}\cr
\code{\link{Rhea2cName}}\cr
}
\examples{
data(example) # load example data

pC = example$parsed_ChEBI # sample ChEBI
Rc = example$Rhea_conv # sample Rhea
print(Rc)

R2KEGG = Rhea2KEGG(pC, Rc) # ChEBI ID to KEGG ID
print(R2KEGG)
}
