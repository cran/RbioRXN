\name{Rhea2cName}
\alias{Rhea2cName}
\title{
Convert ChEBI ID into compound name in reaction equation
}
\description{
BioCyc reaction equation is composed of reaction participants expressed with ChEBI ID. This function converts ChEBI ID into counterpart compound name
}
\usage{
Rhea2cName(parsed_ChEBI, eq)
}
\arguments{
  \item{parsed_ChEBI}{
parse.ChEBI function should be preceded and provided as input
}
  \item{eq}{
Reaction equation to be converted
}
}
\value{
Reaction equation expressed with compound name
}
\references{
http://www.ebi.ac.uk/
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\seealso{
\code{\link{Rhea2KEGG}}\cr
\code{\link{get.ChEBI}}\cr
\code{\link{parse.ChEBI}}
}
\examples{
data(example) # load example data

pC = example$parsed_ChEBI # sample ChEBI
Rc = example$Rhea_conv # sample Rhea data
print(Rc)

R2cName = Rhea2cName(pC, Rc) # ChEBI ID to compound name
print(R2cName)
}
