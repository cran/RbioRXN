\name{Rhea.is.generic}
\alias{Rhea.is.generic}
\title{
Check generic reaction
}
\description{
Generic reaction is a reaction having at least one compound class as reaction participant.
}
\usage{
Rhea.is.generic(parsed_ChEBI, equation)
}
\arguments{
  \item{parsed_ChEBI}{
parse_ChEBI function should be preceded and provided as input
}
  \item{equation}{
Reaction equation to be checked
}
}
\value{
logical TRUE/FALSE
}
\references{
http://www.ebi.ac.uk/rhea/
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\seealso{
\code{\link{BioCyc.is.generic}}
}
\examples{
data(example) # load example data

pC = example$parsed_ChEBI # sample ChEBI
Rg = example$Rhea_generic # sample Rhea generic reaction

Rhea.is.generic(pC, Rg[,'Eq.ChEBI'])
}
