\name{Rhea.instantiate}
\alias{Rhea.instantiate}
\title{
Instantiate Rhea generic reaction
}
\description{
Generic reaction is a reaction having at least one compound class. This function gains instance reactions from generic reaction
}
\usage{
Rhea.instantiate(parsed_Rhea, parsed_ChEBI, Rhea_ID, multicore = 1)
}
\arguments{
  \item{parsed_Rhea}{
parse.Rhea function should be preceded and provided as input
}
  \item{parsed_ChEBI}{
parse.ChEBI function should be preceded and provided as input
}
  \item{Rhea_ID}{
Rhea reaction ID to be instantiated (generic reaction ID)
}
  \item{multicore}{
Using multicore for parallel computing (doMC package)
}
}
\value{
\item{ID}{Instance reaction ID (e.g., 10588_1, 10588_2)}
\item{Instance reaction}{Generated instance reactions}
}
\references{
http://www.ebi.ac.uk/rhea/
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\seealso{
\code{\link{get.Rhea}}\cr
\code{\link{parse.Rhea}}\cr
\code{\link{get.ChEBI}}\cr
\code{\link{parse.ChEBI}}\cr
\code{\link{Rhea.is.generic}}
}
\examples{
data(example) # load example data

## Rhea instantiation (Rhea.instantiate(parsed_Rhea,parsed_ChEBI,Rhea_ID,multicore=1))
pC = example$parsed_ChEBI # sample ChEBI
Rg = example$Rhea_generic # sample Rhea generic reaction

data(thermo)
instanceR = Rhea.instantiate(Rg, pC, Rg[1,'ID'], multicore=1)
print(instanceR)
}
