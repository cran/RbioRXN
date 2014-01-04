\name{Rhea.massBal}
\alias{Rhea.massBal}
\title{
Check mass balance
}
\description{
Check if input reactions are mass-balanced. This function counts atoms in the left and right side of reaction equation based on molecular formula. Note that some reaction participants don't have molecular formula (e.g., polymer or generic compound), so these reactions cannot be checked.
}
\usage{
Rhea.massBal(parsed_ChEBI, equation)
}
\arguments{
  \item{parsed_ChEBI}{
parse.ChEBI function should be preceded and provided as input
}
  \item{equation}{
Reaction equation
}
}
\value{
Character TRUE/FALSE/cannot be checked
}
\references{
CHNOSZ package: http://cran.r-project.org/web/packages/CHNOSZ/index.html \cr
http://www.ebi.ac.uk/rhea/
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\seealso{
\code{\link{parse.ChEBI}}
}
\examples{
data(example) # load example data

pC = example$parsed_ChEBI # sample ChEBI

## Rhea check mass balance (Rhea.massBal(parsed_ChEBI,equation))
Rm = example$Rhea_massBal # sample Rhea data

data(thermo) # this is for package(CHNOSZ)
checkR = Rhea.massBal(pC, Rm)
print(checkR)
}
