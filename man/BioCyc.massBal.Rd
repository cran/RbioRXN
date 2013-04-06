\name{BioCyc.massBal}
\alias{BioCyc.massBal}
\title{
Check mass balance
}
\description{
Check if input reactions are mass-balanced. This function counts atoms in a reaction equation based on molecular formula. Note that some reaction participants don't have molecular formula (e.g., polymer or generic compound), so these reactions cannot be checked.
}
\usage{
BioCyc.massBal(parsed_MetaCyc.c, equation)
}
\arguments{
  \item{parsed_MetaCyc.c}{
parse.MetaCyc.c function should be preceded and provided as input
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
http://biocyc.org
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\seealso{
\code{\link{parse.MetaCyc.c}}
}
\examples{
data(example) # load example data

pMc = example$parsed_MetaCyc.c # sample MetaCyc compound
Mm = example$MetaCyc_massBal # sample MetaCyc reaction

data(thermo) # this is for the package 'CHNOSZ'
BioCyc.massBal(pMc, Mm)
}
