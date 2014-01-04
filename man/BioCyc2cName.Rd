\name{BioCyc2cName}
\alias{BioCyc2cName}
\title{
Convert BioCyc ID into ChEBI ID in reaction equation
}
\description{
BioCyc reaction equation is composed of reaction participants expressed with BioCyc ID. This function converts BioCyc ID into compound name
}
\usage{
BioCyc2cName(parsed_MetaCyc.c, eq)
}
\arguments{
  \item{parsed_MetaCyc.c}{
parse.MetaCyc.c function should be preceded and provided as input
}
  \item{eq}{
Reaction equation to be converted
}
}
\value{
Reaction equation expressed with compound name
}
\references{
http://www.metacyc.org/
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\seealso{
\code{\link{BioCyc2PubChem}}\cr
\code{\link{BioCyc2ChEBI}}\cr
\code{\link{BioCyc2KEGG}}\cr
\code{\link{parse.MetaCyc.c}}
}
\examples{
data(example) # load sample MetaCyc data

pMc = example$parsed_MetaCyc.c # sample parsed MetaCyc compound

Mc = example$MetaCyc_conv # sample MetaCyc compound
print(Mc)

B2cName = BioCyc2cName(pMc, Mc) # BioCyc ID to compound name
print(B2cName)
}
