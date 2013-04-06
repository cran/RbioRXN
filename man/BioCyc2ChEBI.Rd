\name{BioCyc2ChEBI}
\alias{BioCyc2ChEBI}
\title{
Convert BioCyc ID into ChEBI ID in reaction equation
}
\description{
BioCyc reaction equation is composed of reaction participants expressed with BioCyc ID. This function converts BioCyc ID into its counterpart ChEBI ID. Note that some reactions connot be converted because of imcomplete crosslink between the databases.
}
\usage{
BioCyc2ChEBI(parsed_MetaCyc.c, eq)
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
Reaction equation expressed with ChEBI ID.
}
\references{
http://www.metacyc.org/
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\seealso{
\code{\link{BioCyc2cName}}\cr
\code{\link{BioCyc2PubChem}}\cr
\code{\link{BioCyc2KEGG}}\cr
\code{\link{parse.MetaCyc.c}}
}
\examples{
data(example) # load sample MetaCyc data

pMc = example$parsed_MetaCyc.c # sample parsed MetaCyc compound data

Mc = example$MetaCyc_conv # sample MetaCyc reaction
print(Mc)

B2ChEBI = BioCyc2ChEBI(pMc, Mc) # BioCyc ID to ChEBI ID
print(B2ChEBI)
}
