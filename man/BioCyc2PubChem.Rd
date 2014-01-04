\name{BioCyc2PubChem}
\alias{BioCyc2PubChem}
\title{
Express BioCyc reaction equation with PubChem ID
}
\description{
BioCyc reaction equation is composed of reaction participants expressed with BioCyc ID. This function converts BioCyc ID into counterpart PubChem ID. Note that some reactions connot be converted because of imcomplete crosslink between the databases.
}
\usage{
BioCyc2PubChem(parsed_MetaCyc.c, eq)
}
\arguments{
  \item{parsed_MetaCyc.c}{
parse_MetaCyc.c function should be preceded and provided as input
}
  \item{eq}{
Reaction equation to be converted
}
}
\value{
Reaction equation expressed with PubChem ID.
}
\references{
http://www.metacyc.org/
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\seealso{
\code{\link{BioCyc2cName}}\cr
\code{\link{BioCyc2ChEBI}}\cr
\code{\link{BioCyc2KEGG}}\cr
\code{\link{parse.MetaCyc.c}}
}
\examples{
data(example) # load sample MetaCyc data

pMc = example$parsed_MetaCyc.c # sample parsed MetaCyc compound

Mc = example$MetaCyc_conv # sample MetaCyc reaction
print(Mc)

B2PubChem = BioCyc2PubChem(pMc, Mc) # BioCyc ID to PubChem ID
print(B2PubChem)
}
