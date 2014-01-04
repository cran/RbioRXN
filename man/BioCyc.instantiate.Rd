\name{BioCyc.instantiate}
\alias{BioCyc.instantiate}
\title{
Instantiate BioCyc generic reaction (MetaCyc, BioCyc, EcoCyc)
}
\description{
Generic reaction is a reaction having at least one compound class. This function gains instance reactions from generic reaction
}
\usage{
BioCyc.instantiate(parsed_MetaCyc.r, parsed_MetaCyc.c, BioCyc_ID)
}
\arguments{
  \item{parsed_MetaCyc.r}{
parse.MetaCyc.r function should be preceded and provided as input
}
  \item{parsed_MetaCyc.c}{
parse.MetaCyc.c function should be preceded and provided as input
}
  \item{BioCyc_ID}{
BioCyc reaction ID to be instantiated (generic reaction ID)
}
}
\value{
\item{ID}{Instance reaction ID (e.g., RXN-11376_1, RXN-11376_2)}
\item{Instance reaction}{Generated instance reactions}
}
\author{
Byoungnam Min, Kyoung Heon Kim, In-Geol Choi <mbnmbn00@gmail.com>
}
\seealso{
\code{\link{parse.MetaCyc.c}}\cr
\code{\link{parse.MetaCyc.r}}\cr
\code{\link{BioCyc.is.generic}}
}
\examples{
data(example) # load example data

## MetaCyc instantiation
pMc = example$parsed_MetaCyc.c
Mg = example$MetaCyc_generic

data(thermo) # this is for the package 'CHNOSZ'
instanceM = BioCyc.instantiate(Mg, pMc, Mg[1,'ID'])
print(instanceM)
}
