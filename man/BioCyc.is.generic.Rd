\name{BioCyc.is.generic}
\alias{BioCyc.is.generic}
\title{
Check generic reaction
}
\description{
Generic reaction is a reaction having at least one compound class as reaction participant (e.g., a primary alcohol + NADP -> an aldehyde + H(+) + NADPH)
}
\usage{
BioCyc.is.generic(parsed_MetaCyc.c, equation)
}
\arguments{
  \item{parsed_MetaCyc.c}{
parse_MetaCyc.c function should be preceded and provided as input
}
  \item{equation}{
Reaction equation to be checked
}
}
\value{
logical TRUE/FALSE
}
\references{
http://www.metacyc.org/
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\seealso{
\code{\link{BioCyc.instantiate}}
}
\examples{
data(example) # load example data

## MetaCyc instantiation (BioCyc.instantiate(parsed_MetaCyc.r,parsed_MetaCyc.c,BioCyc_ID,multicore))
pMc = example$parsed_MetaCyc.c # sample MetaCyc compound
Mg = example$MetaCyc_generic # sample generic reaction

BioCyc.is.generic(pMc, Mg[,'Equation'])
}
