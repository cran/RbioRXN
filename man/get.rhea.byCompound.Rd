\name{get.rhea.byCompound}
\alias{get.rhea.byCompound}
\title{
Download Rhea reaction by compound
}
\description{
Download and parse Rhea database with given keyword
}
\usage{
get.rhea.byCompound(rheaCompound)
}
\arguments{
  \item{rheaCompound}{ e.g., caffeine }
}

\value{
data.frame() object where data fields include
  \item{rheaId }{Rhea ID}
  \item{reactionType }{BiochemicalReaction or transportReaction}
  \item{equationWithCommonName }{Textual reaction equation (compound standard names)}
  \item{sameParticipant }{Reactions with same reactants and products, but different directionalities}
  \item{mapped }{}
  \item{formuled }{}
  \item{polymerization }{}
  \item{chemicallyBalanced }{}
  \item{status }{}
  \item{transport }{}
  \item{direction }{}
  \item{classOfReactions }{At least one reaction participant is class compound described in ChEBI}
  \item{equationWithChebi }{Reaction equation with ChEBI}
  \item{equationParticipant }{}
  \item{crossLinks }{if reaction has cross-link to other databases, new columns are generated such as "MetaCyc", "KEGG" etc}
}
\references{
Rhea database: http://www.ebi.ac.uk/rhea/home.xhtml
}
\author{
Byoungnam Min <mbnmbn00@gmail.com>, Byeonghyuk Park, Kyoung Heon Kim and In-Geol Choi 
}
\examples{
# Not run:

##  get.rhea.byCompound("caffeine")
}
\keyword{ Rhea }
