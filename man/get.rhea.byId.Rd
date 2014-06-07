\name{get.rhea.byId}
\alias{get.rhea.byId}
\title{
Download Rhea by Rhea ID
}
\description{
Download and parse Rhea with given Rhea ID
}
\usage{
get.rhea.byId(rheaId)
}
\arguments{
  \item{rheaId}{ Single Rhea ID or multiple Rhea IDs in vector }
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
  \item{classOfReactions }{At least one participant is class compounds described in ChEBI}
  \item{equationWithChebi }{Reaction equation with ChEBI}
  \item{equationParticipant }{}
  \item{crossLinks }{If the reaction has cross-link to other databases, new columns are generated such as "MetaCyc", "KEGG" etc}
}
\references{
Rhea database: http://www.ebi.ac.uk/rhea/home.xhtml
}
\author{
Byoungnam Min <mbnmbn00@gmail.com>, Byeonghyuk Park, Kyoung Heon Kim and In-Geol Choi 
}
\examples{
# Not run:

##  get.rhea.byId(c('10000', '10377'))
}
\keyword{ Rhea }
