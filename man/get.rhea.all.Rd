\name{get.rhea.all}
\alias{get.rhea.all}
\title{
Download current Rhea in data.frame() object
}
\description{
Download current Rhea in BioPAX level 2 and parse to generate data.frame() object
}
\usage{
get.rhea.all()
}
\value{
data.frame() object where data fields include
  \item{rheaId }{Rhea ID}
  \item{reactionType }{BiochemicalReaction or transportReaction}
  \item{equationWithCommonName }{Textual reaction equation (compound standard names)}
  \item{sameParticipant }{Reactions with the same reactants and products, but different directionalities}
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
  \item{crossLinks }{If reaction has cross-link to other databases, new columns are generated such as "MetaCyc", "KEGG" etc}
}
\references{
Rhea database: http://www.ebi.ac.uk/rhea/home.xhtml
}
\author{
Byoungnam Min <mbnmbn00@gmail.com>, Byeonghyuk Park, Kyoung Heon Kim and In-Geol Choi 
}
\examples{
# Not run: it takes about 20 minutes

##  rhea = get.rhea.all()
##  dim(rhea)              # The number of downloaded reactions
##  head(rhea)

}
\keyword{ Rhea } 
