\name{parse.Rhea}
\alias{parse.Rhea}
\title{
Parse Rhea data (rhea.owl in BioPAX level 2 format)
}
\description{
Parse reaction ID and equation from Rhea data in BioPAX level2 format
}
\usage{
parse.Rhea(owl)
}
\arguments{
  \item{owl}{
get.Rhea function should be precded and provided as input
}
}
\value{
data frame with reaction properties
}
\references{
http://www.ebi.ac.uk/rhea/
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\seealso{
\code{\link{get.Rhea}}
}
\examples{
# Not run:

## Rhea = get.Rhea()
## parsed_Rhea = parse.Rhea(Rhea)

## head(parsed_Rhea)
## dim(parsed_Rhea) # get dimension of data frame
}
