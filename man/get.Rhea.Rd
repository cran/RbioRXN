\name{get.Rhea}
\alias{get.Rhea}
\title{
Get Rhea data in BioPAX level 2 format
}
\description{
Get Rhea data in BioPAX level 2 format (ftp://ftp.ebi.ac.uk/pub/databases/rhea/biopax/rhea-biopax.owl.gz)
}
\usage{
get.Rhea()
}
\value{
Character vector represented in RDF/XML syntax
}
\references{
http://www.ebi.ac.uk/rhea/download.xhtml;jsessionid=1A621DF38331F644AD1730C328457654
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\seealso{
\code{\link{parse.Rhea}}
}
\examples{
# Not run:

## Rhea = get.Rhea()

## length(Rhea)
## head(Rhea)
}
