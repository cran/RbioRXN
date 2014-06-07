\name{get.metacyc.reaction.byId}
\alias{get.metacyc.reaction.byId}
\title{
Download given MetaCyc reaction via BioCyc web services
}
\description{
This function downloads and parses given MetaCyc reaction in BioPAX level 3 via MetaCyc web services (http://biocyc.org/web-services.shtml), and generate data.frame() object
}
\usage{
get.metacyc.reaction.byId(metacycId)
}
\arguments{
  \item{metacycId }{MetaCyc ID: single or multiple entries in vector}
}
\value{
data.frame() object where data fields include
  \item{MetaCyc }{MetaCyc ID}
  \item{name }{Reaction equation with chemical standard name}
  \item{ec_number }{E.C. number}
  \item{name_id }{Reaction equation with MetaCyc compound ID}
}
\references{
MetaCyc download page: http://biocyc.org/download.shtml \cr
MetaCyc web services: http://biocyc.org/web-services.shtml
}
\author{
Byoungnam Min <mbnmbn00@gmail.com>, Byeonghyuk Park, Kyoung Heon Kim and In-Geol Choi 
}
\examples{
# Not run:

##  get.metacyc.reaction.byId("RXN-13970")
##  get.metacyc.reaction.byId("SERINE--PYRUVATE-AMINOTRANSFERASE-RXN")
##  get.metacyc.reaction.byId("RXN-13759")
##  get.metacyc.reaction.byId("RXN-14352")
##  get.metacyc.reaction.byId("TRANS-RXN-220")
##  get.metacyc.reaction.byId("RXN-12549")
}
\keyword{ Metacyc web services } 
\keyword{ BioPAX level 3 }
