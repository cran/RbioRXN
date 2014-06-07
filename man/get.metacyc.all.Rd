\name{get.metacyc.all}
\alias{get.metacyc.all}
\title{
Download current MetaCyc via BioCyc Web Services
}
\description{
This function first retrieves all MetaCyc reaction IDs and download BioPAX level 3 for each reaction, then parse to generate data.frame() object. Note that because this function interacts with BioCyc web service for every reaction, it will take many hours depending on your network speed. Thus, we strongly recommend to use parse.metacyc.reaction() and parse.metacyc.compound() functions in this package after manually downloading MetaCyc from http://biocyc.org/download.shtml
}
\usage{
get.metacyc.all()
}
\value{
R data.frame() object where data fields indclude
  \item{MetaCyc }{MetaCyc ID}
  \item{name }{Reaction equation with chemical standard name}
  \item{ec_number }{E.C number}
  \item{name_id }{Reaction equation with MetaCyc compound ID}
}
\references{
MetaCyc download page: http://biocyc.org/download.shtml \cr
MetaCyc web services: http://biocyc.org/web-services.shtml
}
\author{
Byoungnam Min <mbnmbn00@gmail.com>, Byeonghyuk Park, Kyoung Heon Kim and In-Geol Choi
}
\note{
We strongly recommend using parse.metacyc.reaction() and parse.metacyc.reaction() functions because of long running process. These two functions need "dat" files downloaded from MetaCyc. See http://biocyc.org/download.shtml
}
\seealso{
\code{\link{parse.metacyc.reaction}} 
\code{\link{parse.metacyc.compound}}
}
\examples{
# Not run:

##  metacycAll = get.metacyc.all()
##  head(metacycAll)
}
\keyword{ Metacyc web services } 
\keyword{ BioPAX level 3 }
