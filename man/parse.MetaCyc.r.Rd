\name{parse.MetaCyc.r}
\alias{parse.MetaCyc.r}
\title{
Parse MetaCyc reaction data (reactions.dat in attribute-value format)
}
\description{
Parse MetaCyc reactions.dat file in attribute-value format (BioCyc or EcoCyc file also applicable). The output includes reaction ID, enzyme, EC number, pathway, equation
}
\usage{
parse.MetaCyc.r(file)
}
\arguments{
  \item{file}{
Path to reactions.dat from MetaCyc (BioCyc or EcoCyc data also applicable)
}
}
\value{
Data frame with reaction properties
}
\references{
http://metacyc.org/
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\examples{
# sample reactions.dat 
## getting full data needs license; http://metacyc.org/download-flatfiles.shtml
url = 'http://bioinformatics.ai.sri.com/ptools/flatfile-samples/reactions.dat'
tmpdest = tempfile(pattern = "reactions") 
download.file(url, destfile=tmpdest)

direction = "LEFT-TO-RIGHT" # sample data is old, so don't mind this code.
sample.reactions = parse.MetaCyc.r(file=tmpdest)
print(sample.reactions)
}
