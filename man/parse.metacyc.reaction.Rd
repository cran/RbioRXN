\name{parse.metacyc.reaction}
\alias{parse.metacyc.reaction}
\title{
Parse MetaCyc reaction data (reactions.dat in attribute-value format)
}
\description{
Parse MetaCyc reactions.dat file in attribute-value format. The output includes reaction ID, enzyme, EC number, pathway and equation
}
\usage{
parse.metacyc.reaction(datPath)
}
\arguments{
  \item{datPath}{
Path to reactions.dat from MetaCyc (BioCyc or EcoCyc data also applicable)
}
}
\value{
data.frame() object where data fields include
  \item{UNIQUE.ID }{}
  \item{TYPES }{}
  \item{ATOM.MAPPINGS }{}
  \item{CITATIONS }{}
  \item{COMMENT }{}
  \item{CREDITS }{}
  \item{DBLINKS }{}
  \item{IN.PATHWAY }{}
  \item{LEFT }{}
  \item{PHYSIOLOGICALLY.RELEVANT. }{}
  \item{REACTION.DIRECTION }{}
  \item{RIGHT }{}
  \item{CANNOT.BALANCE. }{}
  \item{ENZYMATIC.REACTION }{}
  \item{ORPHAN. }{}
  \item{SPONTANEOUS. }{}
  \item{RXN.LOCATIONS }{}
  \item{EC.NUMBER }{}
  \item{COMMON.NAME }{}
  \item{TEMPLATE.FILE }{}
  \item{SYNONYMS }{}
  \item{SYSTEMATIC.NAME }{}
  \item{DELTAG0 }{}
  \item{COMMENT.INTERNAL }{}
  \item{ENZYMES.NOT.USED }{}
  \item{REACTION.LIST }{}
  \item{MEMBER.SORT.FN }{}
  \item{STD.REDUCTION.POTENTIAL }{}
  \item{PREDECESSORS }{}
  \item{PRIMARIES }{}
  \item{SIGNAL }{}
  \item{equation }{}
  \item{LIGAND.RXN }{}
  \item{PLANTCYC }{}
  \item{METACYC }{}
  \item{RHEA }{}
  \item{UNIPROT }{}
  \item{PIR }{}
  \item{LIGAND }{}
  \item{BRENDA }{}
  \item{superNode }{}
}
\author{
Byoungnam Min <mbnmbn00@gmail.com>, Byeonghyuk Park, Kyoung Heon Kim and In-Geol Choi 
}
\examples{
url = "http://bioinformatics.ai.sri.com/ptools/flatfile-samples/reactions.dat"
tmpdest = tempfile(pattern = "reactions")
download.file(url, destfile=tmpdest)

direction = "LEFT-TO-RIGHT" # sample data is old, so don't mind this code.
sample.reactions = parse.metacyc.reaction(datPath=tmpdest)
print(sample.reactions)
}
\keyword{ MetaCyc }
