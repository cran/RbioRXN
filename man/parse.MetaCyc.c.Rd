\name{parse.MetaCyc.c}
\alias{parse.MetaCyc.c}
\title{
parse MetaCyc compound data (compounds.dat in attribute-value format)
}
\description{
Parse MetaCyc compounds.dat file in attribute-value format (BioCyc and EcoCyc data also applicable). The output includes compound name, synonyms, ChEBI ID, PubChem ID, KEGG ID, CAS ID, types, molecular formula, SMILES code and InChI
}
\usage{
parse.MetaCyc.c(file)
}
\arguments{
  \item{file}{
Path to compounds.dat file from MetaCyc (BioCyc and EcoCyc data also applicable)
}
}
\value{
Data frame with chemical properties
}
\references{
http://metacyc.org/
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi <mbnmbn00@gmail.com>
}
\examples{
# sample compounds.dat, full data needs license; http://metacyc.org/download-flatfiles.shtml
url = 'http://bioinformatics.ai.sri.com/ptools/flatfile-samples/compounds.dat'
tmpdest = tempfile(pattern = "compounds")
download.file(url, destfile=tmpdest)

sample.compounds = parse.MetaCyc.c(file=tmpdest)
print(sample.compounds)
}
