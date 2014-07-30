\name{parse.metacyc.compound}
\alias{parse.metacyc.compound}
\title{
Parse MetaCyc compound data (compounds.dat in attribute-value format)
}
\description{
Parse MetaCyc compounds.dat file in attribute-value format. The output includes compound name, synonyms, ChEBI ID, PubChem ID, KEGG ID, CAS ID, types, molecular formula, SMILES code and InChI
}
\usage{
parse.metacyc.compound(datPath)
}
\note{
BioCyc and EcoCyc also applicable
}
\arguments{
  \item{datPath}{
Path to compounds.dat file from MetaCyc (BioCyc and EcoCyc also applicable)
}
}
\value{
data.frame() object where data fields include
  \item{UNIQUE.ID }{}
  \item{TYPES }{}
  \item{COMMON.NAME }{}
  \item{CHEMICAL.FORMULA }{}
  \item{INCHI }{}
  \item{MOLECULAR.WEIGHT }{}
  \item{MONOISOTOPIC.MW }{}
  \item{NON.STANDARD.INCHI }{}
  \item{SMILES }{}
  \item{ATOM.CHARGES }{}
  \item{SYNONYMS }{}
  \item{CREDITS }{}
  \item{SUPERATOMS }{}
  \item{REGULATES }{}
  \item{CITATIONS }{}
  \item{COMMENT.INTERNAL }{}
  \item{COMMENT }{}
  \item{SYSTEMATIC.NAME }{}
  \item{HAS.NO.STRUCTURE. }{}
  \item{COMPONENT.OF }{}
  \item{ABBREV.NAME }{}
  \item{TAUTOMERS }{}
  \item{PKA1 }{}
  \item{COFACTORS.OF }{}
  \item{PKA2 }{}
  \item{PKA3 }{}
  \item{N.1.NAME.1 }{}
  \item{N.NAME }{}
  \item{N.1.NAME }{}
  \item{SPECIES }{}
  \item{IN.MIXTURE }{}
  \item{INTERNALS.OF.GROUP }{}
  \item{COMPONENTS }{}
  \item{GROUP.COORDS.2D }{}
  \item{STRUCTURE.GROUPS}{}
  \item{STRUCTURE.LINKS}{}
  \item{RADICAL.ATOMS}{}
  \item{CFG.ICON.COLOR}{}
  \item{GROUP.INTERNALS}{}
  \item{COFACTORS.OR.PROSTHETIC.GROUPS.OFCHEMSPIDER}{}
  \item{PUBCHEM}{}
  \item{CHEBI}{}
  \item{LIGAND.CPD}{}
  \item{NCI}{}
  \item{CAS}{}
  \item{X.Wikipedia.}{}
  \item{UM.BBD.CPD}{}
  \item{KNAPSACK}{}
  \item{BRENDA.COMPOUND}{}
  \item{LIPID_MAPS}{}
  \item{LIGAND}{}
}
\references{
MetaCyc: http://metacyc.org/
}
\author{
Byoungnam Min <mbnmbn00@gmail.com>, Byeonghyuk Park, Kyoung Heon Kim and In-Geol Choi 
}
\examples{
## parse.metacyc.compound: Not run
#url = 'http://bioinformatics.ai.sri.com/ptools/flatfile-samples/compounds.dat'
#tmpdest = tempfile(pattern = "compounds")
#download.file(url, destfile=tmpdest)

#sample.compounds = parse.metacyc.compound(datPath=tmpdest)
}
\keyword{ MetaCyc }
