\name{RbioRXN-package}
\alias{RbioRXN-package}
\alias{RbioRXN}
\docType{package}
\title{
RbioRXN - Get/process Rhea and MetaCyc (BioCyc, EcoCyc) biochemical reaction data
}
\description{
RbioRXN is an R package that facilitates gaining and processing biochemical reaction data such as Rhea, MetaCyc, BioCyc and EcoCyc. The package provides various functions to download and parse raw data, instantiate generic reaction, convert compound ID into the counterpart ID in other checmial databases (e.g. KEGG, PubChem, and ChEBI) and check mass balance. The package aims to construct an integrated metabolic network and genome-scale metabolic model.
}
\details{
\tabular{ll}{
Package: \tab RbioRXN\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2014-02-24\cr
License: \tab GPL (>=2)\cr
Depends: hash, gdata, CHNOSZ, foreach
}
}
\author{
Byoungnam Min, Kyoung Heon Kim and In-Geol Choi \cr
Maintainer: Byoungnam Min <mbnmbn00@gmail.com>
}
\references{
http://compbio.korea.ac.kr
}
\examples{
## Not run:

# Get/parse Rhea and ChEBI data (full data)
### Rhea = get.Rhea()
### parsed_Rhea = parse.Rhea(Rhea)

### ChEBI = get.ChEBI()
### parsed_ChEBI = parse.ChEBI(ChEBI)

## End(Not run)

# Get/parse MetaCyc reaction and compound data (sample, getting full data needs license)
## compounds.dat 
url = 'http://bioinformatics.ai.sri.com/ptools/flatfile-samples/compounds.dat'
tmpdest = tempfile(pattern = "compounds")
download.file(url, destfile=tmpdest)

sample.compounds = parse.MetaCyc.c(file=tmpdest)

## reactions.dat
url = 'http://bioinformatics.ai.sri.com/ptools/flatfile-samples/reactions.dat'
tmpdest = tempfile(pattern = "reactions")
download.file(url, destfile=tmpdest)

direction = "LEFT-TO-RIGHT" # sample data is old, so don't mind this code.
sample.reactions = parse.MetaCyc.r(file=tmpdest)

# Instantiate
data(example) # load example data

## Rhea instantiation (Rhea.instantiate(parsed_Rhea,parsed_ChEBI,Rhea_ID))
pC = example$parsed_ChEBI # sample ChEBI
Rg = example$Rhea_generic # sample Rhea generic reaction

data(thermo) # this is for the package 'CHNOSZ'
instanceR = Rhea.instantiate(Rg, pC, Rg[1,'ID'])
print(instanceR)

## MetaCyc instantiation (BioCyc.instantiate(parsed_MetaCyc.r,parsed_MetaCyc.c,BioCyc_ID))
pMc = example$parsed_MetaCyc.c # sample MetaCyc compound
Mg = example$MetaCyc_generic # sample generic reaction

instanceM = BioCyc.instantiate(Mg, pMc, Mg[1,'ID'])
print(instanceM)

# Rhea.is.generic (Rhea.is.generic(parsed_ChEBI, equation))
Rhea.is.generic(pC, Rg[,'Eq.ChEBI'])

# BioCyc.is.generic (BioCyc.is.generic(parsed_MetaCyc.c, equation))
BioCyc.is.generic(pMc, Mg[,'Equation'])

# Conversion compound ID into counterpart ID in other chemical database

## Rhea conversion (Rhea2KEGG, Rhea2cName)
Rc = example$Rhea_conv # sample Rhea data
print(Rc)

R2KEGG = Rhea2KEGG(pC, Rc) # ChEBI ID to KEGG ID
print(R2KEGG)

R2cName = Rhea2cName(pC, Rc) # ChEBI ID to compound name
print(R2cName)

## MetaCyc conversion (BioCyc2KEGG, BioCyc2PubChem, BioCyc2ChEBI, BioCyc2cName)
Mc = example$MetaCyc_conv # sample MetaCyc data
print(Mc)

B2KEGG = BioCyc2KEGG(pMc, Mc) # BioCyc ID to KEGG ID
print(B2KEGG)

B2ChEBI = BioCyc2ChEBI(pMc, Mc) # BioCyc ID to ChEBI ID
print(B2ChEBI)

B2PubChem = BioCyc2PubChem(pMc, Mc) # BioCyc ID to PubChem ID
print(B2PubChem)

B2cName = BioCyc2cName(pMc, Mc) # BioCyc ID to compound name
print(B2cName)

# Check mass balance 

## Rhea check mass balance (Rhea.massBal(parsed_ChEBI,equation))
Rm = example$Rhea_massBal # sample Rhea data
checkR = Rhea.massBal(pC, Rm)
print(checkR)

## MetaCyc check mass balance (BioCyc.massBal(parsed_MetaCyc.c, equation))
Mm = example$MetaCyc_massBal # sample Rhea data
checkM = BioCyc.massBal(pMc, Mm)
print(checkM)
}
