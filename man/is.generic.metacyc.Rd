\name{is.generic.metacyc}
\alias{is.generic.metacyc}
\title{
Check given reaction is generic reaction
}
\description{
Generic reaction is a reaction with at least one compound class as reaction participant (e.g., a primary alcohol + NADP -> an aldehyde + H(+) + NADPH). This function check if given reaction is generic reaction
}
\usage{
is.generic.metacyc(equation,compound_df,id_col,smiles_col,parent_col,direction_type)
}
\arguments{
  \item{equation}{
Reaction equation usually from get.metacyc.all() function.
}
  \item{compound_df}{
Chemial table (data.frame) containing SMILES code of each compound in reaction equation
}
  \item{id_col}{
In chemical table, choose the column name which contains chemical id (usually first column)
}
  \item{smiles_col}{
In chemical table, choose the column name which contains SMILES code
}
  \item{parent_col}{
In chemical table, choose the column name which contains TYPES
}
  \item{direction_type}{
In reaction equation, provide reaction direction charaters. Default is c(" <=> ", " => ")
}
}
\value{
Logical TRUE/FALSE
}
\references{
MetaCyc: http://metacyc.org
}
\author{
Byoungnam Min <mbnmbn00@gmail.com>, Byeonghyuk Park, Kyoung Heon Kim and In-Geol Choi 
}
\examples{
data(example) # Load sample data

metacyc.c = example$metacyc.c
metacyc.r.1 = example$metacyc_generic$equation
metacyc.r.2 = example$metacyc_massbal
metacyc.r = c(metacyc.r.1, metacyc.r.2)

names(metacyc.c) # Print column names

# Check generic reaction
is.generic.metacyc(metacyc.r,metacyc.c,'UNIQUE.ID','SMILES','TYPES',c(' <=> ', ' => '))
}
\keyword{ MetaCyc }
\keyword{ Generic reaction }
