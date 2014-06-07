\name{check.mass.balance}
\alias{check.mass.balance}
\title{
Check given biochemical reaction is mass-balanced
}
\description{
Using molecular formula of reaction participants, check given reaction equation is mass-balanced. Note that some reaction participants don't have molecular formula (e.g., polymer or generic compound), so these reactions cannot be checked.
}
\usage{
check.mass.balance(equation, chemical_table, id_col="chebi", formula_col="formula", 
direction_type=c(" <=> ", " => ", " <\\\\?> "))
}
\arguments{
  \item{equation}{
Reaction equation usually from get.rhea.all() or get.metacyc.all() function.
}
  \item{chemical_table}{
Chemial table (data.frame) containing molecular formula
}
  \item{id_col}{
In chemical table, choose the column name which contains chemical ID (usually first column)
}
  \item{formula_col}{
In chemical table, choose the column name which contains molecular formula
}
  \item{direction_type}{
In reaction equation, provide reaction direction symbols. Default value is c(" <=> ", " => ", " <\\\\?> ")
}
}
\value{
Character values (not logical)
  \item{ TRUE }{}
  \item{ FALSE }{}
  \item{ Not available }{}
}
\author{
Byoungnam Min <mbnmbn00@gmail.com>, Byeonghyuk Park, Kyoung Heon Kim and In-Geol Choi 
}
\examples{
data(example)

chebi = example$chebi # sample ChEBI

## Check mass balance for example reactions
Rm = example$rhea_massbal # sample Rhea data
check.mass.balance(Rm, chebi, "chebi", "formula", c(" <\\\\?> ", " <=> ", " => "))
}
\keyword{ Biochemical reaction mass balance }
