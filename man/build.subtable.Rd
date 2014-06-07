\name{build.subtable}
\alias{build.subtable}
\title{
Build sub-table for the entries with multiple values
}
\description{
When you download biochemical reaction data with functions such as get.kegg.all(), get.metacyc.all(), get.rhea.all(), etc, some entries have multiple values such as "synonyms". This function allows to build sub-table that separates multiple values in different rows. 
}
\usage{
build.subtable(table, column1, column2, separator = "///")
}
\arguments{
  \item{table}{
data.frame object where multiple entries exist
}
  \item{column1}{
It is usually first column name for unique ID
}
  \item{column2}{
The column name where the column contains multiple values (e.g., synonyms)
}
  \item{separator}{
Separator symbol (default is "///")
}
}
\value{
data.frame() object
}
\author{
Byoungnam Min <mbnmbn00@gmail.com>, Byeonghyuk Park, Kyoung Heon Kim and In-Geol Choi
}
\examples{
data(example) # load example data

chebi = example$chebi # sample ChEBI
chebi.sub = build.subtable(table = chebi, column1 = "chebi", column2 = "synonyms", separator="///")
head(chebi.sub)
}
