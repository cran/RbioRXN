get.ChEBI <-
function() {
  url = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl'
  tmpdest = tempfile(pattern = "chebi")
  download.file(url, destfile=tmpdest)
  chebi = readLines(tmpdest)
  return(chebi)
}
