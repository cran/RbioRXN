get.Rhea <-
function() {
  url = 'ftp://ftp.ebi.ac.uk/pub/databases/rhea/biopax/rhea-biopax.owl.gz'
  tmpdest = tempfile(pattern = "rhea")
  download.file(url, destfile=tmpdest)
  rhea = readLines(tmpdest)
  return(rhea)
}
