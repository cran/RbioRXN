Rhea2KEGG <-
function(parsed_ChEBI, eq) {
  h.KEGG = hash()
  ind = grep('^$', parsed_ChEBI[,'KEGG'])
  h.KEGG[parsed_ChEBI[-ind,'ChEBI']] = parsed_ChEBI[-ind,'KEGG']
  result = conv(eq, h.KEGG)
  
  return(result)
}
