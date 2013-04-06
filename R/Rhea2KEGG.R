Rhea2KEGG <-
function(parsed_ChEBI, eq) {
  h.KEGG = hash()
  ind = grep('^$', parsed_ChEBI[,'KEGG'])
  .set(h.KEGG, keys = parsed_ChEBI[-ind,'ChEBI'], values = parsed_ChEBI[-ind,'KEGG'])
  result = conv(eq, h.KEGG)
  
  return(result)
}
