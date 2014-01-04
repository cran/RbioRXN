BioCyc2KEGG <-
function(parsed_MetaCyc.c, eq) {
  h.KEGG = hash()
  ind = grep('^$', parsed_MetaCyc.c[,'KEGG'])
  .set(h.KEGG, keys = parsed_MetaCyc.c[-ind,'BioCyc'], values = parsed_MetaCyc.c[-ind,'KEGG'])
  result = conv(eq, h.KEGG)
  
  return(result)
}
