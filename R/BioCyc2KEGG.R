BioCyc2KEGG <-
function(parsed_MetaCyc.c, eq) {
  h.KEGG = hash()
  ind = grep('^$', parsed_MetaCyc.c[,'KEGG'])
  h.KEGG[parsed_MetaCyc.c[-ind,'BioCyc']] = parsed_MetaCyc.c[-ind,'KEGG']
  result = conv(eq, h.KEGG)
  
  return(result)
}
