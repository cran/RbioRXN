BioCyc2PubChem <-
function(parsed_MetaCyc.c, eq) {
  h.PubChem = hash()
  ind = grep('^$', parsed_MetaCyc.c[,'PubChem'])
  h.PubChem[parsed_MetaCyc.c[-ind,'BioCyc']] = parsed_MetaCyc.c[-ind,'PubChem']
  result = conv(eq, h.PubChem)
  
  return(result)
}
