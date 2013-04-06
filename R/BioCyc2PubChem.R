BioCyc2PubChem <-
function(parsed_MetaCyc.c, eq) {
  h.PubChem = hash()
  ind = grep('^$', parsed_MetaCyc.c[,'PubChem'])
  .set(h.PubChem, keys = parsed_MetaCyc.c[-ind,'BioCyc'], values = parsed_MetaCyc.c[-ind,'PubChem'])
  result = conv(eq, h.PubChem)
  
  return(result)
}
