BioCyc2ChEBI <-
function(parsed_MetaCyc.c, eq) {
  h.ChEBI = hash()
  ind = grep('^$', parsed_MetaCyc.c[,'ChEBI'])
  h.ChEBI[parsed_MetaCyc.c[-ind,'BioCyc']] = parsed_MetaCyc.c[-ind,'ChEBI']
  result = conv(eq, h.ChEBI)
  
  return(result)
}
