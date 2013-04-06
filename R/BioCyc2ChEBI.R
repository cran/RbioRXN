BioCyc2ChEBI <-
function(parsed_MetaCyc.c, eq) {
  h.ChEBI = hash()
  ind = grep('^$', parsed_MetaCyc.c[,'ChEBI'])
  .set(h.ChEBI, keys = parsed_MetaCyc.c[-ind,'BioCyc'], values = parsed_MetaCyc.c[-ind,'ChEBI'])
  result = conv(eq, h.ChEBI)
  
  return(result)
}
