BioCyc2cName <-
function(parsed_MetaCyc.c, eq) {
  h.name = hash()
  .set(h.name, keys = parsed_MetaCyc.c[,'BioCyc'], values = parsed_MetaCyc.c[,'name'])
  result = conv(eq, h.name)
  
  return(result)
}
