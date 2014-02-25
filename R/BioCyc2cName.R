BioCyc2cName <-
function(parsed_MetaCyc.c, eq) {
  h.name = hash()
  h.name[parsed_MetaCyc.c[,'BioCyc']] = parsed_MetaCyc.c[,'name']
  result = conv(eq, h.name)
  
  return(result)
}
