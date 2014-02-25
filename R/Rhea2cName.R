Rhea2cName <-
function(parsed_ChEBI, eq) {
  h.name = hash()
  h.name[parsed_ChEBI[,'ChEBI']] = parsed_ChEBI[,'name']
  result = conv(eq, h.name)
  
  return(result)
}
