Rhea2cName <-
function(parsed_ChEBI, eq) {
  h.name = hash()
  .set(h.name, keys = parsed_ChEBI[,'ChEBI'], values = parsed_ChEBI[,'name'])
  result = conv(eq, h.name)
  
  return(result)
}
