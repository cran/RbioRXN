BioCyc.massBal <-
function(parsed_MetaCyc.c, equation) {
  h.formula = hash()
  ### edit moluecular formula
  parsed_MetaCyc.c[,'molecular formula'] = gsub('CL', 'Cl', parsed_MetaCyc.c[,'molecular formula']) 
  parsed_MetaCyc.c[,'molecular formula'] = gsub('BR', 'Br', parsed_MetaCyc.c[,'molecular formula'])
  parsed_MetaCyc.c[,'molecular formula'] = gsub('MG', 'Mg', parsed_MetaCyc.c[,'molecular formula'])
  parsed_MetaCyc.c[,'molecular formula'] = gsub('AS', 'As', parsed_MetaCyc.c[,'molecular formula'])
  parsed_MetaCyc.c[,'molecular formula'] = gsub('COBALT', 'Co', parsed_MetaCyc.c[,'molecular formula'])
  parsed_MetaCyc.c[,'molecular formula'] = gsub('FE', 'Fe', parsed_MetaCyc.c[,'molecular formula'])
  parsed_MetaCyc.c[,'molecular formula'] = gsub('SE', 'Se', parsed_MetaCyc.c[,'molecular formula'])
  parsed_MetaCyc.c[,'molecular formula'] = gsub('CR', 'Cr', parsed_MetaCyc.c[,'molecular formula'])
  parsed_MetaCyc.c[,'molecular formula'] = gsub('NA', 'Na', parsed_MetaCyc.c[,'molecular formula'])
  parsed_MetaCyc.c[,'molecular formula'] = gsub('MO', 'Mo', parsed_MetaCyc.c[,'molecular formula'])
  
  ind = grep("^$", parsed_MetaCyc.c[,'molecular formula'])
  if(length(ind) > 0) {
    h.formula[parsed_MetaCyc.c[-ind,'BioCyc']] = parsed_MetaCyc.c[-ind,'molecular formula']
  } else {
    h.formula[parsed_MetaCyc.c[,'BioCyc']] = parsed_MetaCyc.c[,'molecular formula']
  }
  
  
  R.formula = conv(equation, h.formula)
  result = mass.bal(R.formula)
  
  return(result)
}
