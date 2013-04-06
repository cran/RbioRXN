BioCyc.instantiate <-
function(parsed_MetaCyc.r, parsed_MetaCyc.c, BioCyc_ID, multicore=1) {
  # Make hash on compound class (compound class - instances)
  print('Making hash compound class')
  types = parsed_MetaCyc.c[,'TYPES']
  types = unlist(strsplit(types, '///')) 
  types = unique(types)
  types = types[-1] # remove "Compounds"
  
  h.class = hash() # key-type, value-instatiated compound ID
  
  for(i in 1:length(types)) { 
    ind = grep(types[i], parsed_MetaCyc.c[,'TYPES'])
    .set(h.class, keys = types[i], values = parsed_MetaCyc.c[ind,'BioCyc'])
  }
  
  # Make hash on compound molecular formula
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
  
  ind = grep("^$", parsed_MetaCyc.c[,'molecular formula'])
  if(length(ind) > 0) {
  .set(h.formula, keys = parsed_MetaCyc.c[-ind,'BioCyc'], values = parsed_MetaCyc.c[-ind,'molecular formula'])
  } else {
    .set(h.formula, keys = parsed_MetaCyc.c[,'BioCyc'], values = parsed_MetaCyc.c[,'molecular formula'])
  }
  
  # Instantiate
  print('Instantiating')
  
  R = keys(h.class)  
  
  registerDoMC(multicore)
  tmp = foreach(i = 1:length(BioCyc_ID)) %dopar% {
    instantiate(BioCyc_ID[i], parsed_MetaCyc.r[parsed_MetaCyc.r$ID == BioCyc_ID[i],'Equation'], h.class, h.formula, R)
  }
  
  ID = unlist(tmp)[which(grepl('ID', names(unlist(tmp))))]
  Equation = unlist(tmp)[which(grepl('Equation', names(unlist(tmp))))]
  tmp2 = cbind(ID, Equation)
  rownames(tmp2) = 1:nrow(tmp2)
  tmp3 = as.data.frame(tmp2, stringsAsFactors= FALSE)
  return(tmp3)
}
