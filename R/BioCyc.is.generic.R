BioCyc.is.generic <-
function(parsed_MetaCyc.c, equation) {
  # define regular expression
  regexp = "(^\\(?[0-9]*n?\\+?[0-9]*\\)? )(.*)"
  regexp.io = '(.+)([\\(\\[](out)?(in)?\\]?\\)?)'
  
  # Make hash on compound class (compound class - instances)
  types = parsed_MetaCyc.c[,'TYPES']
  types = unlist(strsplit(types, '///')) 
  types = unique(types)
  types = types[-1] # remove "Compounds"
  
  h.class = hash() # key-type, value-instatiated compound ID
  
  for(i in 1:length(types)) { 
    ind = grep(types[i], parsed_MetaCyc.c[,'TYPES'])
    h.class[types[i]] = parsed_MetaCyc.c[ind,'BioCyc']
  }
  
  R = keys(h.class)  

  # Check if equation is generic
  eq = gsub('\\|', '', equation)
  
  # define direction
  candidate = c(' <=> ', ' => ', ' <\\?> ', ' -> ')
  direction = c()
  for(i in 1:length(candidate)) {
    if(length(which(grepl(candidate[i], eq))) > 0) {
      direction = c(direction, candidate[i])
    }
  }
  To.direction = gsub('\\\\', '', direction)
  
  con1 = matrix(0, length(eq), length(direction))
  for(i in 1:length(eq)) {
    for(j in 1:length(direction)) {
      con1[,j] = grepl(direction[j], eq)
    }
  }
  
  # remove coefficient
  eq2 = c()
  for(i in 1:length(eq)) {
    if(grepl(regexp, eq[i])) {
      for(j in 1:length(direction)) {
        if(con1[i,j]) {
          l = unlist(strsplit(eq[i],split=direction[j]))
          m = unlist(strsplit(l[1], split = " \\+ "))
          n = unlist(strsplit(l[2], split = ' \\+ '))
          m.1 = sub(regexp, "\\2", m)
          m.1 = sub(regexp.io, '\\1', m.1)
          n.1 = sub(regexp, "\\2", n)
          n.1 = sub(regexp.io, '\\1', n.1)
          eq2 = c(eq2, paste(paste(m.1, collapse = " + "), paste(n.1, collapse = " + "), sep = To.direction[j]))
        }
      }
    } else {
      eq2 = c(eq2, eq)
    }
  }
  
  result = c()
  for(i in 1:length(eq2)) {
    for(j in 1:length(direction)) {
      if(con1[i,j]) {
        ll = unlist(strsplit(eq2[i],split=direction[j]))
        m = unlist(strsplit(ll[1], split = " \\+ "))
        n = unlist(strsplit(ll[2], split = " \\+ "))
        m2 = sub(regexp, '\\2', m)
        n2 = sub(regexp, '\\2', n)
        m2 = sub(regexp.io, '\\1', m2)
        n2 = sub(regexp.io, '\\1', n2)
        part = c(m2, n2)
        result = c(result, length(which(part[which(part %in% R)] %in% names(h.class))) > 0 && length(which(part[which(part %in% R)] %in% names(h.class) == FALSE)) == 0)
      }
    }
  }
  return(result)
}
