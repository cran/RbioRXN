mass.bal <-
function(eq) {
  eval = rep('FALSE', length(eq)) # result
  
  ind = which(grepl('cannot convert this equation because of incomplete hash', eq) == FALSE)
  eval[ind] == 'cannot be checked'

  ## define regular exp
  regexp = "(^\\(?[0-9]*n?\\+?[0-9]*\\)? )(.+)"
  regexp.io = '(.+)([\\(\\[](out)?(in)?\\]?\\)?)'
  
  ## remove localization tag
  eq = gsub('\\(in\\)', '', eq)
  eq = gsub('\\[in\\]', '', eq) 
  eq = gsub('\\(out\\)', '', eq)
  eq = gsub('\\[out\\]', '', eq)
  
  reactant.for = list()
  product.for = list()
  for(i in ind) {
    if(grepl('cannot convert this equation because of incomplete hash', eq[i])) {
      evel[i] = 'cannot be checked'
    }      
      l = unlist(strsplit(eq[i], split=" <=> "))
      l = unlist(strsplit(l, split=" => "))
      l = unlist(strsplit(l, split=" -> "))
      l = unlist(strsplit(l, split=' <\\?> '))
      m = unlist(strsplit(l[1], split = " \\+ "))
      n = unlist(strsplit(l[2], split = " \\+ "))
      m.1 = c()
      for(j in 1:length(m)) {
        if(grepl(regexp, m[j])) {
          coefficient = sub(regexp, '\\1', m[j])
          coefficient = as.numeric(coefficient)
          tmp = sub(regexp, '\\2', m[j])
          m.1 = c(m.1, rep(tmp, coefficient))
        } else {m.1 = c(m.1, m[j])}
      }
      n.1 = c()
      for(j in 1:length(n)) {
        if(grepl(regexp, n[j])) {
          coefficient = sub(regexp, '\\1', n[j])
          coefficient = as.numeric(coefficient)
          tmp = sub(regexp, '\\2', n[j])
          n.1 = c(n.1, rep(tmp, coefficient))
        }
        else {n.1 = c(n.1, n[j])}
      }
      reactant.for[[i]] = m.1
      product.for[[i]] = n.1
    }
    
    for(i in ind) {
      reactant.for[[i]] = makeup(reactant.for[[i]], sum = TRUE)
      product.for[[i]] = makeup(product.for[[i]], sum = TRUE)
    }
    
    for(i in ind) {
      reactant.for[[i]] = as.matrix(reactant.for[[i]])
      reactant.for[[i]] = reactant.for[[i]][order(row.names(reactant.for[[i]])),]
      product.for[[i]] = as.matrix(product.for[[i]])
      product.for[[i]] = product.for[[i]][order(row.names(product.for[[i]])),]
    }
    
    for(i in ind){
      if(length(which(names(reactant.for[[i]]) %in% names(product.for[[i]]) == FALSE)) == 0 && length(which(names(product.for[[i]]) %in% names(reactant.for[[i]])==FALSE)) == 0) {
        for(k in names(reactant.for[[i]])) {
          tmp = reactant.for[[i]][k] == product.for[[i]][k]
          if(tmp == FALSE) {
            break
          }
          else if(k == names(tail(reactant.for[[i]], 1))) {eval[i] = TRUE}
        }
      }
  }
  return(eval)
}
