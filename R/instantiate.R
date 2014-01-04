instantiate <-
function(id, eq, h.class, h.formula, R) {
  eq = gsub('\\|', '', eq)
  
  # define regular expression
  regexp = "(^\\(?[0-9]*n?\\+?[0-9]*\\)? )(.+)"
  regexp.io = '(.+)([\\(\\[](out)?(in)?\\]?\\)?)'
  regexp.n = "([0-9]*n{1} )(.+)" # n as reaction coefficient
  
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
  
  id.r = c()
  eq.r = c()
  for(i in 1:length(eq)){
    if(grepl(regexp.n, eq[i])) {
      eq.n = c()
      id.n = c()
      for(j in 1:length(direction)) {
        if(con1[i,j]) {
          for(k in 1:3) {
            l = unlist(strsplit(eq[i],split=direction[j]))
            m = unlist(strsplit(l[1], split = " \\+ "))
            n = unlist(strsplit(l[2], split = ' \\+ '))
            m.1 = sub(regexp, "\\2", m)
            n.1 = sub(regexp, "\\2", n)
            for(o in 1:length(m)) { 
              if(grepl(regexp, m[o])) { # process coefficeint
                coefficient = sub(regexp, "\\1", m[o])
                if(coefficient == "n ") {
                  coefficient = paste(k, " ", sep="")
                }
                else {
                  coefficient = gsub('n ', '', coefficient)
                  coefficient = as.numeric(coefficient) * k
                  coefficient = paste(coefficient, " ", sep="")
                }
              }
              else {coefficient = ""}
              m[o] = paste(coefficient, m.1[o], sep="")
            }
            for(o in 1:length(n)) { 
              if(grepl(regexp, n[o])) { # process coefficeint
                coefficient = sub(regexp, "\\1", n[o])
                if(coefficient == "n ") {
                  coefficient = paste(k, " ", sep="")
                }
                else {
                  coefficient = gsub('n ', '', coefficient)
                  coefficient = as.numeric(coefficient) * k
                  coefficient = paste(coefficient, " ", sep="")
                }
              }
              else {coefficient = ""}
              n[o] = paste(coefficient, n.1[o], sep="")
            }
            eq.n = c(eq.n, paste(paste(m, collapse = " + "), paste(n, collapse = " + "), sep = To.direction[j]))
            id.n = c(id.n, id[i])
          }
        }
      }
      eq.r = c(eq.r, eq.n)
      id.r = c(id.r, id.n)
    } else {
      eq.r = c(eq.r, eq[i])
      id.r = c(id.r, id[i])
    }
  }
  
  eq = eq.r
  id = id.r 
  
  con2 = matrix(0, length(eq), length(direction))
  for(i in 1:length(eq)) {
    for(j in 1:length(direction)) {
      con2[,j] = grepl(direction[j], eq)
    }
  }
  
  instant.eq = c()
  instant.id = c()
  for(i in 1:length(eq)) {
    for(j in 1:length(direction)) {
      if(con2[i,j]) {
        ll = unlist(strsplit(eq[i],split=direction[j]))
        m = unlist(strsplit(ll[1], split = " \\+ "))
        n = unlist(strsplit(ll[2], split = " \\+ "))
        m2 = sub(regexp, '\\2', m)
        n2 = sub(regexp, '\\2', n)
        m2 = sub(regexp.io, '\\1', m2)
        n2 = sub(regexp.io, '\\1', n2)
        part = c(m2, n2)
        if(length(which(part[which(part %in% R)] %in% names(h.class))) > 0
           && length(which(part[which(part %in% R)] %in% names(h.class) == FALSE)) == 0) {
          c.class = part[which(part %in% R)][which(part[which(part %in% R)] %in% names(h.class))]
          reactant = list()
          product = list()
          for(o in 1:length(m)) { 
            if(grepl(regexp, m[o])) { # process coefficeint
              coefficient = sub(regexp, "\\1", m[o])
            }  else {coefficient = ""}
            if(grepl(regexp.io, m[o])) { 
              in.out = sub(regexp.io, "\\2", x = m[o])
            }  else {in.out = ""}
            if(m2[o] %in% c.class) { 
              tmp = h.class[[m2[o]]]          
              tmp = paste(coefficient, tmp, in.out, sep="")
              reactant[[o]] = tmp
            }  else {
              reactant[[o]] = m[o]
            }
          }
          for(o in 1:length(n)) { 
            if(grepl(regexp, n[o])) { # process coefficeint
              coefficient = sub(regexp, "\\1", n[o])
            }  else {coefficient = ""}
            if(grepl("\\(.+\\)", n[o])) { # process (in), (out)
              in.out = sub(regexp.io, "\\2", x = n[o])
            }  else {
              in.out = ""
            }
            if(n2[o] %in% c.class) { 
              tmp = h.class[[n2[o]]]          
              tmp = paste(coefficient, tmp, in.out, sep="")
              product[[o]] = tmp
            }  else {
              product[[o]] = n[o]
            }
          }
          if(length(reactant) > 1) {
            len = length(reactant) - 1
            for(k in 1:len) {
              reactant[[1]] = outer(reactant[[1]], reactant[[k+1]], paste, sep=" + ")
              dim(reactant[[1]]) = NULL
            }
          }
          if(length(product) > 1) {
            len = length(product) - 1
            for(k in 1:len) {
              product[[1]] = outer(product[[1]], product[[k+1]], paste, sep=" + ")
              dim(product[[1]]) = NULL
            }
          }
          temp = outer(reactant[[1]], product[[1]], paste, sep = To.direction[j])
          dim(temp) = NULL
          
          instant.for = conv(temp, h.formula)
          ind = c()
          for(l in 1:length(instant.for)) {
            if(instant.for[l] != 'cannot convert this equation because of incomplete hash') {
              tmp = mass.bal(instant.for[l])
              if(tmp == TRUE) {
                ind = c(ind, l)
              }
            }
          }
          if(length(ind) > 0) {
            instant.eq = c(instant.eq, temp[ind])
            for(l in 1:length(temp[ind])) {
              instant.id = c(instant.id, paste(id[i], l, sep="_"))
            }
          } 
        }
        else {
          instant.id = c(instant.id, id[i])
          instant.eq = c(instant.eq, eq[i])
        }
      }
    }
  }
  result = list()
  result[['ID']] = instant.id
  result[['Equation']] = instant.eq
  
  return(result)
}
