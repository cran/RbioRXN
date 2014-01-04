conv <-
function(eq, h) {
  regexp = "(^\\(?[0-9]*n?\\+?[0-9]*\\)? )(.+)" # regular expression
  regexp.io = '(.+)([\\(\\[](out)?(in)?\\]?\\)?)'
  result = numeric(length(eq))
  
  # Direction
  candidate = c(' <=> ', ' => ', ' <\\?> ', ' -> ')
  direction = c()
  for(i in 1:length(candidate)) {
    if(length(which(grepl(candidate[i], eq))) > 0) {
      direction = c(direction, candidate[i])
    }
  }
  To.direction = gsub('\\\\', '', direction)
  
  con1 = matrix(0, length(eq), length(direction))
  for(i in 1:length(direction)) {
    con1[,i] = grepl(direction[i], eq)
  }
  
  for(i in 1:length(eq)) {
    for(j in 1:length(direction)) {
      if(con1[i,j]) {
        l = unlist(strsplit(eq[i], split=direction[j]))
        m = unlist(strsplit(l[1], split = " \\+ "))
        n = unlist(strsplit(l[2], split = " \\+ "))
        m2 = sub(regexp, '\\2', m)
        n2 = sub(regexp, '\\2', n)
        m2 = sub(regexp.io, '\\1', m2)
        n2 = sub(regexp.io, '\\1', n2)
        if(length(which(!c(m2,n2) %in% names(h)))==0) {
          for(o in 1:length(m)) { 
            if(grepl(regexp, m[o])) { # process coefficeint
              coefficient = sub(regexp, "\\1", m[o])
            }
            else {coefficient = ""}
            if(grepl(regexp.io, m[o])) { # process (in), (out)
              in.out = sub(regexp.io, "\\2", x = m[o])
            }  else {in.out = ""}
            tmp = h[[m2[o]]]
            m[o] = paste(coefficient, tmp, in.out, sep="")
          }
          n.1 = sub(regexp, "\\2", n)
          for(o in 1:length(n)) { 
            if(grepl(regexp, n[o])) { # process coefficeint
              coefficient = sub(regexp, "\\1", n[o])
            }
            else {coefficient = ""}
            if(grepl("\\(.+\\)", n[o])) { # process (in), (out)
              in.out = sub(regexp.io, "\\2", x = n[o])
            }
            else {in.out = ""}
            tmp = h[[n2[o]]]
            n[o] = paste(coefficient, tmp, in.out, sep="")
          }
          result[i] = paste(paste(m, collapse = " + "), paste(n, collapse = " + "), sep = To.direction[j])
        } else {
          result[i] = 'cannot convert this equation because of incomplete hash'
        }
      }
    }
  }
  return(result)
}
