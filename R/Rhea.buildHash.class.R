Rhea.buildHash.class <-
function(c.table, R) {
  
  # make hash
  h.class = hash() 
  for(i in 1:length(c.table$parent)){
    tmp = unlist(strsplit(c.table$parent[i], "///"))
    if(length(tmp) > 0) {
      for(j in 1:length(tmp)) {
        if(tmp[j] %in% R) {
          .set(h.class, keys = tmp[j], values = c(h.class[[tmp[j]]], c.table[i,'ChEBI']))
        }
      }  
    }
  }

  check = c(0,1)
  k = 2
  while(check[k] != check[k-1]) {
    count = 0
    for(i in keys(h.class)) {
      for(j in h.class[[as.character(i)]]) {
        tmp = h.class[[as.character(j)]]
        if(j %in% R && length(tmp)>0) {
          h.class[[as.character(i)]] = c(h.class[[as.character(i)]], tmp)
          ind = grep(glob2rx(j), h.class[[as.character(i)]])
          h.class[[as.character(i)]] = h.class[[as.character(i)]][-ind]
        }
        else if(length(tmp) == 0 && c.table[c.table$ChEBI == j,'SMILES'] == "") {
          ind = grep(glob2rx(j), h.class[[as.character(i)]])
          h.class[[as.character(i)]] = h.class[[as.character(i)]][-ind]
        }
        else if(length(tmp) == 0 && grepl('\\*', c.table[c.table$ChEBI == j,'SMILES'])) {
          ind = grep(glob2rx(j), h.class[[as.character(i)]])
          h.class[[as.character(i)]] = h.class[[as.character(i)]][-ind]
        }
      }
      count = count +1
    }
    v.count = values(h.class)
    v.count = unlist(v.count)
    check = c(check, length(v.count))
    k = k+1
  }
  
  for(i in keys(h.class)) {
    h.class[[as.character(i)]] = unique(h.class[[as.character(i)]])
  }
  
  for(i in keys(h.class)) { 
    if(length(h.class[[as.character(i)]]) == 0) {
      del(as.character(i), h.class)
    }
  }
  return(h.class)
}
