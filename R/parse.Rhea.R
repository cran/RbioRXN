parse.Rhea <-
function(owl) {
  
  rhea.reaction = c()
  
  con1 = grepl('<bp:biochemicalReaction rdf:about="#', owl)
  con2 = grepl("<bp:transport", owl)
  con3 = grepl(pattern = "<bp:NAME rdf:datatype", owl)
  con4 = grepl('</bp:biochemicalReaction>', owl)
  con5 = grepl('</bp:transport>', owl)
  con6 = grepl('</bp:transportWithBiochemicalReaction>', owl)
  
  for(i in 1:length(owl)) {
    if(con1[i] | con2[i]) {
      id = unlist(strsplit(owl[i],split='"'))[2]
      id = gsub('#', '', id)
    }
    else if(con3[i]) {
      regexp = '(<bp:NAME rdf:datatype = "http://www.w3.org/2001/XMLSchema#string">)(.*)(</bp:NAME>)'
      equation = sub(pattern = regexp, replacement = "\\2", x = owl[i])
    }
    else if(con4[i] | con5[i] | con6[i]) {
      int = c(id, equation)
      rhea.reaction = rbind(rhea.reaction, int)
      id = ""
      equation = ""
    }
  }

  rhea.reaction = trim(rhea.reaction)
  
  # replace html code
  rhea.reaction[,2] = gsub(pattern = "&gt;", replacement = ">", x = rhea.reaction[,2])
  rhea.reaction[,2] = gsub(pattern = "&lt;", replacement = "<", x = rhea.reaction[,2])
  rhea.reaction[,2] = gsub(pattern = "&apos;", replacement = "'", x = rhea.reaction[,2])
  
  # Sorting
  rhea.reaction = rhea.reaction[order(rhea.reaction[,1]),]
  colnames(rhea.reaction) = c("ID", "Equation")
  rownames(rhea.reaction) = 1:nrow(rhea.reaction)
  
  h.ChEBI = hash() # ChEBI

  con1 = grepl('<bp:smallMolecule rdf:about="#compound:', owl)
  con2 = grepl(' <bp:NAME rdf:datatype = "http://www.w3.org/2001/XMLSchema#string">', owl)
  con3 = grepl(' <bp:CHEMICAL-FORMULA rdf:datatype = "http://www.w3.org/2001/XMLSchema#string">', owl)
  con4 = grepl('<bp:XREF rdf:resource="#CHEBI:', owl)
  con5 = grepl('</bp:smallMolecule>', owl)
  
  for(i in 1:length(owl)) {
    if(con1[i]) {
      regexp = '(<bp:smallMolecule rdf:about="#compound:)(.+)(">)'
      id = sub(pattern = regexp, replacement = "\\2", x = owl[i])
    }
    else if(con2[i]) {
      regexp = '( <bp:NAME rdf:datatype = "http://www.w3.org/2001/XMLSchema#string">)(.*)(</bp:NAME>)'
      name = sub(pattern = regexp, replacement = '\\2', x = owl[i])
      name = gsub("&apos;","'",name)
      name = gsub("&gt;", ">",name)
      name = gsub("&lt;", "<",name)
    }
    else if (con3[i]) {
      regexp = '( <bp:CHEMICAL-FORMULA rdf:datatype = "http://www.w3.org/2001/XMLSchema#string">)(.*)(</bp:CHEMICAL-FORMULA>)'
      formula = sub(pattern = regexp, replacement = '\\2', x = owl[i])
    }
    else if(con4[i]) {
      regexp = '( <bp:XREF rdf:resource="#CHEBI:)(.*)(" />)'
      chebi = sub(pattern = regexp, replacement = '\\2', x = owl[i])
    }
    else if(con5[i]) {
      h.ChEBI[name] = chebi
    }
  }
  
  # Rhea to ChEBI
  direction = c(" => ", " <=> ", " <\\?> ")
  To.direction = c(" => ", " <=> ", " <?> ")
  
  ## to chebi
  con1.1 = grepl(direction[1], rhea.reaction[,2])
  con1.2 = grepl(direction[2], rhea.reaction[,2])
  con1.3 = grepl(direction[3], rhea.reaction[,2])
  con1 = cbind(con1.1, con1.2, con1.3)
  
  regexp = "(^\\(?[0-9]*n?\\+?[0-9]*\\)? )(.+)"
  regexp2 = "(.+)(\\(.+\\))"
  
  Eq.ChEBI = numeric(length(rhea.reaction[,2]))
  for(i in 1:nrow(rhea.reaction)) {
    for(j in 1:length(direction)) {
      if(con1[i,j]) {
        l = unlist(strsplit(rhea.reaction[i,2], split=direction[j]))
        m = unlist(strsplit(l[1], split = " \\+ "))
        m.1 = sub(regexp, "\\2", m)
        for(o in 1:length(m)) { 
          if(grepl(regexp, m[o])) { # process coefficeint
            coefficient = sub(regexp, "\\1", m[o])
          }
          else {coefficient = ""}
          if(grepl("\\(in\\)", m[o]) || grepl("\\(out\\)", m[o])) { # process (in), (out)
            in.out = sub(regexp2, "\\2", x = m[o])
          }
          else {in.out = ""}
          m.1[o] = gsub("\\(in\\)", "", m.1[o])
          m.1[o] = gsub("\\(out\\)", "", m.1[o])
          tmp = h.ChEBI[[m.1[o]]]
          m[o] = paste(coefficient, tmp, in.out, sep="")
        }
        n = unlist(strsplit(l[2], split = " \\+ "))
        n.1 = sub(regexp, "\\2", n)
        for(o in 1:length(n)) { 
          if(grepl(regexp, n[o])) { # process coefficeint
            coefficient = sub(regexp, "\\1", n[o])
          }
          else {coefficient = ""}
          if(grepl("\\(in\\)", n[o]) || grepl("\\(out\\)", n[o])) { # process (in), (out)
            in.out = sub(regexp2, "\\2", x = n[o])
          }
          else {in.out = ""}
          n.1[o] = gsub("\\(in\\)", "", n.1[o])
          n.1[o] = gsub("\\(out\\)", "", n.1[o])
          tmp = h.ChEBI[[n.1[o]]]          
          n[o] = paste(coefficient, tmp, in.out, sep="")
        }
        Eq.ChEBI[i] = paste(paste(m, collapse = " + "), paste(n, collapse = " + "), sep = To.direction[j])
      }  
    }
  }
  tmp = cbind(rhea.reaction, Eq.ChEBI)
  tmp = as.data.frame(tmp, stringsAsFactors=FALSE)
  return(tmp)
}
