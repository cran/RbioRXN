parse.MetaCyc.r <-
function(file) {
  txt = readLines(file)
  
  ## remove HTML code
  txt = gsub("&alpha;", "alpha", txt)
  txt = gsub("&beta;", "beta", txt)
  txt = gsub("&delta;", "delta", txt)
  txt = gsub("&Delta;", "Delta", txt)
  txt = gsub("&gamma;", "gamma", txt)
  txt = gsub("&omega;", "omega", txt)
  txt = gsub("&mdash;", "-", txt)
  txt = gsub("&Psi;", "Psi", txt)
  txt = gsub("&psi;", "psi", txt)
  txt = gsub("&chi;", "chi", txt)
  txt = gsub("&xi;", "xi", txt)
  txt = gsub("&zeta;", "zeta", txt)
  txt = gsub("&harr;", "<->", txt)
  txt = gsub("&pi;", "pi", txt)
  txt = gsub("&tau;", "tau", txt)
  txt = gsub("&lambda;", "lambda", txt)
  txt = gsub("&amp;", "&", txt)
  txt = gsub("&kappa;", "kappa", txt)
  txt = gsub("&iota;", "iota", txt)
  txt = gsub("<sup>", "", txt)
  txt = gsub("</sup>", "", txt)
  txt = gsub("<sub>", "", txt)
  txt = gsub("</sub>", "", txt)
  txt = gsub("<SUP>", "", txt)
  txt = gsub("</SUP>", "", txt)
  txt = gsub("<SUB>", "", txt)
  txt = gsub("</SUB>", "", txt)
  txt = gsub("<i>", "", txt)
  txt = gsub("<I>", "", txt)
  txt = gsub("</I>", "", txt)
  txt = gsub("</i>", "", txt)
  txt = gsub("<em>", "", txt)
  txt = gsub("</em>", "", txt)
  txt = gsub("<small>", "", txt)
  txt = gsub("</small>", "", txt)
  txt = gsub("&larr;", "<-", txt)
  txt = gsub("&rarr;", "->", txt)
  txt = gsub("&epsilon;", "epsilon", txt)
  
  right = c()
  left = c()
  MetaCyc.r = c()
  pathway = c()
  unique.id = ""
  ec.num = ""
  pathway = c()
  eq = ""
  name = ""
  
  con1 = grepl("^UNIQUE-ID", txt)
  con2 = grepl("^COMMON-NAME - ", txt)
  con3 = grepl("^EC-NUMBER - ", txt)
  con4 = grepl("^IN-PATHWAY - ", txt)
  con5 = grepl("^LEFT - ", txt)
  con6 = grepl("^\\^COMPARTMENT", txt)
  con7 = grepl("^\\^COEFFICIENT", txt)
  con8 = grepl("^RIGHT - ",txt)
  con9 = grepl("^REACTION-DIRECTION - ", txt)
  con10 = grepl("^//$", txt)
  
  for(i in 1:length(txt)) {
    if(con1[i]) {
      unique.id = unlist(strsplit(txt[i], " - "))[2]
    }
    else if(con2[i]) {
      name = unlist(strsplit(txt[i], " - "))[2]
    }
    else if(con3[i]) {
      ec.num = unlist(strsplit(txt[i], " - "))[2]
    }
    else if(con4[i])  {
      pathway = c(pathway, unlist(strsplit(txt[i], " - "))[2])
    }
    else if(con5[i]) {
      temp = unlist(strsplit(txt[i], " - "))[2]
      if(con6[i+1] || con6[i+2]) { # compartment
        j = grep("^\\^COMPARTMENT", c(txt[i+1],txt[i+2]))
        in.out = unlist(strsplit(txt[i+j], " - CCO-"))[2]
        in.out = paste("[", tolower(in.out), "]", sep="")
        temp = paste(temp, in.out, sep="")
        if(con7[i+1]) {
          coeff = unlist(strsplit(txt[i+1], " - "))[2]
          temp = paste(coeff, temp, sep=" ")
        }
        left = c(left, temp)
      }
      else {
        left = c(left, temp)
      }
    }
    else if(con8[i]) {
      temp = unlist(strsplit(txt[i], " - "))[2]
      if(con6[i+1] || grepl("^\\^COMPARTMENT", txt[i+2])) { # compartment
        j = grep("^\\^COMPARTMENT", c(txt[i+1], txt[i+2]))
        in.out = unlist(strsplit(txt[i+j], " - CCO-"))[2]
        in.out = paste("[", tolower(in.out), "]", sep="")
        temp = paste(temp, in.out, sep="")
        if(con7[i+1]) {
          coeff = unlist(strsplit(txt[i+1], " - "))[2]
          temp = paste(coeff, temp, sep=" ")
        }
        right = c(right, temp)
      }
      else {
        right = c(right, temp)
      }
    }
    else if(con9[i]) {
      direction = unlist(strsplit(txt[i], " - "))[2]
    }
    else if(con10[i]) {
      if(grepl("LEFT-TO-RIGHT", direction)) {
        left2 = paste(left, collapse=" + ")
        right2 = paste(right, collapse=" + ")
        eq = paste(left2, right2, sep=" -> ")
      }
      else if(grepl("RIGHT-TO-LEFT", direction)) {
        left2 = paste(left, collapse=" + ")
        right2 = paste(right, collapse=" + ")
        eq = paste(right2, left2, sep=" -> ")
      }
      else if(grepl("REVERSIBLE", direction)) {
        left2 = paste(left, collapse=" + ")
        right2 = paste(right, collapse=" + ")
        eq = paste(right2, left2, sep=" <=> ")
      }
      intermediate = c(unique.id, name, ec.num, paste(pathway, collapse=", "), eq)
      MetaCyc.r = rbind(MetaCyc.r, intermediate)
      name = ""
      left = c() # initialization
      right = c()
      unique.id = ""
      ec.num = ""
      pathway = c()
      eq = ""
    }
  }
  ind1 = grep('^ <=> $', MetaCyc.r[,5]) # 2 reactions have no information.
  ind1 = c(ind1, grep('^ -> $', MetaCyc.r[,5])) # 2 reactions have no information.
  MetaCyc.r = MetaCyc.r[-ind1,] # 10,649 reactions
  
  colnames(MetaCyc.r) = c("ID", "Enzyme", "E.C number", "Pathway", "Equation")
  rownames(MetaCyc.r) = 1:nrow(MetaCyc.r)
  MetaCyc.r = as.data.frame(MetaCyc.r, stringsAsFactors=FALSE)
  return(MetaCyc.r)
}
