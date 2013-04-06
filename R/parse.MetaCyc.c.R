parse.MetaCyc.c <-
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
  
  con1 = grepl("UNIQUE-ID - ", txt)
  con2 = grepl("DBLINKS - \\(CHEBI", txt)
  con3 = grepl("DBLINKS - \\(PUBCHEM", txt)
  con4 = grepl("DBLINKS - \\(LIGAND-CPD", txt)
  con5 = grepl("DBLINKS - \\(CAS", txt)
  con6 = grepl("TYPES - ", txt)
  con7 = grepl('CHEMICAL-FORMULA - \\(', txt)
  con8 = grepl('SMILES - ', txt)
  con9 = grepl('COMMON-NAME - ', txt)
  con10 = grepl('SYNONYMS - ', txt)
  con11 = grepl('INCHI - ', txt)
  con12 = grepl("^//$", txt)
  
  MetaCyc.c = c()
  cas = ""
  name = ""
  types = c()
  formula = c()
  synonyms = c()
  formula.result = ""
  chebi = ""
  pubchem = ""
  InChI = ""
  
  for(i in 1:length(txt)) {
    if(con1[i]) {
      regexp = "(UNIQUE-ID - )(.+)"
      uniqueId = sub(regexp, "\\2", x = txt[i])
    }
    else if(con2[i]) {
      regexp = '(.+ ")(.+)(".+)'
      chebi = sub(regexp, "\\2", x=txt[i])
    }
    else if(con3[i]) {
      regexp = '(.+ ")(.+)(" .+)'
      pubchem = sub(regexp, "\\2", x = txt[i])
    }
    else if(con4[i]) {
      regexp = '(.+ ")(.+)(".+)'
      kegg = sub(regexp, "\\2", x = txt[i])
    }
    else if(con5[i]) {
      regexp = '(.+ ")(.+)(".+)'
      cas = sub(regexp, "\\2", x = txt[i])        
    }
    else if(con6[i]) {
      regexp = '(TYPES - )(.+)'
      types = c(types, sub(regexp, '\\2', x = txt[i]))
    }
    else if(con7[i]) {
      regexp = '(.+)( - \\()(.+)( )([0-9]+)(\\))'
      tem1 = sub(pattern = regexp, replacement = '\\3', x = txt[i])
      tem2 = sub(pattern = regexp, replacement = '\\5', x = txt[i])
      tem3 = c(tem1, tem2)
      formula = rbind(formula, tem3)
    }
    else if(con8[i]) {
      regexp = '(SMILES - )(.+)'
      SMILES = sub(regexp, '\\2', x = txt[i])
    }
    else if(con9[i]) {
      regexp = '(COMMON-NAME - )(.+)'
      name = sub(regexp, '\\2', x = txt[i])
    }
    else if(con10[i]) {
      regexp = '(SYNONYMS - )(.+)'
      synonyms = c(synonyms, sub(regexp, '\\2', x = txt[i]))
    }
    else if(con11[i]) {
      regexp = '(INCHI - )(.+)'
      InChI = sub(regexp, '\\2', x = txt[i])
    }
    else if(con12[i]) {
      if(length(formula) > 2) {
        formula = formula[order(formula[,1]),]
        for(j in 1:length(formula[,1])) {
          formula.part = paste(formula[j,], collapse = "")
          formula.result = paste(formula.result, formula.part, sep="")
        }
      }
      else if(length(formula) == 2) {
        formula.result = paste(formula, collapse="")
      }
      temp = c(uniqueId, name, paste(synonyms, collapse='///'), chebi, pubchem, kegg, cas, paste(types, collapse='///'), formula.result, SMILES, InChI)
      MetaCyc.c = rbind(MetaCyc.c, temp)
      formula.result = ""
      formula = c()
      name = ""
      synonyms = c()
      types = c()
      temp = c()
      uniqueId = ""
      chebi = ""
      pubchem = ""
      kegg = ""
      cas = ""
      SMILES = ""
    }
  }
  
  colnames(MetaCyc.c) = c("BioCyc", 'name', 'synonyms', "ChEBI", "PubChem", "KEGG", "CAS", "TYPES", "molecular formula", "SMILES", 'InChI')
  rownames(MetaCyc.c) = 1:nrow(MetaCyc.c)
  MetaCyc.c = as.data.frame(MetaCyc.c, stringsAsFactors=FALSE)
  
  return(MetaCyc.c)
}
