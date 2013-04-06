parse.ChEBI <-
function(owl) {
    
  con1 = grepl('    <owl:Class rdf:about="http://purl.obolibrary.org/obo/CHEBI_',owl)
  con2 = grepl('        <rdfs:label rdf:datatype="http://www.w3.org/2001/XMLSchema#string">', owl)
  con3 = grepl('        <obo2:Synonym rdf:datatype="http://www.w3.org/2001/XMLSchema#string">', owl)
  con4 = grepl('        <obo2:SMILES rdf:datatype="http://www.w3.org/2001/XMLSchema#string">', owl)
  con5 = grepl('        <obo2:InChI rdf:datatype="http://www.w3.org/2001/XMLSchema#string">', owl)
  con6 = grepl('        <obo2:xref rdf:datatype="http://www.w3.org/2001/XMLSchema#string">KEGG COMPOUND:C', owl)
  con7 = grepl('        <rdfs:subClassOf rdf:resource="http://purl.obolibrary.org/obo/CHEBI_', owl)
  con8 = grepl("</owl:Class>",owl)
  
  id = ""
  name = ""
  synonym = c()
  SMILE = ""
  InChi = ""
  KEGG = ""
  ChEBI = c()
  parent = c()
  for (i in 1:length(owl)) {
    if(con1[i]) { # ID
      regexp = '(    <owl:Class rdf:about=\"http://purl.obolibrary.org/obo/CHEBI_)(.*)(">)'
      id = sub(pattern = regexp, replacement = "\\2", x = owl[i])
      id = trim(id)
    }
    else if(con2[i]) { # Name
      regexp = '(        <rdfs:label rdf:datatype=\"http://www.w3.org/2001/XMLSchema#string\">)(.*)(</rdfs:label>)'
      name = sub(pattern = regexp, replacement = "\\2", x = owl[i])
      name = trim(name)
    }
    else if(con3[i]) { # synonym
      regexp = '(        <obo2:Synonym rdf:datatype="http://www.w3.org/2001/XMLSchema#string">)(.*)(</obo2:Synonym>)'
      synonym = c(synonym, sub(pattern = regexp, replacement = '\\2', x = owl[i]))
      synonym = trim(synonym)
    }
    else if(con4[i]) { # SMILE
      regexp = '(        <obo2:SMILES rdf:datatype=\"http://www.w3.org/2001/XMLSchema#string\">)(.*)(</obo2:SMILES>)'
      SMILE = sub(pattern = regexp, replacement = '\\2', x = owl[i])
      SMILE = trim(SMILE)
    }
    else if(con5[i]) { # InChi
      regexp = '(        <obo2:InChI rdf:datatype=\"http://www.w3.org/2001/XMLSchema#string\">)(.*)(</obo2:InChI>)'
      InChi = sub(pattern = regexp, replacement = '\\2', x = owl[i])
      InChi = trim(InChi)
    }
    else if(con6[i]) { # KEGG
      regexp = '(        <obo2:xref rdf:datatype=\"http://www.w3.org/2001/XMLSchema#string\">KEGG COMPOUND:)(.*)(</obo2:xref>)'
      KEGG = sub(pattern = regexp, replacement = '\\2', x = owl[i])
      KEGG = trim(KEGG)
    }
    else if(con7[i]) { # subClassOf
      regexp = '(        <rdfs:subClassOf rdf:resource="http://purl.obolibrary.org/obo/CHEBI_)(.*)("/>)'
      parent = c(parent, sub(pattern = regexp, replacement = '\\2', x = owl[i]))
      parent = trim(parent)
    }
    else if(con8[i]) {
      data = c(id, name, paste(synonym, collapse="///"), SMILE, InChi, KEGG, paste(parent, collapse="///"))
      ChEBI = rbind(ChEBI, data)
      id = ""
      name = ""
      synonym = c()
      parent = c()
      SMILE = ""
      InChi = ""
      KEGG = ""
    }
  }
  
  ChEBI = gsub("&#39;", "'", ChEBI) # convert HTML code
  
  colnames(ChEBI) = c("ChEBI", "name", "synonyms", "SMILES", "InChI", "KEGG", "parent")
  rownames(ChEBI) = 1:nrow(ChEBI)
  ChEBI = as.data.frame(ChEBI, stringsAsFactors=FALSE)
  return(ChEBI)
}
