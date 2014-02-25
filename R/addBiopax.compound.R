addBiopax.compound = function(df, compound, parsed_ChEBI) {
	# compound_1
	for(i in compound) {
		class = "smallMolecule"
		id = paste("#compound:", i, sep='')
		property = "NAME"
		property_attr = "rdf:datatype"
		property_attr_value = "http://www.w3.org/2001/XMLSchema#string"
		property_value = parsed_ChEBI[parsed_ChEBI$ChEBI == i, 'name']
		
		compound_1 = c(class, id, property, property_attr, property_attr_value, property_value)
		df = rbind(df, compound_1)
	}
	
	# Compound_2
	for(i in compound) {
		class = "smallMolecule"
		id = paste("#compound:", i, sep='')
		property = "DATA-SOURCE"
		property_attr = "rdf:resource"
		property_attr_value = "#dataSource:ChEBI"
		property_value = ""
		
		compound_2 = c(class, id, property, property_attr, property_attr_value, property_value)
		df = rbind(df, compound_2)
	}

	# Compound_3
	for(i in compound) {
		class = "smallMolecule"
		id = paste("#compound:", i, sep='')
		property = "CHEMICAL-FORMULA"
		property_attr = "rdf:datatype"
		property_attr_value = "http://www.w3.org/2001/XMLSchema#string"
		property_value = parsed_ChEBI[parsed_ChEBI$ChEBI == i, 'formula']
		
		compound_3 = c(class, id, property, property_attr, property_attr_value, property_value)
		df = rbind(df, compound_3)
	}
	
	# Compound_4
	for(i in compound) {
		class = "smallMolecule"
		id = paste("#compound:", i, sep="")
		property = "XREF"
		property_attr = "rdf:resource"
		property_attr_value = paste("#CHEBI:", i, sep="")
		property_value = ""
		
		compound_4 = c(class, id, property, property_attr, property_attr_value, property_value)
		df = rbind(df, compound_4)
	}
	
	colnames(df) = c("class", "id", "property", "property_attr", "property_attr_value", "property_value")
	rownames(df) = 1:nrow(df)
	return(df)
}

