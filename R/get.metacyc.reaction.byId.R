get.metacyc.reaction.byId <-
function(metacycId) {
    if(length(metacycId) == 0) {
        message('Please enter more than one MetaCyc ID')
    }
  
    result_df = c()
    for(i in metacycId) {
        cat('processing',i,'\n')
        urlBase = 'http://websvc.biocyc.org/META/pathway-biopax?type=3&object=%s'
        url = sprintf(urlBase, i)
    
        h = basicTextGatherer()
        try(curlPerform(url = url, writefunction = h$update))
        if(h$value() != '') {
            biopax = h$value()
            biopax = unlist(strsplit(biopax, '\n'))
		
            parsedBiopax = tryCatch({
              .parse.biopax(biopax)
            }, error = function(cond) {
              message(sprintf('WARNING: RbioRXN could not parse reaction %s. It\'s going to be empty data frame', metacycId))
              MetaCyc = i
              result = data.frame(MetaCyc, stringsAsFactors=F)
              return(result)
            }, warning = function(cond) {
              message(sprintf('WARNING: RbioRXN could not parse reaction %s. It\'s going to be empty data frame', metacycId))
              MetaCyc = i
              result = data.frame(MetaCyc, stringsAsFactors=F)
              return(result)
            })
            result_df = rbind.fill(result_df, data.frame(parsedBiopax, stringsAsFactors=F))
        }
    }
    result_df[is.na(result_df)] = ''
    return(result_df)
}
