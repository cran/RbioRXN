get.metacyc.all <-
function() {
	url = 'http://websvc.biocyc.org/xmlquery?[x:x<-meta^^reactions]'

  message('Download MetaCyc reaction list. It takes a while')
	h = basicTextGatherer()
	curlPerform(url = url, writefunction = h$update)
	
	xml = h$value()
	xml = unlist(strsplit(xml, '\n'))
	
	index = grep('^  <Reaction ID', xml)
	reg_ex = "(.*META:)(.*)(' orgid.+)"
	reactionIds = sub(reg_ex, '\\2', xml[index])
  
  message(sprintf('%s reaction entries are being downloaded', length(reactionIds)))
	metacycDf = .parse.metacyc.biopax(reactionIds)
	return(metacycDf)
}
