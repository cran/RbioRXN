Rhea.massBal <-
function(parsed_ChEBI, equation) {
  regexp = "(^\\(?[0-9]*n?\\+?[0-9]*\\)? )(.*)"
  
  R.participant = unlist(strsplit(equation, " => "))
  R.participant = unlist(strsplit(R.participant, " <=> "))
  R.participant = unlist(strsplit(R.participant, " <\\?> "))
  R.participant = unlist(strsplit(R.participant, ' \\+ '))
  R.participant = sub(regexp, '\\2', R.participant)
  R.participant = gsub('\\(.+\\)', "", R.participant)
  R.participant = unique(R.participant)
  R.participant = R.participant[which(R.participant %in% parsed_ChEBI$ChEBI)]
  
  h.formula = hash()
  for(i in 1:length(R.participant)) {
    if(grepl('^$', parsed_ChEBI[parsed_ChEBI$ChEBI == R.participant[i],'InChI']) == FALSE) {
      regexp.p = '(.*)(/p[-+]*[0-9]*)(.*)'
      regexp.h = '(.*)(H[0-9]*)(.*)'
      if(grepl(regexp.p, parsed_ChEBI[parsed_ChEBI$ChEBI == R.participant[i],'InChI'])) {
        tmp = sub(regexp.p, '\\2', parsed_ChEBI[parsed_ChEBI$ChEBI == R.participant[i],'InChI'])
        tmp2 = unlist(strsplit(tmp, '/p'))[2]
        tmp3 = as.numeric(tmp2)
        form = unlist(strsplit(parsed_ChEBI[parsed_ChEBI$ChEBI == R.participant[i],'InChI'], '/'))[2]
        form.h = sub(regexp.h, '\\2', form)
        form.h2 = unlist(strsplit(form.h, 'H'))[2]
        if(length(form.h2) > 0) {
          n = as.numeric(form.h2)
        } else {
          n = 1
        }
        n2 = n + tmp3
        tmp = paste(sub(regexp.h, '\\1', form), 'H', n2, sub(regexp.h, '\\3', form), sep = '')
      } else {
        tmp = unlist(strsplit(parsed_ChEBI[parsed_ChEBI$ChEBI == R.participant[i],'InChI'], '/'))[2]
      }
      .set(h.formula, keys = parsed_ChEBI[parsed_ChEBI$ChEBI == R.participant[i],'ChEBI'], values = tmp)
    }
  }
  # exception
  h.formula[['15378']] = 'H'
  
  R.formula = conv(equation, h.formula)
  result = mass.bal(R.formula)
  
  return(result)
}
