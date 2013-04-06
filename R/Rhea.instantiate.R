Rhea.instantiate <-
function(parsed_Rhea, parsed_ChEBI, Rhea_ID, multicore = 1) {
  regexp = "(^\\(?[0-9]*n?\\+?[0-9]*\\)? )(.*)"
  
  R.participant = unlist(strsplit(parsed_Rhea[parsed_Rhea$ID %in% Rhea_ID,'Eq.ChEBI'], " => "))
  R.participant = unlist(strsplit(R.participant, " <=> "))
  R.participant = unlist(strsplit(R.participant, " <\\?> "))
  R.participant = unlist(strsplit(R.participant, ' \\+ '))
  R.participant = sub(regexp, '\\2', R.participant)
  R.participant = gsub('\\(.+\\)', "", R.participant)
  R.participant = unique(R.participant)
  R.participant = R.participant[which(R.participant %in% parsed_ChEBI$ChEBI)]
  
  SMILES = c()
  for(i in 1:length(R.participant)) {
    SMILES = c(SMILES, parsed_ChEBI[parsed_ChEBI$ChEBI == R.participant[i], 'SMILES'])
  }
  ind = grep('\\*', SMILES)
  ind = c(ind, grep('^$', SMILES))
  R = R.participant[ind]
  print('making hash on compound class')
  h.class = Rhea.buildHash.class(parsed_ChEBI, R) # gain compound class-instances hash from parsed_ChEBI
  
  print('making hash on molecular formula')
  
  # make hash (ChEBI ID - molecular formula)
  participant = c(R.participant, unlist(values(h.class)))
  participant = unique(participant)
  h.formula = hash()
  for(i in 1:length(participant)) {
    if(grepl('^$', parsed_ChEBI[parsed_ChEBI$ChEBI == participant[i],'InChI']) == FALSE) {
      regexp.p = '(.*)(/p[-+]*[0-9]*)(.*)'
      regexp.h = '(.*)(H[0-9]*)(.*)'
      if(grepl(regexp.p, parsed_ChEBI[parsed_ChEBI$ChEBI == participant[i],'InChI'])) {
        tmp = sub(regexp.p, '\\2', parsed_ChEBI[parsed_ChEBI$ChEBI == participant[i],'InChI'])
        tmp2 = unlist(strsplit(tmp, '/p'))[2]
        tmp3 = as.numeric(tmp2)
        form = unlist(strsplit(parsed_ChEBI[parsed_ChEBI$ChEBI == participant[i],'InChI'], '/'))[2]
        form.h = sub(regexp.h, '\\2', form)
        if(form.h == 'H') {
          form.h = 'H1'
        }
        form.h2 = unlist(strsplit(form.h, 'H'))[2]
        if(length(form.h2) > 0) {
          n = as.numeric(form.h2)
        } else {
          n = 1
        }
        n2 = n + tmp3
        tmp = paste(sub(regexp.h, '\\1', form), 'H', n2, sub(regexp.h, '\\3', form), sep = '')
        tmp = gsub('H0', '', tmp)
      } else {
        tmp = unlist(strsplit(parsed_ChEBI[parsed_ChEBI$ChEBI == participant[i],'InChI'], '/'))[2]
      }
      .set(h.formula, keys = parsed_ChEBI[parsed_ChEBI$ChEBI == participant[i],'ChEBI'], values = tmp)
    }
  }
  # exception
  h.formula[['15378']] = 'H'
  
  # Instantiate
  print('Instantiating')
  
  registerDoMC(multicore)
  tmp = foreach(i = 1:length(Rhea_ID)) %dopar% {
    #for(i in 1:length(Rhea_ID)) {
    instantiate(Rhea_ID[i], parsed_Rhea[parsed_Rhea$ID == Rhea_ID[i],'Eq.ChEBI'], h.class, h.formula, R)
  }
  ID = unlist(tmp)[which(grepl('ID', names(unlist(tmp))))]
  Eq.ChEBI = unlist(tmp)[which(grepl('Equation', names(unlist(tmp))))]
  tmp2 = cbind(ID, Eq.ChEBI)
  rownames(tmp2) = 1:nrow(tmp2)
  tmp3 = as.data.frame(tmp2, stringsAsFactors=FALSE)
  return(tmp3)
}
