## funny story: I exported ALL the loci... only needed to export the ones we
##   mapped. Fixing that now. In next round, only send the ones we need.

tempFiles <- dir('allLocRax', patt = 'bestTree') %>%
  gsub(pattern = 'RAxML_bestTree.', replacement = '', fixed = TRUE)

removeFiles <- setdiff(tempFiles, oaks.map$uniques$name)
writeLines(paste('rm RAxML_bestTree.', removeFiles, sep = ''), 'allLocRax/removeEm.sh')
