# (outgroups,(allOhterQuercus,((Ponticae), ((Virentes), ((Dumosae), (Stellatae + Mexicanae, (Prinoideae, ((Albae), (Roburoids)))))))));

# 0. Create list of tip vectors
taxSetList <- list(
  Ponticae = radMetaOut$code[which(radMetaOut$section == 'Ponticae' & radMetaOut$singleTip)],
  Virentes = radMetaOut$code[which(radMetaOut$section == 'Virentes' & radMetaOut$singleTip)],
  Dumosae = radMetaOut$code[which(radMetaOut$clade == 'Dumosae' & radMetaOut$singleTip)],
  Stellatae.Leucomexicana = radMetaOut$code[which(radMetaOut$clade %in% c('Stellatae', 'Leucomexicana', 'Tx white oaks') & radMetaOut$singleTip)],
  Prinoids = radMetaOut$code[which(radMetaOut$clade == 'Prinoids' & radMetaOut$singleTip)],
  Albae = radMetaOut$code[which(radMetaOut$clade == 'Albae' & radMetaOut$singleTip)],
  Roburoids = radMetaOut$code[which(radMetaOut$clade == 'Roburoids' & radMetaOut$singleTip)]
)
taxSetList$allOtherQuercus <- setdiff(radMetaOut$code[intersect(grep('Quercus', radMetaOut$tip), which(radMetaOut$singleTip))], unlist(taxSetList))
taxSetList$outgroups <- radMetaOut$code[intersect(grep('Quercus', radMetaOut$tip, invert = T), which(radMetaOut$singleTip))]

# 1. Fill them with commas...
taxSet <- sapply(taxSetList, paste, collapse = ',')

# 2. ... then concatenate them with parentheses...
taxConstraint <- paste(
  '(', taxSet['outgroups'],
  ',(', taxSet['allOtherQuercus'],
  ',((', taxSet['Ponticae'], ')',
  ',((', taxSet['Virentes'], ')',
  ',((', taxSet['Dumosae'], ')',
  ',(', taxSet['Stellatae.Leucomexicana'],
  ',(', taxSet['Prinoids'],
  ',((', taxSet['Albae'], ')',
  ',(', taxSet['Roburoids'],
  ')))))))));',
  sep = ''
)
tr.constraint <- read.tree(text = taxConstraint)
tr.constraint.labeled <- tr.constraint
tr.constraint.labeled$tip.label <- radMetaOut[tr.constraint$tip.label, 'tip']

# e. ... and write to the hard-drive!
write.tree(tr.constraint, '../OUT/ANALYSIS.PRODUCTS/taxaConstraint.tre')
write.tree(tr.constraint.labeled, '../OUT/ANALYSIS.PRODUCTS/taxaConstraintLabeled.tre')
