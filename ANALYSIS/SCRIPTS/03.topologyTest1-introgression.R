## set up topology test for locus-by-locus analyses
## 2019-01-31

library(magrittr) # b/c you always need magrittr!
source('../SCRIPTS/99.makeConstraint.R')
consDir <- format(Sys.time(), "locConstraints_%Y-%m-%d")
raxDir <- format(Sys.time(), "locRaxml_%Y-%m-%d")
if(!consDir %in% dir()) dir.create(consDir)
if(!raxDir %in% dir()) dir.create(raxDir)
## here, we set up trees to address the question of which loci are introgressed
## with respect to dumosae and albae

## constraints:
##  albae.mono = clades constrained if albae+roburoids is monophyletic
##  albae.intr = clades constrained if roburoids are introgressed with ponticae
##  macro.mono = clades constrained if macrocarpae excludes dumosae
##  macro.intr = clades constrained if macrocarpa and dumosae are intertwined
## for testing, ask whether lnL for the monophyletic constraint > lnL for the introgressed contstraint(s)

## n.b. simplistic constraints, but they give us a way of classing loci
constraints <- list(
  albae.mono = list(albae = unlist(taxa[c('albae', 'roburoid')]),
                    ponticae = taxa$ponticae),
  albae.intr = list(unlist(taxa[c('pontica', 'roburoid')])),
  dumos.mono = list(dumosae = taxa$dumosae,
                    macrocarpae = taxa$macrocarpae),
  dumos.int1 = list(unlist(taxa[c('dumosae', 'macrocarpa')])),
  dumos.int2 = list(unlist(taxa[c('macrocarpae', 'lobata')]))
)

constraintTreesFull <- sapply(constrraxDiraints, makeConstraint,
                              taxa = row.names(rads.subset$indsMatFull))
sapply(names(constraintTreesFull),
  function(x) writeLines(constraintTreesFull[x], paste(x, 'tre', sep = '.'))
)

consBatchOut <- '#!/bin/sh'
for(i in lociToWrite) {
  inds <- which(rads.subset$indsMatFull[, i]) %>% names
  temp <- sapply(constraints, makeConstraint,
                  taxa = names(which(rads.subset$indsMatFull[, i])))
  for(j in names(constraints)) {
    writeLines(temp[[j]], paste(consDir, '/', i, '_', j, '.tre', sep = ''))
    consBatchOut <-
      c(consBatchOut,
        paste('raxmlHPC-PTHREADS-AVX -T 4 -f d -m GTRCAT -p 12345 -s ',
        paste('../', locDir, '/', i, '.phy', sep = ''),
        ' -g ', paste('../', consDir, '/', i, '_', j, '.tre', sep = ''),
        ' -n ', paste(i, '_', j, '.tre', sep = ''),
        sep = ''))
  } # close j
} # close i

consBatchOut <- c(consBatchOut, 'echo finished with phylogenies, cleaning up')
consBatchOut <- c(consBatchOut, paste('\nrm ', '../', locDir, '/*.phy.reduced', sep = ''))
consBatchOut <- c(consBatchOut, paste('rm RAxML_info*', sep = ''))
consBatchOut <- c(consBatchOut, paste('rm RAxML_log*', sep = ''))
consBatchOut <- c(consBatchOut, paste('rm RAxML_parsimony*', sep = ''))
consBatchOut <- c(consBatchOut, paste('rm RAxML_result*', sep = ''))

writeLines(consBatchOut, paste(raxDir, 'AAA.locConstraints.rax.sh', sep = '/'))
