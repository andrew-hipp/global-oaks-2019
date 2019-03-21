raxDirUnconstrained <- gsub('axml', 'axml-unconstrained', raxDir)
if(!raxDirUnconstrained %in% dir()) {
  dir.create(raxDirUnconstrained)
    i = dir(locDir)
    out <- '#!/bin/sh'
    out <- c(out, paste('raxmlHPC-PTHREADS-AVX -f d -m GTRCAT -p 12345 -s ',
                        paste('..', locDir, i, sep = '/'),
                        ' -n ', gsub('.phy', '.tre', i, fixed = T),
                        sep = ''))
    out <- c(out, 'echo finished with phylogenies, cleaning up')
    out <- c(out, paste('\nrm ', '../', locDir, '/*.phy.reduced', sep = ''))
    out <- c(out, paste('rm RAxML_info*', sep = ''))
    out <- c(out, paste('rm RAxML_log*', sep = ''))
    out <- c(out, paste('rm RAxML_parsimony*', sep = ''))
    out <- c(out, paste('rm RAxML_result*', sep = ''))
    writeLines(out, paste(raxDirUnconstrained, 'AAA.loc.rax.unconstrained.sh', sep = '/'))
}
