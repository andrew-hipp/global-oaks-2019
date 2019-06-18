## compare log likelihood for all trees for each locus

raxLnLdir <- paste(raxDir, 'lnL', sep = "_")
if(!raxLnLdir %in% dir()) dir.create(raxLnLdir)

files.rax <- dir(raxDir, patt = 'tre')
files.loc <- paste(sapply(strsplit(files.rax, 'bestTree.|_'), '[', 3), 'phy', sep = '.')

out <- c('#!/bin/sh',
  paste('raxmlHPC-PTHREADS-AVX -f e -t ',
        '../', raxDir, '/', files.rax,
        ' -m GTRGAMMA -s ../', locDir, '/', files.loc,
        ' -n ', sapply(strsplit(files.rax, 'bestTree.|.tre'), '[', 2), '_lnL.txt',
        sep = '')
      )

out <- c(out, 'echo finished with phylogenies, cleaning up')
out <- c(out, paste('\nrm ', '../', locDir, '/*.phy.reduced', sep = ''))
out <- c(out, paste('rm RAxML_info*', sep = ''))
out <- c(out, paste('rm RAxML_binaryModelParameters*', sep = ''))
out <- c(out, paste('rm RAxML_result*', sep = ''))

writeLines(out, paste(raxLnLdir, '/AAA.getLnL.sh', sep =''))
