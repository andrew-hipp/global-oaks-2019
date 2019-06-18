## write out loci and generate tree for each
## 2019-02-19: revised to just be for the single-tips tree, so we can map these back
##   only loci with at least 15 individuals in tree
##   note: we need rooted tree! so make sure there is a non-quercus in each

library(ape)
library(RADami)
library(magrittr)
library(parallel)

nIndsMin <- 10
ncores <- 4
writeLoci = TRUE

#the following line substitutes corrected tip names for extraction codes
rads.subset$tips <- radMetaOut[row.names(rads.subset$indsMat), 'tip']
rads.subset$tipsToNames <-
  setNames(row.names(rads.subset$indsMat), rads.subset$tips)

allLocDir <- format(Sys.time(), 'allLociPhy_%Y-%m-%d')
tr.write.loci <-
  (colSums(rads.subset$indsMat) >= nIndsMin) %>%
  which %>%
  names %>%
  intersect(oaks.map$uniques$name)

##HERE I NEED TO SUBSET BY LOCI THAT HAVE AN OUTGROUP
tr.write.loci <- apply(rads.subset$indsMat[grep('Quercus', rads.subset$tips, invert = T), tr.write.loci], 2, any) %>%
  which %>%
  names

tr.write.og <- apply(rads.subset$indsMat[ , tr.write.loci], 2, function(x) {
  rads.subset$tipsToNames[
    grep("Quercus", rads.subset$tips[which(x)], invert = T, value = T)[1]
    ]
  })
names(tr.write.og) <- tr.write.loci

# 1. write loci to a folder
if(!allLocDir %in% dir()) dir.create(allLocDir)
out <- character(0)

  mclapply(tr.write.loci, function(i) {
    message(paste('doing', i))
    indsTemp <- rads.subset$tipsToNames[rads.subset$tips[which(rads.subset$indsMat[, i])]] %>%
      as.character
    ogTemp <- rads.subset$tipsToNames[
      grep('Quercus', indsTemp, invert = T, value = T)[1]
      ]
    rad2phy(rads.subset$mat,
      loci = i, padding = 100,
      inds = indsTemp,
      outfile = paste(allLocDir, '/', i, '.phy', sep = ''),
      logfile = 'temp.discard.log'
      )
     }, mc.cores = ncores
    )

out <- paste('raxmlHPC-PTHREADS-AVX -T 2 -f d -m GTRCAT -p 12345 -s ',
         '../', allLocDir, '/', tr.write.loci, '.phy',
         " -o '", tr.write.og, "'",
         ' -n ', tr.write.loci,
         ' &',
          sep = ''
          )

# 2. write a batchfile to generate a tree for each and cleanup

out <- rbind(matrix(out, nrow = ncores), 'rm RAxML_info* RAxML_log* RAxML_parsimonyTree* RAxML_result* &', 'wait')
out <- as.character(out)
out <- c('!#/bin/sh', out)
writeLines(out, paste(allLocDir, '/aaa.runRax.sh', sep = ''))
