library(RADami)
library(openxlsx)
library(magrittr)

#source('https://raw.githubusercontent.com/andrew-hipp/RADami/master/R/read.pyRAD.R')

minInds <- 15 # the minimum number of taxa per loci

# 1. read RAD data, maintain original tips
if(!exists('rads')) {
  rads <- read.pyRAD('../../../RAD.DAT/oaksall_v1_2.m15.loci')
  rads.mat <- rad2mat(rads)
  rads.indsMat <- rads$radSummary$inds.mat
  a <- row.names(rads.indsMat)
  a.fix <- grep('.', a, fixed = T)
  a[a.fix] <-
    sapply(strsplit(a[a.fix], '.', fixed = T), '[', i=1)
  a <- gsub('_2012|_2010', '', a)
  row.names(rads.indsMat) <-
    row.names(rads.mat) <-
    a
  rm(a, a.fix)
}

# clean up the codes now, not above, as this ties it up
# with the costly step of reading the RADs
row.names(rads.mat) <- tidyName(row.names(rads.mat), case = 'upper')
row.names(rads.indsMat) <- tidyName(row.names(rads.indsMat), case = 'upper')

# 2. subset data by individuals
rads.subset <- list(
  indsMatFull = rads.indsMat[radMetaOut$code, ],
  matFull =     rads.mat[radMetaOut$code, ],
  indsMat =     rads.indsMat[radMetaOut$code[which(radMetaOut$singleTip)], ],
  mat =         rads.mat[radMetaOut$code[which(radMetaOut$singleTip)], ]
  )

# 3. subset data by loci, just down to loci present in >= 15 in the full matrix
rads.subset.loci <- names(which(colSums(rads.subset$indsMatFull) >= minInds))
rads.subset$indsMatFull <- rads.subset$indsMatFull[, rads.subset.loci]
rads.subset$indsMat <- rads.subset$indsMat[, rads.subset.loci]
rads.subset$matFull <- rads.subset$matFull[, rads.subset.loci]
rads.subset$mat <- rads.subset$mat[, rads.subset.loci]

class(rads.subset$mat) <- class(rads.subset$matFull) <- 'rad.mat'
