## checking whether the loci favoring more introgression really show an
##   introgression tree and vice versa

library(RADami)
library(magrittr)

indsThresh = 4
locThresh = 20
dir.out <- 'introgressionCheck.2019-02-19'

if(!dir.out %in% dir()) dir.create(dir.out)

lociList <- list(
  albaeMono = row.names(lnL_mat)[which(lnL_mat$albae.ratio > 0)],
  albaeIntr = row.names(lnL_mat)[which(lnL_mat$albae.ratio < 0)],
  dumosMono = row.names(lnL_mat)[which(lnL_mat$dumos.ratioC > 0)],
  dumosIntr = row.names(lnL_mat)[which(lnL_mat$dumos.ratioC > 0)]
)

for(i in names(lociList)) {
  print(paste('doing', i))
  rad2phy(rads.subset$matFull,
          loci = lociList[[i]],
          inds = apply(rads.subset$indsMatFull[, lociList[[i]]], 1, sum) %>%
                '>='(locThresh) %>%
                which %>%
                names,
          padding = 100,
          outfile = paste(dir.out, '/', i, '.phy', sep = ''),
          logfile = 'temp.log')
}
