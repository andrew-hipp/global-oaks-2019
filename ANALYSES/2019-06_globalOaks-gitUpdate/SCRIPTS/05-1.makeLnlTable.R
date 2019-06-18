## take all the lnL log files and stuff the lnL into a locus-by-constraint table
library(tidyr)

## 1. read in log-likelihoods using 2nd element of read.table, create named vector
logFiles <- dir(raxLnLdir, patt = 'log')
lnL_vect <- sapply(logFiles, function(x) read.table(paste(raxLnLdir, '/', x, sep = ''))[[2]])

## 2. extract locus# and constraint from each file name, set as columns of data.frame
temp <- strsplit(names(lnL_vect), '_|log.') %>% simplify2array
lnL_mat <- data.frame(
  lnL = as.numeric(lnL_vect),
  locus = temp[3, ],
  constraint = temp[4, ]
)

## 3. populate table
lnL_mat <- spread(lnL_mat, constraint, lnL)
row.names(lnL_mat) <- lnL_mat$locus
lnL_mat$locus <- NULL

## 4. likelihood ratios of each topology
lnL_mat$albae.ratio <- lnL_mat$albae.mono - lnL_mat$albae.intr
lnL_mat$dumos.ratio1 <- lnL_mat$dumos.mono - lnL_mat$dumos.int1
lnL_mat$dumos.ratio2 <- lnL_mat$dumos.mono - lnL_mat$dumos.int2
lnL_mat$dumos.ratioC <- apply(lnL_mat[c('dumos.ratio1', 'dumos.ratio2')], 1, function(x) x[which(abs(x) == max(abs(x)))][1])
lnL_mat$monoSum <- lnL_mat$albae.ratio + lnL_mat$dumos.ratioC

## 5. write table to file
write.csv(lnL_mat, 'lnL_matrix.csv')
