library(ape)
library(openxlsx)
library(phytools)
library(parallel)
library(magrittr)
library(phyloch) # if necessary, install from source: http://www.christophheibl.de/phyloch_1.5-3.tar.gz
data(strat2012) # A data frame containing the stratigraphic chart by issued in 2012 by the International Commission on Stratigraphy.

## ORIGINAL ANALYSES WAS CONDUCTED WITH :
#do.cvLambda <- TRUE
#lambdaVect <- c(0, 0.05, 0.1, 0.5, 1, 5, 10) # smoothing parameter for PL

## but for computational reasons on reanalysis, you can use this:
do.cvLambda <- FALSE # as you don't need to redo cross-validation
lambdaVect <- c(1,10) # as 1 was the lambda identified by cross-validation
                      # and 10 is used in the supplemental BAMM plot
calibOptions <- c('taxaGrepCrown', 'taxaGrepStem')

tr.singletons.pruned <- drop.tip(tr.singletons, grep('Quercus|Notholithocarpus', tr.singletons$tip.label, invert = T))
tr.calib <- read.xlsx('../DATA/calibrations-v3.xlsx', 1)
tr.calib.dataFrame <- list( taxaGrepCrown = NA,
                            taxaGrepStem = NA )
for(i in calibOptions) {
  tr.calib.dataFrame[[i]] <- data.frame(
    node = sapply(tr.calib[[i]][!is.na(tr.calib[[i]])], function(x) {
      findMRCA(tr.singletons.pruned, grep(x, tr.singletons.pruned$tip.label, value = T))
    }),
    age.min = tr.calib$minAge[!is.na(tr.calib[[i]])],
    age.max = tr.calib$maxAge[!is.na(tr.calib[[i]])],
    soft.bounds = NA
  )
  pdf(paste('../OUT/SUPPLEMENTS/FIG.S.singletons.', i, '.pdf', sep = ''), 15, 20)
  plot(tr.singletons.pruned, cex = 0.5)
  nodelabels(text = paste(tr.calib$nodeLabel[!is.na(tr.calib[[i]])], ' - ', tr.calib$minAge[!is.na(tr.calib[[i]])], ' Ma', sep = ''), node = tr.calib.dataFrame[[i]]$node)
  dev.off()
}

# this step doesn't vary by the calibrations
if(do.cvLambda) {
  tr.singletons.correlated.cv <- mclapply(lambdaVect, function(x) {
    chronopl(tr.singletons.pruned, x, CV = T)
    },
    mc.cores = length(lambdaVect)
  )
  lambdaVect.min <-
    sapply(tr.singletons.correlated.cv, function(x) attributes(x)$D2) %>%
    t %>%
    apply(., 1, sum) %>%
    abs %>%
    (function(x) which(x == min(x))) %>%
    lambdaVect[.]

  paste('Min CV: lambda =', lambdaVect.min) %>%
    writeLines(con = '../OUT/ANALYSIS.PRODUCTS/cvLambda.txt')
   # write minimum CV to file 'cvLambda.txt'
} # close do.cvLambda

# here, format so that the two calibration options (crown vs stem) are separated
tr.singletons.relaxed <- tr.singletons.correlated <-
  vector('list', length(lambdaVect))
names(tr.singletons.relaxed) <- as.character(lambdaVect)
for(i in lambdaVect) {
  print(i)
  tr.singletons.relaxed[[as.character(i)]] <- lapply(calibOptions, function(j) {
    chronos(tr.singletons.pruned, model = 'relaxed', calibration = tr.calib.dataFrame[[j]], lambda = i)
  })
  tr.singletons.correlated[[as.character(i)]] <- lapply(calibOptions, function(j) {
    chronos(tr.singletons.pruned, model = 'correlated', calibration = tr.calib.dataFrame[[j]], lambda = i)
  })
  names(tr.singletons.relaxed[[as.character(i)]]) <-
    names(tr.singletons.correlated[[as.character(i)]]) <-
    calibOptions
  } # close i

for(i in as.character(lambdaVect)) {
  for(j in calibOptions) {
    tr.singletons.relaxed[[i]][[j]] <-
      drop.tip(tr.singletons.relaxed[[i]][[j]],
               grep('Quercus', tr.singletons.relaxed[[i]][[j]]$tip.label, invert = T))
    tr.singletons.correlated[[i]][[j]] <-
      drop.tip(tr.singletons.correlated[[i]][[j]],
               grep('Quercus', tr.singletons.correlated[[i]][[j]]$tip.label, invert = T))

    write.tree(tr.singletons.relaxed[[i]][[j]], paste('../OUT/ANALYSIS.PRODUCTS/tr.singletons.relaxed', i, j, 'tre', sep = '.'))
    write.tree(tr.singletons.correlated[[i]][[j]], paste('../OUT/ANALYSIS.PRODUCTS/tr.singletons.correlated', i, j, 'tre', sep = '.'))

    pdf(paste('../OUT/ANALYSIS.PRODUCTS/tr.singletons.relaxed', i, j, 'pdf', sep = '.'), 15, 30)
    plot(tr.singletons.relaxed[[i]][[j]], cex = 0.7)
    add.geoscale(tr.singletons.relaxed[[i]][[j]])
    axisGeo(GTS = strat2012)
    dev.off()

    pdf(paste('../OUT/ANALYSIS.PRODUCTS/tr.singletons.correlated', i, j, 'pdf', sep = '.'), 15, 30)
    plot(tr.singletons.correlated[[i]][[j]], cex = 0.7)
    add.geoscale(tr.singletons.correlated[[i]][[j]])
    axisGeo(GTS = strat2012)
    dev.off()
  }
}

tr.singletons.relaxed.stats <-
  tr.singletons.correlated.stats <-
  list(taxaGrepCrown = NA, taxaGrepStem = NA )
for(j in calibOptions) {
  tr.singletons.relaxed.stats[[j]] <-
    sapply(tr.singletons.relaxed[as.character(lambdaVect)], function(x) attributes(x[[j]])$'PHIIC')
  tr.singletons.correlated.stats[[j]] <-
    sapply(tr.singletons.correlated[as.character(lambdaVect)], function(x) attributes(x[[j]])$'PHIIC')

  write.csv(tr.singletons.relaxed.stats[[j]], '../OUT/SUPPLEMENTS/TABLE.S.tr.singletons.relaxed.PhiIC.csv')
  write.csv(tr.singletons.correlated.stats[[j]], '../OUT/SUPPLEMENTS/TABLE.S.tr.singletons.correlated.PhiIC.csv')
}
