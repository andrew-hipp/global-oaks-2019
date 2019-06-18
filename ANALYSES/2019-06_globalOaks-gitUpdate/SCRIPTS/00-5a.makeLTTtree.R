source('../SCRIPTS/99.add.random.R')
library(ape)
library(openxlsx)

dat.taxTable <- read.xlsx('../DATA/cladeNumbers-v2.xlsx', 1, rowNames = T)
lambdaVect = c('1', '10')
calibVect = c('taxaGrepCrown', 'taxaGrepStem')
trList.LTT <- list(lambda = vector('list', length(lambdaVect)),
                   calibPoint = vector('list', length(calibVect))
                 ) # close list
names(trList.LTT[['lambda']]) <- lambdaVect
names(trList.LTT[['calibPoint']]) <- calibVect

for(lambda in lambdaVect) {
  for(calibPoint in calibVect) {
    trTemp <- tr.singletons.correlated[[lambda]][[calibPoint]]
    class(trTemp) <- 'phylo'
    for(i in which(dat.taxTable$samplingFraction < 1)) {
      trTemp <- add.random.startNode(trTemp,
                                     tips = paste(dat.taxTable$shortName[i],
                                                  seq(dat.taxTable$missing[i]),
                                                  sep = '.'
                                                  ),
                                     startNodeGrep = dat.taxTable$grepInclude[i])
      } # close i
    tr.singletons.correlated.pruneSet <- list(
      sect_Lobatae = extract.clade(trTemp,
        findMRCA(trTemp, grep('Quercus_agrifolia|Quercus_emoryi', trTemp$tip.label))),
      sect_Quercus = extract.clade(trTemp,
          findMRCA(trTemp, grep('Quercus_lobata|Quercus_arizonica', trTemp$tip.label))),
      sect_Ilex = extract.clade(trTemp,
        findMRCA(trTemp, grep('Quercus_franchetii|Quercus_rehderiana', trTemp$tip.label))),
      sect_Cerris = extract.clade(trTemp,
        findMRCA(trTemp, grep('Quercus_variabilis|Quercus_cerris', trTemp$tip.label))),
      sect_Cyclobalanopsis = extract.clade(trTemp,
        findMRCA(trTemp, grep('Quercus_gilva|Quercus_phanera', trTemp$tip.label)))
      )

    class(tr.singletons.correlated.pruneSet) <- 'multiPhylo'
    trList.LTT[[lambda]][[calibPoint]] <- tr.singletons.correlated.pruneSet
    write.tree(trTemp, paste('../OUT/ANALYSIS.PRODUCTS/lttTree', lambda, calibPoint, 'tre', sep = '.'))

    pdf(paste('../OUT/ANALYSIS.PRODUCTS/lttTree', lambda, calibPoint, 'pdf', sep = '.'), 12, 30)
    plot(trTemp, cex = 0.5)
    dev.off()
  } # close calibPoint
}# close lambda
