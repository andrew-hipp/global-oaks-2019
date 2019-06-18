library(ape)
library(Quartet)
library(magrittr)
source('../SCRIPTS/99.taxaTest.R')

tr <- read.tree('../TREES/raxml.2019-03-13/RAxML_bipartitions.2019-03-13.singles.v4.tre') %>% ladderize

taxa <- list(
  dumosae = 'lobata|garryana|berberidifolia|durata|cornelius|pacifica|dumosa|douglasii|tucker',
  lobata = 'lobata',
  pontica = 'pontica',
  ponticae = 'pontica|sadleriana',
  macrocarpae = 'macrocarpa|bicolor|lyrata',
  macrocarpa = 'macrocarpa',
  roburoid = 'serrata|mongolica|yunnanensis|Quercus_dentata|fabri|aliena|griffith|hartwissiana|robur|canariensis|pyrenaica|petraea|lusitanica|dalechampii|macranthera|frainetto|faginea|pubescens|vulcanica|pinnatifida|infectoria|cedrorum|boissieri',
  albae = 'alba|michauxi|montana'
)
#the following line substitutes corrected tip names for extraction codes
rads.subset$fullTips <- radMetaOut[row.names(rads.subset$indsMatFull), 'tip']
taxa <- sapply(taxa, grep, x = rads.subset$fullTips, value = T)

if(any(sapply(taxa, length) == 1)) {
  taxaSums <- cbind(
    sapply(names(taxa)[sapply(taxa, length) > 1], function(x) {
      colSums(rads.subset$indsMatFull[match(taxa[[x]], rads.subset$fullTips), ])
    }
  ),
    sapply(names(taxa)[sapply(taxa, length) == 1], function(x) {
      as.integer(rads.subset$indsMatFull[match(taxa[[x]], rads.subset$fullTips), ])
    }
  ) # close sapply
) # close cbind
} else taxaSums <- sapply(names(taxa), function(x) {
  colSums(rads.subset$indsMatFull[match(taxa[[x]], rads.subset$fullTips), ])
  })

taxaSums <- cbind(taxaSums,
  dumTest = apply(taxaSums[, c('lobata', 'dumosae', 'macrocarpa', 'macrocarpae')], 1,
                  function(x) all(x >= c(1, 2, 1, 2))),
  robTest = apply(taxaSums[, c('roburoid', 'albae', 'pontica')], 1,
                  function(x) all(x >= 1))
  )
taxaSums <- cbind(taxaSums, rob.dumTest = and(taxaSums[, 'dumTest'], taxaSums[, 'robTest']))
