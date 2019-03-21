library(ape)
library(Quartet)
library(magrittr)
source('../SCRIPTS/99.taxaTest.R')

if(!exists('tr')) tr <- read.tree('../TREES/raxml.2018-01-23/RAxML_bipartitions.RAxML_result.oaksall_v1_2.m15.2018-01-23.singles-v3_biparts.tre') %>% ladderize

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
#taxa <- sapply(taxa, grep, x = tr$tip.label, value = T)
taxa <- sapply(taxa, grep, x = row.names(rads.subset$indsMatFull), value = T)

#if(all(sapply(taxa, is.monophyletic, phy = tr))) message('Taxon sets are all monophyletic!') else
#message('Some taxon set is not monophyletic... \n***BEWARE!!!!***')

if(any(sapply(taxa, length) == 1)) {
  taxaSums <- cbind(
    sapply(names(taxa)[sapply(taxa, length) > 1], function(x) colSums(rads.subset$indsMatFull[taxa[[x]], ])),
    sapply(names(taxa)[sapply(taxa, length) == 1], function(x) as.integer(rads.subset$indsMatFull[taxa[[x]], ]))
  )
} else taxaSums <- sapply(names(taxa), function(x) colSums(rads.subset$indsMatFull[taxa[[x]], ]))

taxaSums <- cbind(taxaSums,
  dumTest = apply(taxaSums[, c('lobata', 'dumosae', 'macrocarpa', 'macrocarpae')], 1,
                  function(x) all(x >= c(1, 2, 1, 2))),
  robTest = apply(taxaSums[, c('roburoid', 'albae', 'pontica')], 1,
                  function(x) all(x >= 1))
  )
taxaSums <- cbind(taxaSums, rob.dumTest = and(taxaSums[, 'dumTest'], taxaSums[, 'robTest']))
