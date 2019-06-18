library(ape)
library(RADami)
library(magrittr)
library(phytools)
source('../SCRIPTS/99.split.plotTree.R')

tr.big <- read.tree('../TREES/big.tree_2018-07/tr.trial.v7-gtr.tre.export-withDeletions_noQuotes.tre')
tr.big$tip.label <- sapply(strsplit(tr.big$tip.label, '_|_', fixed = T), tail, n=1) %>%
  tidyName(case = 'upper')
tr.big2 <- drop.tip(tr.big, setdiff(tr.big$tip.label, radMetaOut$code))
tr.big2$tip.label <- radMetaOut[tr.big2$tip.label, 'tip']
write.tree(tr.big2, '../OUT/SUPPLEMENTS/DATA.S.full.tips.tre')
pdf('../OUT/SUPPLEMENTS/FIG.S.fullTipsTree.pdf', 10, 50)
plot(tr.big2, show.node.label = TRUE, cex = 0.3, edge.width = 0.5)
dev.off()

ppN = 7
split.plotTree(tr.big2,file = '../OUT/SUPPLEMENTS/FIG.S.fullTipsTree-split.pdf',
              splits = seq(from = 1/ppN, to = 1-1/ppN, by = 1/ppN),
              fsize = 0.5, lwd = 0.5)
