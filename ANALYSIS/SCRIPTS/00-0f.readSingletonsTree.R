library(ape)
library(magrittr)

path <- '../TREES/raxml.2019-03-13/'

# 1. read in singletons tree
tr.singletons <- read.tree(dir(path, patt = 'bipartitions.2', full = T)) %>%
  ladderize

# 2. relabel
tr.singletons$tip.label <- radMetaOut[tr.singletons$tip.label, 'tip']
pdf('../OUT/ANALYSIS.PRODUCTS/tr.singletons.pdf', 10, 30)
plot(tr.singletons, cex = 0.5, show.node.label = T)
dev.off()
