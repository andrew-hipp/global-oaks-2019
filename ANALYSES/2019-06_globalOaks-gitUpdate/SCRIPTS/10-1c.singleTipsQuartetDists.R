## quartet dists betw locus trees and pruned global tree
library(Quartet)
library(magrittr)
steelThresh = 0.01 # similarities below this level get dropped out

tr.singletons.rawTips <- read.tree(dir(path, patt = 'bipartitions.2', full = T)) %>%
  ladderize

# 1. read in trees and squash down the branches that are < 1E-05
if(!exists('tr.singletonsByLoc')) {
  tr.singletonsByLoc <- lapply(dir(allLocDir, patt = 'bestTree', full = T), function(x) {
    read.tree(x) %>% di2multi(tol = 1E-05) })
  names(tr.singletonsByLoc) <- dir(allLocDir, patt = 'bestTree', full = F) %>%
    gsub(pattern = 'RAxML_bestTree.', replacement = '', fixed = TRUE)
  }

# 2. eliminate trees with < 2 nodes
tr.singletonsByLoc.dropped <- sapply(tr.singletonsByLoc, '[[', 'Nnode') %>% '==' (1) %>% which
tr.singletonsByLoc.dropped <- tr.singletonsByLoc[tr.singletonsByLoc.dropped]
keepsies <- sapply(tr.singletonsByLoc, '[[', 'Nnode') %>% '>'(1) %>% which
tr.singletonsByLoc <- tr.singletonsByLoc[keepsies]
rm(keepsies)

## ... and get stats:
sapply(tr.singletonsByLoc, '[[', 'Nnode') %>% sort(decreasing = T) %>% mean
# [1] 4.476828
sapply(tr.singletonsByLoc, '[[', 'Nnode') %>% sort(decreasing = T) %>% sd
# [1] 1.825909
sapply(tr.singletonsByLoc, '[[', 'Nnode') %>% sort(decreasing = T) %>% median
# [1] 4
sapply(tr.singletonsByLoc, '[[', 'Nnode') %>% sort(decreasing = T) %>% range
# [1]  2 15


# 3. write squashed trees back to disc for phyparts, which I want NOT to use
#    spurious bipartitions.
if(!'squashedTrees' %in% dir(allLocDir)) {
  dir.create(paste(allLocDir, 'squashedTrees/', sep = '/'))
  lapply(names(tr.singletonsByLoc),
          function(i) {
            write.tree(tr.singletonsByLoc[[i]],
                       paste(allLocDir, '/squashedTrees/RAxML_bestTree.squashed.',
                              i, '.tre', sep = ''))
                     } # close fct(i)
                   ) # close lapply
                 } # close if

# 3. get sharedQuartetStatus and define similarity as # of nodes resolved the same
#    over the total number of nodes resolved in both trees --
#    ignores tree size and bipartitions in one tree but not the other
if(!exists('tr.singletonsByLoc.q')) {
  tr.singletonsByLoc.q <- SharedQuartetStatus(tr.singletonsByLoc, tr.singletons.rawTips)
  tr.singletonsByLoc.s <- tr.singletonsByLoc.q[, 's'] /
                          (tr.singletonsByLoc.q[, 's'] + tr.singletonsByLoc.q[, 'd'])
 }

message(paste('undefined distances:', sum(is.na(tr.singletonsByLoc.s))))
