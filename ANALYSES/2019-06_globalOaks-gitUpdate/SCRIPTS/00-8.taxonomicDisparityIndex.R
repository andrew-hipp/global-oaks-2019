## taxonomic disparity for big tree
library(ape)
if(!exists('summary.by.elements')) {
  source('https://raw.githubusercontent.com/andrew-hipp/morton/master/R/label.elements.R')
  source('https://raw.githubusercontent.com/andrew-hipp/morton/master/R/summary.by.elements.R')
  source('https://raw.githubusercontent.com/andrew-hipp/morton/master/R/tips.expected.R')
  }

tr.big2.summary <-
  summary.by.elements(drop.tip(tr.big2, grep('Quercus', tr.big2$tip.label, invert = T)),
                      delim= '[ _|]', fixed = F, returnNum=1:2)
  tr.big2.summary$disparity.mat <-
    tr.big2.summary$disparity.mat[grep('Quercus sp.|Quercus new', row.names(tr.big2.summary$disparity.mat), invert = T), ]
write.csv(tr.big2.summary$disparity.mat,
          '../OUT/SUPPLEMENTS/TABLE.S.tr.big2.TDI.matrix.csv')
tr.big2.summary$TDI.clean <- as.data.frame(
  tr.big2.summary$disparity.mat
     )
tr.big2.summary$TDI.clean <- tr.big2.summary$TDI.clean[
  grep('|', row.names(tr.big2.summary$TDI.clean), fixed = T, invert = T),
]
tr.big2.summary$stats <- c(
  'All-tips tree stats',
  '-------------------',
  paste('Total Quercus spp:', dim(tr.big2.summary$TDI.clean)[1]),
  paste('Total Quercus with 1 sample:', sum(tr.big2.summary$TDI.clean$count == 1)),
  paste('Total Quercus with >1 sample:', sum(tr.big2.summary$TDI.clean$count > 1)),
  paste('Mean number of samples for Quercus with > 1 sample:',
    round(mean(tr.big2.summary$TDI.clean$count[which(tr.big2.summary$TDI.clean$count > 1)]), 2),
    '+/-',
    round(sd(tr.big2.summary$TDI.clean$count[which(tr.big2.summary$TDI.clean$count > 1)]), 2),
    '(sd)'),
  paste('Total Quercus with TDI > 0:', sum(tr.big2.summary$TDI.clean$disparity > 0)),
  paste('Total Quercus with TDI 10 or more:', sum(tr.big2.summary$TDI.clean$disparity >= 10))
)

writeLines(tr.big2.summary$stats, '../OUT/ANALYSIS.PRODUCTS/treeSummaryStats.tr.big2.txt')
