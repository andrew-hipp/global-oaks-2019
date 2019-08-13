## check on model-averaged diversification rates of clades for Tansley Review
## Kremer and Hipp, 2019, New Phytologist
## ahipp@mortonarb.org (2019-08-13)

library(BAMMtools)
library(overlapping)
library(magrittr)

pN <- function(x) paste(round(mean(x), 3), ' (',
    paste(round(quantile(x, c(0.025, 0.975)), 3), collapse = '-'),
    ')', sep = '')
tr.Tansley <- read.tree(
  '../OUT/ANALYSIS.PRODUCTS/tr.singletons.correlated.1.taxaGrepCrown.tre'
  )

if(!exists('e.Tansley')) {
  e.Tansley <- getEventData(tr.Tansley, '../BAMM.WORKING/crown-sampleProps/event_data.txt', burnin = 0.25)
  e.Tansley$tip.label <- sapply(strsplit(e.Tansley$tip.label, '|', fixed = T), function(x) {
    paste(x[1], x[length(x)], sep = '|')
    })
  e.Tansley$tip.label <- gsub('|', ' | ', e.Tansley$tip.label, fixed = TRUE)
  e.Tansley$tip.label <- gsub('_', ' ', e.Tansley$tip.label, fixed = TRUE)
  }

n.Tansley <- c(
  roburoid = 'mongolica|faginea',
  albae = 'alba|montana',
  prinoides = 'prinoides|bicolor',
  stellatae = 'sinuata|austrina',
  cerris = 'chenii|trojana',
  ilex = 'franchet|longispica'
)
n.Tansley <- lapply(n.Tansley, grep, x = tr.Tansley$tip.label, value = T)
n.Tansley <- sapply(n.Tansley, getMRCA, phy = tr.Tansley)

div.Tansley <- sapply(names(n.Tansley), function(x) {
  a = getCladeRates(e.Tansley, node = n.Tansley[x])
  a = a$lambda - a$mu
  a
})

pdf('../OUT/divTansely-2.pdf')
plot(density(div.Tansley[, 1], x.lim = range(div.Tansley)), col = cbbPalette[1])
for(i in 2:length(n.Tansley))
  lines(density(div.Tansley[, i]), col = cbbPalette[i])
dev.off()

apply(div.Tansley, 2, pN) %>% as.matrix %>% print

pdf('../OUT/all.overlap.pdf')
ov.Tansley <- overlap(as.list(as.data.frame(div.Tansley)), plot = TRUE)
dev.off()
ov.Tansley$OV %>% sort %>% print


pdf('../OUT/tansleyTreeCheck.pdf', 10, 20)
plot(tr.Tansley, cex = 0.6)
nodelabels(text = names(n.Tansley), node = n.Tansley)
dev.off()
