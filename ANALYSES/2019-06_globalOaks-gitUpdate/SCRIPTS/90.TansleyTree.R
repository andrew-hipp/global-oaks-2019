## tree for Tansley review, Kremer and Hipp 2019, New Phytologist
library(ggtree)
library(phytools)

cbbPalette <- c("#EEEEE9", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

tr.nodes <- list(
  Subgenus = sapply(unique(tr.quickLabeler$subgenus), function (x) {
    getMRCA(tr.calibrated, tr.quickLabeler$tip[tr.quickLabeler$subgenus == x])
  }),
  Section = sapply(unique(tr.quickLabeler$section), function (x) {
    getMRCA(tr.calibrated, tr.quickLabeler$tip[tr.quickLabeler$section == x])
  })
)

trg <- groupClade(tr.calibrated, tr.nodes$Section, group_name = 'Section')
#names(attributes(trg))[names(attributes(trg)) == 'group'] <- "Section"
levels(attr(trg, 'Section'))[levels(attr(trg, 'Section')) == '0'] <- 'Basal lineages'
p <- ggtree(trg, aes(color=Section), layout = 'circular')
p <- p + scale_color_manual(values = cbbPalette)
p <- p + theme(legend.position="right")
p <- p + geom_treescale()

pdf('tansleyTree.v1.pdf')
print(p)
dev.off()
