library(ggplot2)
library(dplyr)
library(ncf)

oaks.map$quartetD <- cbind(oaks.map$uniques[names(tr.singletonsByLoc.s), ],
                          QuartetSimilarity = tr.singletonsByLoc.s,
                          LocusTreeLength = sapply(tr.singletonsByLoc, function(x) length(x$tip.label)))
oaks.map$quartetD$LG <- gsub('chr', '', oaks.map$quartetD$chrom) %>%
                        as.integer
oaks.map$quartetD <- oaks.map$quartetD[!is.na(oaks.map$quartetD$QuartetSimilarity), ]
oaks.map$quartetD <- oaks.map$quartetD[order(oaks.map$quartetD$chrom, oaks.map$quartetD$cleanStart), ]

pdf('../OUT/FIGS.AND.TABLES/FIGURE.chromoRadMap_quartetDistances-v2.pdf', 11, 8.5)
p <- ggplot(oaks.map$quartetD, aes(x = chrom, y = cleanStart,
                                label = name))
#p <- p + scale_color_gradient2(
#                    "Monophyly\nlnL difference",
#                    midpoint=0, low="blue", mid="white",
#                    high="red", space ="Lab" )

p <- p + scale_color_gradient(
                    "RAD locus suport\nfor consensus tree",
                    low="blue", high="yellow", space ="Lab" )

p <- p + geom_point(data = oaks.map$uniques,
                    aes(x = chrom, y = cleanStart),
                    color = 'gray50', alpha = 0.1,
                    pch = '-', size = 9,
                    na.rm = T)
p <- p + geom_point(data=oaks.map$quartetD,
                    aes(x = LG + 0.2, y = (cleanStart), color = QuartetSimilarity),
                    pch = '-', size = 7,
                    na.rm = T)

p <- p + theme_tufte()
p <- p + theme(legend.position = c(0.9, 0.8))
print(p)
dev.off()


pdf('../OUT/FIGS.AND.TABLES/FIGURE.splinesByChromosome_Quartets-v1.pdf', 11, 8.5)
layout(matrix(1:12, 3, 4))
oaks.map$quartetSplines <- vector('list', 12)
names(oaks.map$quartetSplines) <- unique(oaks.map$quartetD$chrom)
for(i in unique(oaks.map$quartetD$chrom)) {
  print(paste('doing', i))
  oaks.map$quartetSplines[[i]] <- spline.correlog(x = oaks.map$quartetD$cleanStart[oaks.map$quartetD$chrom == i], # dist along chromosome
                          y = rep(1, sum(oaks.map$quartetD$chrom == i)), # fill y with 1s
                          z = oaks.map$quartetD$QuartetSimilarity[oaks.map$quartetD$chrom == i]) # which chromosome
  plot(oaks.map$quartetSplines[[i]], main = i)
}
dev.off()
