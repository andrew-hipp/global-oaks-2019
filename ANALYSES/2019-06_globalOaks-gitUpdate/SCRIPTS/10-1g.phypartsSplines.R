library(parallel)
library(ncf)

chromsToDo <- sort(unique(phyParts.df$chrom))
phyParts.splines <- mclapply(chromsToDo, function(i) {
  print(paste('doing splines for', i))
  temp = phyParts.df[phyParts.df$chrom == i, ]
  tempG <- dist(temp$genomeLinStart)
  temp$iBin <- ifelse(temp$info == 'conflict', 0, 1)
  tempI <- dist(temp$iBin)
  temp.spline <- spline.correlog(x = temp$cleanStart, y = rep(1, dim(temp)[1]), z = temp$iBin)
  return(temp.spline)
  }, mc.cores = length(chromsToDo)
)
names(phyParts.splines) <- chromsToDo

pdf('../OUT/FIGS.AND.TABLES/FIGURE.phyparts.splines.pdf', 11, 8.5)
layout(matrix(1:12, 3, 4, byrow = T))
for(i in chromsToDo) {
  plot(phyParts.splines[[i]], main = i)
}
dev.off()
