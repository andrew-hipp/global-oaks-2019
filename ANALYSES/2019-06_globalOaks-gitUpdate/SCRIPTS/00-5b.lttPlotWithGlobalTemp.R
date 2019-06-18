library(ape)
library(openxlsx)
library(phyloch)
library(phytools)

data(strat2012) # A data frame containing the stratigraphic chart by issued in 2012 by the International Commission on Stratigraphy.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

zachos = read.xlsx('../DATA/zachos2001.xlsx')
zachos$Age = -zachos$Age
zachos$Temperature = 16.5 -
                      4.3*zachos$d18Oadj +
                      0.14*(zachos$d18Oadj^2)
for(lambda in c('1', '10')) {
  for(calibPoint in c('taxaGrepCrown', 'taxaGrepStem')) {
    tr.singletons.correlated.pruneSet <- trList.LTT[[lambda]][[calibPoint]]
    tr.singletons.cols <- structure(cbbPalette[1:length(tr.singletons.correlated.pruneSet)],
                                    names = names(tr.singletons.correlated.pruneSet)
                                  )

    pdf(paste('../OUT/SUPPLEMENTS/FIG.S.lttPlot', lambda, calibPoint, 'pdf', sep = '.'))
    layout(matrix(c(rep(2, 7), rep(1, 3)), 10))
    par(mar=c(5.1,4.1,0.1,2.1))
    plot(d18Oadj ~ Age, zachos, pch = "+", cex = 0.5, xlim = c(-60,0),
          xlab = 'Age (Mya)', ylab = expression(paste(delta^18, 'O')), axes = 0,
        ylim = rev(range(zachos$d18Oadj, na.rm = T)))
    axis(2)
    axis(3, at = seq(from = -60, to = 0, by = 1), labels = F)
    par(mar=c(1.0,4.1,4.1,2.1))
    plot(1, xlim = c(-60,0),
        ylim = c(0, max(sapply(tr.singletons.correlated.pruneSet, function(x) length(x$tip.label)))),
              xaxp = c(-60,0,60), cex.axis = 0.000001,
              xlab = '', ylab = 'Cumulative oak species by lineage',
            type = 'n',
          axes = 0)
    for(i in names(tr.singletons.correlated.pruneSet))
      ltt.lines(tr.singletons.correlated.pruneSet[[i]],
                col = tr.singletons.cols[i])
    axis(2, at = seq(from = 0, to = 150, by = 10))
    legend('topleft',
            legend = c(expression(paste('sect. ', italic(Lobatae))),
                       expression(paste('sect. ', italic(Quercus))),
                       expression(paste('sect. ', italic(Ilex))),
                       expression(paste('sect. ', italic(Cerris))),
                       expression(paste('sect. ', italic(Cyclobalanopsis)))
                     ),
            lwd = 2, col = tr.singletons.cols, bty = 'n')

    dev.off()
}}
