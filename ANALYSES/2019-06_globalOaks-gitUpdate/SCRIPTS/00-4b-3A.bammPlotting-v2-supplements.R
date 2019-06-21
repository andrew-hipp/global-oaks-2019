library(BAMMtools)
library(phyloch) # if necessary, install from source: http://www.christophheibl.de/phyloch_1.5-3.tar.gz
library(phytools)
library(magrittr)
data(strat2012) # A data frame containing the stratigraphic chart by issued in 2012 by the International Commission on Stratigraphy.

treePal = c('blue', 'gray', 'red')
treePalLabel = paste(treePal, collapse = '.')
btwBarsAdj = 0.1
colSubg = 'black'
xSubg = 72
colSect = 'black'
xSect = 77
colClade = 'gray50'
xClade = 89

do.treePlot <- TRUE
refreshEvents <- TRUE

  tr.calibrated.supplement <- read.tree(
    '../BAMM.WORKING/crown/tr.singletons.correlated.1.taxaGrepCrown.tre'
    )
  tr.quickLabeler <- radMetaOut[match(tr.calibrated.supplement$tip.label, radMetaOut$tip), c('tip', 'subgenus', 'section', 'clade')]
  write.csv(tr.quickLabeler, '../OUT/ANALYSIS.PRODUCTS/orderedSingletonTreeTips-SUPPLEMENT.csv')
    events <- getEventData(tr.calibrated.supplement, '../BAMM.WORKING/crown/event_data.txt', burnin = 0.05)
    message('I got the event data just fine')
    events$tip.label <- sapply(strsplit(events$tip.label, '|', fixed = T), function(x) {
      paste(x[1], x[length(x)], sep = '|')
      })
    events$tip.label <- gsub('|', ' | ', events$tip.label, fixed = TRUE)
    events$tip.label <- gsub('_', ' ', events$tip.label, fixed = TRUE)

## get probability of changes by node
  events.table <- do.call('rbind', events$eventData)
  events.p <-
    (table(events.table[, 'node']) / length(events$eventData)) %>%
    sort(decreasing = TRUE) # probably of changes by node

    pdf(paste('../OUT/SUPPLEMENTS/FIG.S3a.phyloplot.', treePalLabel, '.pdf', sep = ''), 8.5, 11,
        useDingbats = FALSE)
    eplot = plot.bammdata(events, labels = T,
                          xlim = c(0, 100), cex = 0.25,
                          pal = treePal
                          )
    nodelabels(node = names(events.p)[2:5], pch = 21, bg = 'black', cex = 3)
    nodelabels(node = names(events.p)[2:5],
               text = as.character(round(events.p[2:5], 2)),
               cex = 0.5,
               frame = 'n',
               col = 'white')

    #addBAMMshifts(events, par.reset=FALSE, cex=2)
    addBAMMlegend(eplot, location = c(xmin=2,3,10,65), cex.axis = 0.4)
    text(1, mean(c(10, 65)), "Net diversification rate", cex = 0.55, srt = 90)
    #add.geoscale(tr)
    axisGeo(GTS = strat2012, ages=T, cex = 0.6)
    for(i in unique(tr.quickLabeler$subgenus)[unique(tr.quickLabeler$subgenus) != '']) {
      yTemp = range(which(tr.quickLabeler$subgenus == i))
      segments(y0 = yTemp[1]+btwBarsAdj*2, y1 = yTemp[2]-btwBarsAdj*2, x0 = xSubg, x1 = xSubg, lwd = 4, col = colSubg)
      text(xSubg+1.5, mean(yTemp), labels = i, cex = 0.6, srt = 270, col = colSubg, font = 4)
      #text(xSect+1, mean(yTemp), labels = i, cex = 0.5, adj =0, col = colSect)
    }
    for(i in unique(tr.quickLabeler$section)[unique(tr.quickLabeler$section) != '']) {
      yTemp = range(which(tr.quickLabeler$section == i))
      segments(y0 = yTemp[1]+btwBarsAdj, y1 = yTemp[2]-btwBarsAdj, x0 = xSect, x1 = xSect, lwd = 2, col = colSect)
      #text(xSect+1, mean(yTemp), labels = i, cex = 0.6, srt = 270, col = colSect)
      text(xSect+1, mean(yTemp), labels = i, cex = 0.5, adj =0, col = colSect, font = 3)
    }
    for(i in unique(tr.quickLabeler$clade)[unique(tr.quickLabeler$clade) != '']) {
      yTemp = range(which(tr.quickLabeler$clade == i))
      segments(y0 = yTemp[1]+btwBarsAdj, y1 = yTemp[2]-btwBarsAdj, x0 = xClade, x1 = xClade, lwd = 2, col = colClade)
      #text(xClade+1, mean(yTemp), labels = i, cex = 0.6, srt = 270, col = colClade)
      text(xClade+1, mean(yTemp), labels = i, cex = 0.5, adj=0, col = colClade)
    }
    dev.off()
