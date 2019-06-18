library(BAMMtools)
library(phyloch) # if necessary, install from source: http://www.christophheibl.de/phyloch_1.5-3.tar.gz
library(phytools)
data(strat2012) # A data frame containing the stratigraphic chart by issued in 2012 by the International Commission on Stratigraphy.

btwBarsAdj = 0.1
colSubg = 'black'
xSubg = 72
colSect = 'black'
xSect = 77
colClade = 'gray50'
xClade = 90

do.ltt <- FALSE
do.treePlot <- TRUE
do.rtt <- TRUE
for(calibPoint in c('crown', 'stem', 'crown-sampleProps', 'stem-sampleProps')) {
  tr.calibrated <- read.tree(
    paste('../OUT/ANALYSIS.PRODUCTS/tr.singletons.correlated.1.taxaGrep',
    paste(toupper(substr(calibPoint, 1,1)), substr(calibPoint, 2, nchar(calibPoint)), sep = ''),
    '.tre', sep = '')
  )
  tr.quickLabeler <- radMetaOut[match(tr.calibrated$tip.label, radMetaOut$tip), c('tip', 'subgenus', 'section', 'clade')]
  write.csv(tr.quickLabeler, '../OUT/ANALYSIS.PRODUCTS/orderedSingletonTreeTips-v2.csv')
  events <- getEventData(tr.calibrated, paste('../BAMM.WORKING/', calibPoint, '/event_data.txt', sep = ''), burnin = 0.05)
  events$tip.label <- sapply(strsplit(events$tip.label, '|', fixed = T), function(x) {
    paste(x[1], x[length(x)], sep = '|')
    })

  if(do.treePlot) {
    pdf(paste('../OUT/SUPPLEMENTS/FIG.S.phyloplot.withDiv.', calibPoint, '.pdf',sep = ''), 8.5, 11)
    eplot = plot(events, labels = T, xlim = c(0, 100), cex = 0.25)
    #addBAMMshifts(events, par.reset=FALSE, cex=2)
    addBAMMlegend(eplot, location = c(xmin=2,3,10,65), cex.axis = 0.4)
    #text(2, 66, "Net diversification", cex = 0.4)
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
    if(do.ltt) ltt.lines(tr.calibrated, col = 'gray', lty = 'dashed', backward = FALSE, lwd = 2)
    dev.off()
  }

  if(do.rtt) {
    st <- max(branching.times(tr.calibrated))
    pdf(paste('../OUT/FIGS.AND.TABLES/FIG.S.rtt-comboNetDivAll.', calibPoint, '.pdf', sep = ''))
    #layout(matrix(1:3, 1))
    plotRateThroughTime(events, ratetype = 'netdiv', avgCol="black",
      start.time=st, ylim=c(0,1), cex.axis = 1, intervalCol = 'gray80', intervals = c(0.025, 0.975), opacity = 0.8)
    plotRateThroughTime(events, ratetype = 'netdiv', avgCol="blue",
      node = findMRCA(tr.calibrated, grep('Quercus_alba|Quercus_velutina', tr.calibrated$tip.label, value = T)),
      start.time=st, add = T, cex.axis=1, intervalCol='lightblue', intervals=c(0.025, 0.975), opacity=0.8)
    #text(x=30, y= 0.8, label="American oaks", font=4, cex=2.0, pos=4)
    plotRateThroughTime(events, ratetype = 'netdiv', avgCol="darkorange",
      node = findMRCA(tr.calibrated, grep('Quercus_ilex|Quercus_gilva', tr.calibrated$tip.label, value = T)),
      start.time=st, add = T, intervalCol='orange', intervals=c(0.025, 0.975), opacity=0.8)
    #text(x=30, y= 0.8, label="Eurasian oaks", font=4, cex=2.0, pos=4)
    legend('topleft', bty = 'n', legend=c('Global oaks', 'subg. Quercus', 'subg. Cerris'), lwd = rep(2,3), col = c('black', 'blue', 'darkorange'))
    dev.off()
  }
} # close calibPoint
