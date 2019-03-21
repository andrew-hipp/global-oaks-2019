## map divergent / introgressed loci along genome
## presupposes you have filtered and exported loci, done the topology tests,
##    compared lnl and made the lnL table

library(tidyr)
library(ncf)
library(ggplot2)
library(ggthemes)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

do.splines <- FALSE ## nothing to show... not worth doing
do.splinesAll <- FALSE ## do not do this... completely artifactual
do.chromes <- FALSE ## a nice picture that recapitulates what Andy's shown
lnLthreshold <- 4
prune.by.lnLthreshold <- TRUE

## 1. identify loci of interest, those for which the combined lnL ratio >= 2.0
locToMap.intro <-
  which(abs(lnL_mat$dumos.ratioC) >= lnLthreshold |
        abs(lnL_mat$albae.ratio) >= lnLthreshold) %>%
  row.names(lnL_mat)[.] %>%
  intersect(oaks.map$uniques$name)

message(paste('lnL individuals that map --',
        intersect(row.names(lnL_mat), oaks.map$uniques$name) %>% length
      ))
message(paste('individuals that map and have lnL difference of at least',
        lnLthreshold, '--', length(locToMap.intro)))

## 2. make table for these, with map positions
row.names(oaks.map$uniques) <- oaks.map$uniques$name
locToMap.intro <- cbind(
  oaks.map$uniques[locToMap.intro, ],
  lnL_mat[locToMap.intro, ]
)
locToMap.intro <- locToMap.intro[order(locToMap.intro$chrom, locToMap.intro$cleanStart), ]
locToMap.intro$LG <- gsub('chr', '', locToMap.intro$chrom) %>%
  as.integer

if(prune.by.lnLthreshold) {
  locToMap.intro$dumos.ratioC[abs(locToMap.intro$dumos.ratioC) < lnLthreshold] <- NA
  locToMap.intro$albae.ratio[abs(locToMap.intro$albae.ratio) < lnLthreshold] <- NA
}

locToMap.intro$albaeBinary <- ifelse(locToMap.intro$albae.ratio > 0, 'Roburoid divergence', 'Roburoid introgression')
locToMap.intro$dumosaeBinary <- ifelse(locToMap.intro$dumos.ratioC > 0, 'Dumosae divergence', 'Dumosae introgression')

## 3. aggregate loci < 200 bp from one another
diff(locToMap.intro$cleanStart) %>% abs %>% '<'(100) %>% which
# only 6 pairs to worry about... ignoring for right now... come back to this!

## 4. individual spline correlogralbae.ratioams
if(do.splines) {
  pdf('../OUT/FIGS.AND.TABLES/FIGURE.splinesByChromosome.pdf', 11, 8.5)
  layout(matrix(1:12, 3, 4))
  for(i in unique(locToMap.intro$chrom)) {
    print(paste('doing', i))
    temp <- spline.correlog(x = locToMap.intro$cleanStart[locToMap.intro$chrom == i], # dist along chromosome
                            y = rep(1, sum(locToMap.intro$chrom == i)), # fill y with 1s
                            z = locToMap.intro$monoSum[locToMap.intro$chrom == i]) # which chromosome
    plot(temp, main = i)
  }
  dev.off()
}

if(do.splinesAll) {
  pdf('../OUT/FIGS.AND.TABLES/FIGURE.splinesAll.pdf')
  spline.correlog(x = locToMap.intro$contPoint, # dist along chromosome
                  y = rep(1, dim(locToMap.intro)[1]), # fill y with 1s
                  z = locToMap.intro$monoSum) %>% # which chromosome
  plot
  dev.off()
}

##
## 5. chromosome map
if(do.chromes) {
  pdf('../OUT/FIGS.AND.TABLES/FIGURE.chromoRadMap_lnL4.pdf', 11, 8.5)
  p <- ggplot(locToMap.intro, aes(x = chrom, y = cleanStart,
                                  label = name))
  #p <- p + scale_color_gradient2(
  #                    "Monophyly\nlnL difference",
  #                    midpoint=0, low="blue", mid="white",
  #                    high="red", space ="Lab" )

  #p <- p + scale_color_gradient(
  #                    "Monophyly\nlnL difference",
  #                    low="blue", high="yellow", space ="Lab" )
  p <- p + scale_color_manual(
                      "Monophyly lnL difference",
                      values = c('black', 'red', 'brown', 'orange')
                      )

  p <- p + geom_point(data = oaks.map$uniques,
                      aes(x = chrom, y = cleanStart),
                      color = 'gray50', alpha = 0.1,
                      pch = '-', size = 9,
                      na.rm = T)
  # p <- p + geom_point(pch = "-", size = 12)
  p <- p + geom_point(data=subset(locToMap.intro, !is.na(dumosaeBinary)),
                      aes(x = LG + 0.2, y = (cleanStart), color = dumosaeBinary),
                      pch = 15, size = 2,
                      na.rm = T)

  p <- p + geom_point(data=subset(locToMap.intro, !is.na(albaeBinary)),
                      aes(x = LG - 0.2, y = (cleanStart), color = albaeBinary),
                      pch = 15, size = 2,
                      na.rm = T)

  p <- p + theme_tufte()
  p <- p + theme(legend.position = c(0.9, 0.8))
  print(p)
  dev.off()
}

## 6. Moran's I
