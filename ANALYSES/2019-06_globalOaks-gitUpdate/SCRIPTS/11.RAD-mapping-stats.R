## summarize results of introgression tests, mapping back to RADs, etc.
## 2019-03-17: currently requires v2 results, being run in v2 workspace

library(magrittr)
library(grid)
library(gridExtra)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(IRanges)
library(ape)
library(phytools)
library(scatterpie)

so <- function() source('../SCRIPTS/11.RAD-mapping-stats.R', echo = T)
pN <- function(x) prettyNum(x, big.mark = ',', )
rms <- c('RAD mapping stats', '----------------')

# Read data
final.loci <- dimnames(rads.subset$indsMatFull)[[2]]

trs.final <- list(
  stem.calib = read.tree('../OUT/ANALYSIS.PRODUCTS/tr.singletons.correlated.1.taxaGrepStem.tre'),
  crown.calib = read.tree('../OUT/ANALYSIS.PRODUCTS/tr.singletons.correlated.1.taxaGrepCrown.tre')
)
# 1. Descriptive stats: RADs vs genome
#   a. How many RADs map back based on the actual loci used in the paper
oaks.map.uniques <- oaks.map$uniques[intersect(final.loci, row.names(oaks.map$uniques)), ]
rms <- c(rms, paste('Number of RADs mapping at least once:', length(intersect(row.names(oaks.blast), final.loci))))
rms <- c(rms, paste('Number of RADs mapping uniquely to one of the 12 chromosomes:', dim(oaks.map.uniques)[1]))

#   b. How many of these are in genes? -- REMAKING SOME OF WHAT'S IN 00-7.radMapping.R
oaks.map.ranges <- sapply(chromLabels, function(x) {
	IRanges(start = oaks.map.uniques$cleanStart[oaks.map.uniques$chrom == x],
					end = oaks.map.uniques$cleanEnd[oaks.map.uniques$chrom == x],
          names = oaks.map.uniques$name[oaks.map.uniques$chrom == x])
				}
			)
oaks.map.overlaps <- sapply(chromLabels, function(i) {
	findOverlaps(oaks.map.ranges[[i]], oaks.geneRanges[[i]])
	}
)

oaks.map.uniques$inGene <- F
for(i in chromLabels) {
  range.temp <- oaks.map.ranges[[i]] %>% as.data.frame
  overlap.temp <- oaks.map.overlaps[[i]] %>% as.data.frame
  oaks.map.uniques[range.temp$names[overlap.temp$queryHits], 'inGene'] <- T
}

oaks.map.geneIntersects <- data.frame(
	radsTotal = sapply(chromLabels, function(i) sum(oaks.map.uniques$chrom == i)),
	radsInGenes = oaks.map.overlaps %>% sapply(., length),
	chrL = oaks.chrLengths,
	chromLabel = chromLabels,
	row.names = chromLabels
)
oaks.map.geneIntersects$radsNotInGenes <-
	oaks.map.geneIntersects$radsTotal -
	oaks.map.geneIntersects$radsInGenes

rms <- c(rms,
          paste('By chrom totals:',
                colnames(oaks.map.geneIntersects[-c(3,4)]),
                ':',
                apply(oaks.map.geneIntersects[, -c(3,4)], 2, mean) %>% pN,
                "+/-",
                apply(oaks.map.geneIntersects[, -c(3,4)], 2, sd) %>% pN
              )
            )
rms <- c(rms,
          paste('Unique genes total:',
                sapply(oaks.map.overlaps, function(x) {
                  as.data.frame(x)$subjectHits %>%
                    unique %>%
                    length

                }) %>% sum
              ))
#   c. What is the range in percentages mapping to genes among chromosomes?
rms <- c(rms,
  paste('By chrom percent genes:',
  (oaks.map.geneIntersects[,'radsInGenes'] /
    oaks.map.geneIntersects[, 'radsTotal']) %>% mean %>% round(3),
  "+/-",
  (oaks.map.geneIntersects[,'radsInGenes'] /
    oaks.map.geneIntersects[, 'radsTotal']) %>% sd %>% round(3)
  ) # close paste
) # close c

#   d. plot rads to chromosome piecharts
p <- ggplot(oaks.map.geneIntersects,
						aes(y = radsTotal, x = chrL/20000, label = chromLabel))
#p <- p + geom_point()
#p <- p + geom_smooth(method = 'lm')
p <- p + geom_scatterpie(data = oaks.map.geneIntersects,
												aes(y = radsTotal, x = chrL/20000, r = 90),
												cols = c('radsInGenes', 'radsNotInGenes'), alpha = 1)
p <- p + coord_fixed()
p <- p + geom_text_repel(point.padding = 0.5, min.segment.length = 3)
p <- p + labs(x = "Chromosome length * 20,000 (bp)",
							y = "Number of RAD loci mapped to chromosome")
p <- p + scale_fill_manual(values=c('black', 'gray80'),
							             name="Genomic position\nof RAD loci",
							             labels=c("Overlapping a gene", "Not overlapping genes"))
p <- p + theme(legend.position = c(0.85, 0.15),
							 legend.background = element_blank())

pdf('../OUT/FIGS.AND.TABLES/FIGURE.radsVsChromosomes.pdf')
print(p)
dev.off()

#   e. Distribution of RADs on the chromosome: are there conspicuous hotspots or
#       missing areas?
#   NOT DONE FOR NOW -- WOULD BE MORE INTERESTING IF WE HAD CENTROMERES, BUT WE DON'T

# 2. Introgression analyses
#   a. What is the correlation between the introgression / divergence signal for
#       the two clades investigated?

lnL_mat$dumosBin <- ifelse(lnL_mat$dumos.ratioC > 2, 1, 0)
lnL_mat$albaeBin <- ifelse(lnL_mat$albae.ratio > 2, 1, 0)
#### no binning does any good at this point

temp <- cor.test(lnL_mat$dumos.ratioC, lnL_mat$albae.ratio)
rms <- c(rms, paste("Correlation between Albae and Dumosae introgression by locus: r =",
        round(temp$estimate, 4), ', p =',
        round(temp$p.value, 4), ', df =', temp$parameter))

#   b. Is there a difference in signal between genes / nongenes for each clade?
lnL_mat.byMap <- lnL_mat[intersect(row.names(lnL_mat), row.names(oaks.map.uniques)), ]
lnL_mat.byMap$inGene <- oaks.map.uniques[row.names(lnL_mat.byMap), 'inGene']

temp <- aov(albae.ratio ~ inGene, lnL_mat.byMap) %>% anova
rms <- c(rms, paste("Effect of presence in gene on Albae introgression: F =",
                    round(temp$'F value'[1], 4),
                    '(', paste(temp$'Df', collapse= ','),')',
                    'P =', round(temp$"Pr(>F)"[1], 4)
                  ))

temp <- aov(dumos.ratioC ~ inGene, lnL_mat.byMap) %>% anova
rms <- c(rms, paste("Effect of presence in gene on Dumosae introgression: F =",
                    round(temp$'F value'[1], 4),
                    '(', paste(temp$'Df', collapse= ','),')',
                    'P =', round(temp$"Pr(>F)"[1], 4)
                  ))

rms <- c(rms, paste('total loci with lnL diff >= 2, at least one test:',
                    ((lnL_mat.byMap$albae.ratio %>% abs %>% '>='(2)) | (lnL_mat.byMap$dumos.ratioC %>% abs %>% '>='(2))) %>% sum(na.rm = T))
                  )
rms <- c(rms, paste('total loci with lnL diff >= 2, both tests:',
                    ((lnL_mat.byMap$albae.ratio %>% abs %>% '>='(2)) & (lnL_mat.byMap$dumos.ratioC %>% abs %>% '>='(2))) %>% sum(na.rm = T))
                  )
rms <- c(rms, paste('total loci with lnL diff >= 2, Dumosae:',
                    (lnL_mat.byMap$dumos.ratioC %>% abs %>% '>='(2)) %>% sum(na.rm = T))
                  )

rms <- c(rms, paste('total loci with lnL diff >= 2, Albae:',
                    (lnL_mat.byMap$albae.ratio %>% abs %>% '>='(2)) %>% sum(na.rm = T))
                  )
#   c. Are the divergence signals spatially autocorrelated on the genome?
### -- NO -- SEE SPLINES

# 3 Phylogenetic signal analyses (quartet distances)
#   a. Is phylogenetic signal correlated with being in a gene?
oaks.map.quartetD <- oaks.map$quartetD[intersect(row.names(oaks.map$quartetD), row.names(oaks.map.uniques)),]
oaks.map.quartetD$inGene <- oaks.map.uniques[row.names(oaks.map.quartetD), 'inGene']
temp = aov(QuartetSimilarity ~ inGene, oaks.map.quartetD) %>% anova
rms <- c(rms, paste("Effect of presence in gene on quartet similarity to consensus tree: F=",
                    round(temp$'F value'[1], 4),
                    '(', paste(temp$'Df', collapse= ','),')',
                    'P =', round(temp$"Pr(>F)"[1], 4)
                    ))

temp <- lm(QuartetSimilarity ~ LocusTreeLength, oaks.map.quartetD) %>% summary
rms <- c(rms, paste("Effect of locus taxon sampling on quartet similarity: r2=",
                    round(temp$r.squared, 4),
                    'B1 =', round(temp$coefficients['LocusTreeLength', 'Estimate'], 4),
                    'P =',
                    round(temp$coefficients['LocusTreeLength', 'Pr(>|t|)'], 4)))

#   b. Is phylogenetic signal genomically autocorrelated?
### -- NO -- SEE SPLINES

# 4. Phyparts analyses
#   a. Regress percent loci supporting a node with the node age and clade size
#       (two plots). Label nodes... are there outliers?
      # i. get ages for all relevant nodes on stem and crown calib trees
      phyPartNodes.new <- phyPartNodes
      phyPartNodes.new$nodeCalibStem <- sapply(phyPartNodes.new$nodeCalibGrep, function(x) {
        findMRCA(trs.final$stem.calib, tips = grep(x, trs.final$stem.calib$tip.label, value = T))
      })
      phyPartNodes.new$nodeCalibCrown <- sapply(phyPartNodes.new$nodeCalibGrep, function(x) {
        findMRCA(trs.final$crown.calib, tips = grep(x, trs.final$crown.calib$tip.label, value = T))
      })
      temp <- sapply(phyPartNodes.new$nodeCalibGrep, function(x) {
        findMRCA(trs.final$crown.calib, tips = grep(x, trs.final$crown.calib$tip.label, value = T)) %>%
          getParent(tree = trs.final$crown.calib)
      })
      phyPartNodes.new$nodeCalibCrown.stem[!sapply(temp, is.null)] <- unlist(temp[!sapply(temp, is.null)])
      phyPartNodes.new$depthCalibStem <- max(node.depth.edgelength(trs.final$stem.calib)) -
        node.depth.edgelength(trs.final$stem.calib)[phyPartNodes.new$nodeCalibStem]
      phyPartNodes.new$depthCalibCrown <- max(node.depth.edgelength(trs.final$crown.calib)) -
        node.depth.edgelength(trs.final$crown.calib)[phyPartNodes.new$nodeCalibCrown]
      phyPartNodes.new$depthCalibCrown.stemAge[!is.na(phyPartNodes.new$nodeCalibCrown.stem)] <-
        max(node.depth.edgelength(trs.final$crown.calib)) -
        node.depth.edgelength(trs.final$crown.calib)[phyPartNodes.new$nodeCalibCrown.stem[!is.na(phyPartNodes.new$nodeCalibCrown.stem)]]

      # ii. get clade sizes (samples only) for all relevant nodes
      phyPartNodes.new$DescendantTips <- sapply(phyPartNodes.new$nodeCalibStem, function(x) {
        sum(getDescendants(trs.final$stem.calib, x) <= length(trs.final$stem.calib$tip.label))
      })

      # iii. Calculate support and reject numbers for each node
      phyPartNodes.new$concord <- sapply(phyPartNodes.new$node, function(x) {
        sum(phyParts.df$info[phyParts.df$node == x] == 'concord')
      })
      phyPartNodes.new$conflict <- sapply(phyPartNodes.new$node, function(x) {
        sum(phyParts.df$info[phyParts.df$node == x] == 'conflict')
      })
      phyPartNodes.new$TotalLoci <-
        phyPartNodes.new$concord + phyPartNodes.new$conflict
      phyPartNodes.new$ProportionLociConcordant <-
        phyPartNodes.new$concord /
          (phyPartNodes.new$concord + phyPartNodes.new$conflict)
      phyPartNodes.new$Taxon <- substr(phyPartNodes.new$taxon, 5, nchar(as.character(phyPartNodes.new$taxon)))

      # iv. regressions
      pp1 <- ggplot(phyPartNodes.new, aes(y = ProportionLociConcordant, x = DescendantTips,
                                          label = Taxon))
      pp1 <- pp1 + geom_smooth(method='lm')
      pp1 <- pp1 + ylab('Proportion of loci concordant with crown node')
      pp1 <- pp1 + xlab('Number of tips sampled in clade')
      pp1 <- pp1 + geom_point(aes(size = TotalLoci))
      pp1 <- pp1 + geom_label_repel()
      pp1 <- pp1 + theme(legend.position = 'none')
      pp1 <- pp1 + ylim(c(-0.2,0.9))

      pp2 <- ggplot(phyPartNodes.new, aes(y = ProportionLociConcordant, x = depthCalibCrown,
                    label = Taxon))
      pp2 <- pp2 + geom_smooth(method='lm')
      pp2 <- pp2 + ylab('')
      pp2 <- pp2 + xlab('Crown age of clade (Mya), estimated under crown calibration')
      pp2 <- pp2 + geom_point(aes(size = TotalLoci))
      pp2 <- pp2 + geom_label_repel()
      pp2 <- pp2 + theme(legend.position = c(0.9, 0.1))
      pp2 <- pp2 + ylim(c(-0.2,0.9))

      pdf('../OUT/SUPPLEMENTS/FIGURE.S.phyPartPreds.pdf', 11, 8.5, useDingbats = FALSE)
      grid.arrange(pp1, pp2, nrow = 1)
      dev.off()

#     basic stats
      rms <- c(rms, paste('Phyparts by node,', c('concordant loci:', 'conflicting loci:', 'total loci:', 'proportion loci concordant:'),
              sapply(phyPartNodes.new[, c('concord', 'conflict', 'TotalLoci', 'ProportionLociConcordant')], function(x) paste(round(mean(x, na.rm = T), 4), '+/-', round(sd(x, na.rm = T), 4), '-- range = ',
              paste(round(range(x, na.rm = T), 4), collapse = '-')))
            ))


#   b. what are the simple effects of node depth and clade size on percent supporting?
      temp <- cor.test(phyPartNodes.new$ProportionLociConcordant, phyPartNodes.new$depthCalibCrown)
      rms <- c(rms, paste('Correlation between proportion of loci concordant and crown age: r =',
                      round(temp$estimate, 4), ', p =', round(temp$p.value, 4)))
      temp <- cor.test(phyPartNodes.new$ProportionLociConcordant, phyPartNodes.new$DescendantTips)
      rms <- c(rms, paste('Correlation between proportion of loci concordant and number of descendents: r =',
                      round(temp$estimate, 4), ', p =', round(temp$p.value, 4)))
      temp <- cor.test(phyPartNodes.new$depthCalibCrown, phyPartNodes.new$DescendantTips)
      rms <- c(rms, paste('Why not do multiple regression? expected correlation between descendants and crown age: r =',
      round(temp$estimate, 4), ', p =', round(temp$p.value, 4)))

#   c. Important to note that all these nodes are supported strongly... but there
#       is still substantial genomic conflict around them, due to ILS or introgression

rad.map.stats <- rms
writeLines(rad.map.stats, '../OUT/ANALYSIS.PRODUCTS/rad.map.stats.txt')
write.csv(phyPartNodes.new, '../OUT/SUPPLEMENTS/TABLE.phyParts.csv')
