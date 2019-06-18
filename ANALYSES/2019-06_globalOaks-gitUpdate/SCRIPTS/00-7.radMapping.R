#library(chromoMap)
library(LinkageMapView)
library(magrittr)
library(grid)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(IRanges)
library(scatterpie)

## data columns needed for lmv.linkage.plot:
# 1. Required, linkage group name. This will be the title
# 2. Required, position - must be in numerical order ascending within linkage group name.
# 3. Required, locus - marker name at this position.
# 4. Optional, segcol - color for the line across the chromosome at this marker.

oaks.genomeScaffs <- read.xlsx('../DATA/41477_2018_172_MOESM4_ESM.xlsx', 4)
oaks.genomeScaffs$Pseudomolecule.ID <- paste('chr', oaks.genomeScaffs$Pseudomolecule.ID, sep = '') %>%
	gsub(pattern = 'chr', replacement = 'chr0', x = .) %>%
	gsub(pattern = '01', replacement = '1', x = .)

oaks.genomeScaffs$Pseudomolecule.ID[oaks.genomeScaffs$Pseudomolecule.ID == 'chr1'] <- 'chr01'
oaks.chrLengths <- aggregate(oaks.genomeScaffs$'Scaffold.length.(bp)', list(oaks.genomeScaffs$Pseudomolecule.ID), sum)
oaks.chrLengths <- structure(oaks.chrLengths$x, names = oaks.chrLengths$Group.1)
oaks.geneModels <- read.xlsx('../DATA/41477_2018_172_MOESM3_ESM.xlsx', 1)
names(oaks.geneModels) <-
  read.xlsx('../DATA/41477_2018_172_MOESM3_ESM.xlsx', 2)$cleanLabel
oaks.geneModels$cleanStart <- apply(oaks.geneModels[c('geneStart', 'geneEnd')],1,min)
oaks.geneModels$cleanEnd <- apply(oaks.geneModels[c('geneStart', 'geneEnd')],1,max)
oaks.geneModels$pseudomolecule <- gsub('Chr', 'chr0', oaks.geneModels$pseudomolecule)
oaks.geneModels$pseudomolecule <- gsub('chr010', 'chr10', oaks.geneModels$pseudomolecule)
oaks.geneModels$pseudomolecule <- gsub('chr011', 'chr11', oaks.geneModels$pseudomolecule)
oaks.geneModels$pseudomolecule <- gsub('chr012', 'chr12', oaks.geneModels$pseudomolecule)
oaks.geneRanges <- sapply(unique(grep('chr', oaks.geneModels$pseudomolecule, value = T)), function(x) {
	IRanges(oaks.geneModels$cleanStart[oaks.geneModels$pseudomolecule == x],
					oaks.geneModels$cleanEnd[oaks.geneModels$pseudomolecule == x])
				}
			)
oaks.blast <- cbind(read.delim('../DATA/blastN_oaksall_vs_12chromos_percentID80_eval10-5_alnpercent80_besthit_andrew.txt', row.names = 1),
					read.delim('../DATA/blastN_oaksall_vs_12chromos_percentID80_eval10-5_alnpercent80_count_andrew.txt', row.names = 1)
					)
oaks.map <- list(
  uniques = data.frame(name = row.names(oaks.blast)[oaks.blast$Nber_significant_aln == 1],
                  chrom = gsub('Qrob_Chr', 'chr', oaks.blast$Oak_chromosome_ID[oaks.blast$Nber_significant_aln == 1]),
				  start = oaks.blast$aln_start[oaks.blast$Nber_significant_aln == 1],
					end = oaks.blast$aln_end[oaks.blast$Nber_significant_aln == 1],
				  data = 1
				  )
  )
oaks.map$uniques$cleanStart <- apply(oaks.map$uniques[c('start', 'end')], 1, min)
oaks.map$uniques$cleanEnd <- apply(oaks.map$uniques[c('start', 'end')], 1, max)

chromLabels <- unique(oaks.map$uniques$chrom) %>%
	sort %>%
	grep(pattern = 'chr', value = T)

oaks.map$ranges <- sapply(chromLabels, function(x) {
	IRanges(oaks.map$uniques$cleanStart[oaks.map$uniques$chrom == x],
					oaks.map$uniques$cleanEnd[oaks.map$uniques$chrom == x])
				}
			)
oaks.map$overlaps <- sapply(chromLabels, function(i) {
	findOverlaps(oaks.map$ranges[[i]], oaks.geneRanges[[i]])
	}
)
oaks.map$geneIntersects <- data.frame(
	radsTotal = sapply(chromLabels, function(i) sum(oaks.map$uniques$chrom == i)),
	radsInGenes = oaks.map$overlaps %>% sapply(., length),
	chrL = oaks.chrLengths,
	chromLabel = chromLabels,
	row.names = chromLabels
)
oaks.map$geneIntersects$radsNotInGenes <-
	oaks.map$geneIntersects$radsTotal -
	oaks.map$geneIntersects$radsInGenes

for(i in 1:2) oaks.map$uniques[[i]] <- as.character(oaks.map$uniques[[i]])
oaks.map$uniques <- oaks.map$uniques[grep('chr', oaks.map$uniques$chrom), ]
oaks.map$lmv <- data.frame(oaks.map$uniques[c('chrom', 'start', 'name')])
names(oaks.map$lmv) <- c('lgName', 'position', 'locus')
oaks.map$lmv <- oaks.map$lmv[order(oaks.map$lmv$lgName, oaks.map$lmv$position), ]
#oaks.map$uniques$chrom <- gsub('chr0', 'chr', oaks.map$uniques$chrom)
#chromoMap(oaks.map$uniques,type = "annotation")

oaks.map$dists <- sapply(unique(oaks.map$lmv$lgName), function(x) {
	diff(oaks.map$lmv$position[oaks.map$lmv$lgName == x])
})

pdf('oak.rad.hists.pdf', 11, 8.5)
layout(matrix(1:12, 3))
sapply(names(oaks.map$dists), function(x) hist(log10(oaks.map$dists[[x]]), 100, main = x, xlab = 'log10 - inter-RAD distance'))
dev.off()

pdf('oak.rad.hist.all.pdf')
grid.newpage()
pushViewport(viewport())
oaks.map$dists %>%
	unlist %>%
	log10 %>%
	hist(breaks = 100, xlab = 'log10 - inter-RAD distance, all data', main = '')
vp <- viewport(3, 1500, 3, 1500)
pushViewport(vp)
oaks.map$dists %>%
	unlist %>%
	.[. > 100] %>%
	log10 %>%
	hist(breaks = 100, xlab = 'log10 - inter-RAD distances > 100 bp', main = '')
dev.off()

oaks.numRadsMapped = sapply(oaks.map$dists, length) + 1
pdf('../OUT/SUPPLEMENTS/FIGURE.S.chrom.barplot-horiz.pdf')
barplot(rev(oaks.numRadsMapped), horiz = TRUE,
				las = 2,
				xlab = "Number of RAD-seq loci mapped uniquely to each chromosome")
dev.off()

percentGenomeGenes <-
	(oaks.geneModels$geneL[grep('chr', oaks.geneModels$pseudomolecule)] %>% sum) /
	(oaks.chrLengths %>% sum)
