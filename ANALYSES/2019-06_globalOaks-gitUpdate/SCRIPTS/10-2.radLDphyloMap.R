## Map quartet distance on genome
## hypothesis: genomic autocorrelation of phylogenetic
  ## informativeness, suggesting LD at large scales
  ## and over deep phylogenetic time (~ 50 million yrs)

library(ncf)

# 1. Read in global (concatenation) tree,
#
# 2. Calculate quartet distances between the tree for each
#    locus and the global (consensus) tree, unconstrained
#
# 3. Calculate pairwise distances along genome between RAD
#    markers
#
# 4. Aggregate markers < 200 bp, averaging quartet distances
#
# 5. Iterate over 2 - 3 to ensure that markers are >= 200 bp apart
#
# 6. Map quartet distances on chromosomes
#
# 7. Moran's I for quartet distances: what is the scale of LD at whole-phylogeny scale?
rad.splines <- spline.correlog(x = MAPposition, y = 0, w = qDistance)
