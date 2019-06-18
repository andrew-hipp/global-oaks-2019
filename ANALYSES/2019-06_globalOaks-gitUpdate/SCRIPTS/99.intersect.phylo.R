intersect.phylo <- function(phy1, phy2, taxSet = NA, minTaxa = 4) {
  ## returns the intersection of two phylogenies
  shared.tips <- intersect(phy1$tip.label, phy2$tip.label)
  if(!is.na(taxSet[1])) shared.tips <- intersect(shared.tips, taxSet)
  if(length(shared.tips) < minTaxa) out <- paste('Too short: intersect taxa =', length(shared.tips)) else
    out <- list(drop.tip(phy1, which(!phy1$tip.label %in% shared.tips)),
                drop.tip(phy2, which(!phy2$tip.label %in% shared.tips))
                )
  return(out)
}
