# function to test for number of taxa
taxaTest <- function(taxa, taxList, threshold, simpleOut = TRUE) {
  ## taxa = list of taxa to compare to
  ## taxList = list of character vectors of tips
  ## threshold = integer vector of length = length(taxList)
  out <- sapply(seq(length(taxList)), function(i) {
    length(intersect(taxa, taxList[[i]])) >= threhold[i]
  })
  if(simpleOut) out <- all(out)
  return(out)
}
