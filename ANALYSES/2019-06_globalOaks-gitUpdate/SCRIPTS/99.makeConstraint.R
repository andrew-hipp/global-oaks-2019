makeConstraint <- function(taxa, constraintList) {
  ## just does a simple, unnested constraint, with each vector in constraintList set up as a clade

  # 0. subset constraintList to elements that are in taxa

  #if(!all(unlist(constraintList) %in% taxa)) {
  #  print(c('\n', 'Missing taxa', '-----'))
  #  print(setdiff(unlist(constraintList), taxa))
  #}
  constraintList <- lapply(constraintList, intersect, y = taxa)

  # 1. make each constraint
  outCon <- lapply(constraintList, paste, collapse = ',')
  outCon <- lapply(outCon, function(x) paste('(',x,')', sep = ''))
  outCon <- do.call('paste', c(outCon, sep = ','))

  # 2. make entire vector using the remainder of the taxa
  outCon <- paste(paste(setdiff(taxa, unlist(constraintList)), collapse = ','),
                  outCon,
                  sep = ',')
  outCon <- paste('(', outCon, ');', sep = '')

  # 4. ... and return the result
  return(outCon)
}
