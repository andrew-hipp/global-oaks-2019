## convert raxml-style long-format pairwise distance table into a square distance matrix

rax.conv <- function(file = '../DATA/RAxML_distances.oaksall_v1_2.m15.singles.2019-03-13.dists.txt',
                    convertTable = radMetaOut, # set to NA if you don't have to change the names
                    convertFrom = 'code',
                    convertTo = c('Cleaned_NAMES-USE-THIS', 'Specimen.CODE'),
                    toDelim = '|',
                    diag.fill = 0) {
  dat <- read.table(file, as.is = T)
  allTips <- sort(unique(unlist(dat[1:2])))
  out <- matrix(NA, length(allTips), length(allTips),
                dimnames = list(allTips, allTips))
  for(i in 1:dim(dat)[1]) out[dat[i, 'V1'], dat[i, 'V2']] <-
                          out[dat[i, 'V2'], dat[i, 'V1']] <-
                          dat[i, 'V3']
  diag(out) <- diag.fill
  if(!is.na(convertTable[1])) dimnames(out)[[1]] <- dimnames(out)[[2]] <-
                                apply(convertTable[dimnames(out)[[1]], convertTo],
                                1,
                                paste, collapse = toDelim)
  out
  }

rads.singletons.distML <- rax.conv()
write.csv(rads.singletons.distML, '../OUT/ANALYSIS.PRODUCTS/singletonsML.dists.csv')
