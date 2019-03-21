library(RADami)

lociToWrite <- which(taxaSums[, 'rob.dumTest'] == 1) %>% names
locDir <- format(Sys.time(), "locfiles_%Y-%m-%d")

# only write loci if directory doesn't exist;
# if you want to write files, delete old directory or create new directory name
if(!locDir %in% dir()) {
  dir.create(locDir)
  for(i in lociToWrite) {
    rad2phy(rads.subset$matFull,
            inds = which(rads.subset$indsMatFull[, i]) %>% names,
            loci = i,
            padding = 100,
            outfile = paste(locDir, '/', i, '.phy', sep = ''),
            logfile = NA
          )
        }
      } # close if(!locDir)
