library(magrittr)
library(openxlsx)
library(RADami)

# 1. read full metadata file
colsToUse <- c(
  'Cleaned_NAMES-USE-THIS',
  'Determination_HIST_PIPED',
  'code', 'fastQ.file.name', 'Specimen.CODE',
  'singleTip',
  'subgenus', 'section', 'clade',
  'lociRetained-2018-07',
  'collectionDate', 'collectors', 'CollNum',
  'Georeference_source', 'country-georef', 'state_province_georef', 'County_georef', 'City_georef',
  'site',
  'livingCollectionsNumber', 'OrigSourceInfo', 'Original.Source-.Cultivated',
  'latitude_georef', 'longitude_georef',
  'voucher1', 'MOR_Acc_Num',
  'SRA_bioproject_accession', 'SRA_biosample_accession',
  'SRA_library_ID','SRA_filename', 'SRA_Run_Accession_RADseq', # Run Accession is new as of 8 April 2019
  'SRA_Publication','SRA_SEQUENCES.USED.IN'
)

# 2. add cleaned code
rad.fullDat <- read.xlsx(dir('../DATA/DATABASE', patt = 'xlsx', full = T), 1)
rad.fullDat$'fastQ.file.name' <-
  gsub('_2012.techRep', '', rad.fullDat$'fastQ.file.name', fixed = T)
rad.fullDat$code <- sapply(strsplit(rad.fullDat$'fastQ.file.name', '.', fixed = T), function(x) {
  x <- strsplit(x[1], '.', fixed = T)[[1]][1] %>%
    tidyName(case = 'upper')
  return(x)
  })

rad.fullDat$singleTip <- as.logical(rad.fullDat$singleTip)
rad.fullDat$include2019 <- as.logical(rad.fullDat$include2019)

radMetaOut <- rad.fullDat[which(rad.fullDat$include2019), colsToUse]
radMetaOut$tip <- paste(
  ifelse(is.na(radMetaOut$'Cleaned_NAMES-USE-THIS'), radMetaOut$LABEL_taxonOnly, radMetaOut$'Cleaned_NAMES-USE-THIS'),
  radMetaOut$'country-georef',
  radMetaOut$state_province_georef,
  radMetaOut$CollNum,
  radMetaOut$Specimen.CODE,
  sep = "|"
)
radMetaOut$tip <- gsub(' ', '_', radMetaOut$tip)
radMetaOut$tip <- gsub(',', '.', radMetaOut$tip)

row.names(radMetaOut) <- radMetaOut$code

radMetaOut[radMetaOut %in% c('n/a','NA')] <- ''
radMetaOut[is.na(radMetaOut)] <- ''
write.csv(radMetaOut, '../OUT/SUPPLEMENTS/TABLE.S.samples.csv', fileEncoding = "UTF-8" )
