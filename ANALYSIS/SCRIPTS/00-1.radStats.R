## RAD-seq matrix stats for paper
library(tidyverse)

pr <- function(num) prettyNum(num, digits = 4, big.mark = ',')

rad.stats <- sapply(c('s2', 's3', 's4', 's5'), function(x) {
  temp <- read.table(grep(x, dir('../../../RAD.STATS', full = T), value = T), as.is = T, header = T)
  row.names(temp) <- tidyName(row.names(temp), case = 'upper') %>% make.unique
  row.names(temp) <- gsub('FQ|BARCODESTRIPPED|TECHREP|NAMEFIXED|2012', '', row.names(temp)) %>%
    make.unique
  temp <- temp[grep('.', row.names(temp), invert = T, fixed = T), ]
  temp <- temp[radMetaOut$code, ]
})

attach(rad.stats)
rad.stats.out <- c(
  paste("Total clusters > depth 15:", dim(rads.subset$indsMatFull)[2]),
  paste("Mean clusters per individual in final dataset:", pr(mean(apply(rads.subset$indsMatFull, 1, sum) / dim(rads.subset$indsMatFull)[2])), '+/-',
                                                          pr(sd(apply(rads.subset$indsMatFull, 1, sum) / dim(rads.subset$indsMatFull)[2]))),
  paste("Aligned nucleotides in final dataset:", paste(rads.subset$matFull[1,], collapse = '') %>% nchar),
  paste("Clusters in singletons tre", sum(colSums(rads.subset$indsMat) > 14)),
  paste("Aligned nucleotides in singletons dataset:", paste(rads.subset$mat[1,which(colSums(rads.subset$indsMat) > 14)], collapse = '') %>% nchar),
  "",
  "Following lines are mean followed by sd over all individuals in study",
  paste("Raw reads:", pr(mean(s2$reads_raw)), '+/-', pr(sd(s2$reads_raw))),
  paste("Reads passing filters:", pr(mean(s2$reads_passed_filter)), '+/-', pr(sd(s2$reads_passed_filter))),
  paste("Total clusters per individual:", pr(mean(s3$clusters_total)), '+/-', pr(sd(s3$clusters_total))),
  paste("Mean depth of clusters per individual per locus:", pr(mean(s3$avg_depth_total)), '+/-', pr(sd(s3$avg_depth_total))),
  paste("Mean estimated heterozygosity:", pr(mean(s4$hetero_est)), '+/-', pr(sd(s4$hetero_est))),
  paste("Mean estimated sequencing error:", pr(mean(s4$error_est)), '+/-', pr(sd(s4$error_est))),
  paste('Number of individuals in study:'), # sum(include2019)
  paste('Number of individuals sequenced new for this study:') # grep "Oaks of the World" in metadata table$SRA_Publication
)
detach(rad.stats)
writeLines(rad.stats.out, con = '../OUT/ANALYSIS.PRODUCTS/rad.stats.txt')
