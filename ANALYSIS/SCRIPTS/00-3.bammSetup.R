## just to set up for BAMM -- priors block is pasted in
library(BAMMtools)

## gets BAMM priors on the lowest CV tree
for(j in calibOptions) {
  setBAMMpriors(tr.singletons.correlated[[as.character(lambdaVect.min)]][[j]],
              outfile = paste('../OUT/ANALYSIS.PRODUCTS/lambda',
                              as.character(lambdaVect.min), j,
                              'bammPriorsBlock.txt', sep = '_')
                            ) # close setBAMMpriors
}
## run using:
## /usr/local/bin/bamm -c diversification_correlated.txt
