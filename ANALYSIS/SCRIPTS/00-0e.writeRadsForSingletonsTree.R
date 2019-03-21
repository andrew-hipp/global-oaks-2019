require(RADami)

minInds = 15 # minimum indls per locus
phyFile = '../TREES/raxml.2019-03-13/oaksall_v1_2.m15.singles.2019-03-13.phy'

rad2phy(rads.subset$mat,
        loci = names(which(colSums(rads.subset$indsMat) >= minInds)),
        outfile = phyFile,
        padding = 100
      )

message('run raxml or whatever you like on the phylip file exported')

## first line finds best tree with RELL boots;
## second line writes bootstrap bipartitions to tree file
# raxmlHPC-PTHREADS-AVX -f D -T 8 -p 12345 -m GTRCAT -n oaksall_v1_2.m15.2019-03-13.singles-v4 -s oaksall_v1_2.m15.singles.2019-03-13.phy -g taxaConstraint.tre
# raxmlHPC-PTHREADS-AVX -f b -t RAxML_bestTree.oaksall_v1_2.m15.2019-03-13.singles-v4 -z RAxML_rellBootstrap.oaksall_v1_2.m15.2019-03-13.singles-v4 -m GTRCAT -n 2019-03-13.singles.v4.tre
