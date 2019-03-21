library(ggplot2)
library(magrittr)
library(ggthemes)

do.quickPlot = FALSE

## quick and dirty plot!
if(do.quickPlot) {
  p <- ggplot(phyParts.df, aes(x = genomeLinStart, y = node))
  p <- p + geom_point(aes(color = info), pch = '|')
  pdf('../OUT/FIGS.AND.TABLES/FIGURE.phyParts.by.node-quickAll.pdf', 11, 8.5)
  print(p)
  dev.off()
}

# NOW DO IT RIGHT

# 1. Assign key nodes
phyPartNodes <- matrix(
                  c('01. Quercus', 3,
                    '02. subg. Quercus', 4,
                    '03. subg. Cerris', 186,
                    '04. sect. Quercus', 8,
                    '09. sect. Lobatae', 112,
                    '06. sect. Virentes', 100,
                    '05. sect. Protobalanus', 107,
                    '13. sect. Cerris', 213,
                    '14. sect. Ilex', 188,
                    '15. sect. Cyclobalanopsis', 226,
                    #'sect. Ponticae',
                    '10. Laurifoliae', 159,
                    '11. Rubrae', 174,
                    '12. Mexican red oaks', 116,
                    '08. Mexican white oaks', 12,
                    '07. Roburoid white oaks', 63
                ),
                nrow = 15,
                ncol = 2,
                byrow = T,
                dimnames = list(NULL, c('taxon', 'node'))
              ) %>% as.data.frame

phyPartNodes$node <- phyPartNodes$node %>% as.character

# 2. Replace those nodes in phyParts.df
phyParts.dfPlot <- phyParts.df
for(i in 1:dim(phyPartNodes)[1]) {
  phyParts.dfPlot$node[phyParts.dfPlot$node == phyPartNodes[i, 'node']] <- as.character(phyPartNodes[i, 'taxon'])
}

# 3. Plot just the key nodes
p <- ggplot(subset(phyParts.dfPlot, node %in% phyPartNodes$taxon),
            aes(x = genomeLinStart, y = node)
          )
p <- p + geom_vline(data = chromStarts[chromStarts$Group.1 != 'chr01', ],
                    aes(xintercept = x),
                    lwd = 0.25,
                    color = 'gray80')
p <- p + geom_point(aes(color = info), pch = '|')
p <- p + scale_color_manual('Node support', values = c('black', 'gold'))
p <- p + xlab('') + ylab('')
p <- p + theme_tufte()
p <- p + theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank()
             )
pdf('../OUT/FIGS.AND.TABLES/FIGURE.phyParts.by.node-selected-v1.pdf', 11, 4)
print(p)
dev.off()
