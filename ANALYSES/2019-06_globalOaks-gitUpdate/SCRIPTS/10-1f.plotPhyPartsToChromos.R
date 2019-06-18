library(ggplot2)
library(magrittr)
library(ggthemes)

# 1. Assign key nodes
##  -- this is done by manually inspecting the trees, which is pokey!
phyPartNodes <- matrix(
                  c('16. Quercus', 5, 'Quercus',
                    '15. subg. Quercus', 78, 'Quercus_rubra|Quercus_alba',
                    '14. subg. Cerris', 6, 'Quercus_gilva|Quercus_senescens',
                    '13. sect. Quercus', 89, 'Quercus_durata|Quercus_arizonica',
                    '12. sect. Protobalanus', 180, 'Quercus_cedrosensis|Quercus_chrysolepis',
                    '11. sect. Ponticae', 81, 'Quercus_pontica|Quercus_sadleriana',
                    '10. sect. Virentes', 83, 'Quercus_fusiformis|Quercus_virginiana',
                    '09. Roburoid white oaks', 103, 'Quercus_dentata|Quercus_faginea',
                    '08. Mexican white oaks', 145, 'Quercus_germana|Quercus_arizonica',
                    '07. sect. Lobatae', 184, 'Quercus_parvula|Quercus_radiata',
                    '06. Laurifoliae', 194, 'Quercus_ilicifolia|Quercus_incana',
                    '05. Rubrae', 188, 'Quercus_coccinea|Quercus_acerifolia',
                    '04. Mexican red oaks', 210, 'Quercus_affinis|Quercus_radiata',
                    '03. sect. Cerris', 41, 'Quercus_chenii|Quercus_trojana',
                    '02. sect. Ilex', 54, 'Quercus_franchetii|Quercus_ilex',
                    '01. sect. Cyclobalanopsis', 7, 'Quercus_gilva|Quercus_bella'
                ),
                nrow = 16,
                ncol = 3,
                byrow = T,
                dimnames = list(NULL, c('taxon', 'node', 'nodeCalibGrep'))
              ) %>% as.data.frame

phyPartNodes$node <- as.character(phyPartNodes$node)
phyPartNodes$nodeCalibGrep <- as.character(phyPartNodes$nodeCalibGrep)

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
p <- p + scale_color_manual('Node support', values = c('black', 'orange'))
p <- p + xlab('') + ylab('')
p <- p + theme_tufte()
p <- p + theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank()
             )
pdf('../OUT/FIGS.AND.TABLES/FIGURE.phyParts.by.node-selected-v3.pdf', 7, 4)
print(p)
dev.off()
