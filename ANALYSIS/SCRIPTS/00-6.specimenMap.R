require(maps)
require(mapdata)
library(ggplot2)
library(ggrepel)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

radMetaMap <- radMetaOut
radMetaMap$Latitude <- as.numeric(radMetaMap$latitude_georef)
radMetaMap$Longitude <- as.numeric(radMetaMap$longitude_georef)
radMetaMap$Latitude[is.na(radMetaMap$Latitude)] <-
  as.numeric(radMetaMap$latitude_georef)[is.na(radMetaMap$Latitude)]
radMetaMap$Longitude[is.na(radMetaMap$Longitude)] <-
  as.numeric(radMetaMap$longitude_georef)[is.na(radMetaMap$Longitude)]
names(radMetaMap)[names(radMetaMap) == 'section'] <- 'Section'
names(radMetaMap)[names(radMetaMap) == 'subgenus'] <- 'Subgenus'
radMetaMap <- radMetaMap[grep('Quercus', radMetaMap$tip), ]

global <- map_data("world")
gg1 <- ggplot()
gg1 <- gg1 + geom_polygon(data = global, aes(x=long, y = lat, group = group),
              color = "black", fill = "gray80")
gg1 <- gg1 + coord_fixed(1.3)
gg1 <- gg1 +
  geom_point(data=radMetaMap,
             aes(x=Longitude, y=Latitude,
                 color = Section,
                 shape = Subgenus)
               )
gg1 <- gg1 + xlim(-150, 150) + ylim(0,100)
gg1 <- gg1 + theme(legend.position = 'top')
gg1 <- gg1 + scale_color_manual(values = cbbPalette)
pdf('../OUT/FIGS.AND.TABLES/FIGURE.samplingMap-v3.pdf', 11, 8.5)
print(gg1)
dev.off()
