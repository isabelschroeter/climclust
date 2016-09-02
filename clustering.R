


library(colormap)
library(ggplot2)
library(raster)
library(dplyr)
library(tidyr)
library(rgdal)
library(gridExtra)
library(grid)
library(spdep)

setwd("E:/climate_classification")


# climate data
s <- getData("worldclim", var="bio", res=10)

# cali border
b <- readOGR("E:/california_phylomodelling/Geographic_Subdivisions_of_California_TJMII_v2_060415",
             "Geographic_Subdivisions_of_California_TJMII_v2_060415")
b <- spTransform(b, crs(r))

# munge
r <- mask(crop(s, extent(b)), b)
s <- mask(crop(s, extent(b)), b)
r <- subset(r, c("bio1", "bio12"))
d <- as.data.frame(rasterToPoints(r))
d$bio12 <- log(d$bio12)


for(k in c(5, 10, 15)){
      
      # K-MEANS
      
      d$km <- kmeans(scale(d[,c("bio1", "bio12")]), k)$cluster
      colors <- distant_colors(15)
      d$km_color <- colors[d$km]
      
      
      # SKATER
      
      # construct connectivity graph 
      gabn <- gabrielneigh(as.matrix(d[,c("x", "y")]), nnmult=4)
      nb <- graph2nb(gabn, sym=T)
      #plot(nb, coords, col="red", pch=16, cex=.5)
      
      # weight graph edges by climate dissimiliarity
      costs <- nbcosts(nb, scale(as.matrix(d[,c("bio1", "bio12")])))
      nbw <- nb2listw(nb, costs)
      
      # grow minimum spanning tree within connectivity graph
      mst <- mstree(nbw, ini=1)
      
      # partition tree into clusters
      sk <- skater(mst[,1:2], scale(as.matrix(d[,c("bio1", "bio12")])), 
                   method="euclidean", 
                   ncuts=k-1)
      
      d$sk <- sk$groups
      d$sk_color <- colors[d$km]
      
      
      
      # PLOTS
      
      eb <- element_blank()
      map <- theme(panel.background=eb, panel.grid=eb,
                   axis.ticks=eb, axis.text=eb, axis.title=eb,
                   legend.position="none")
      scatter <- theme(panel.background=eb, panel.grid=eb,
                       axis.ticks=eb, axis.text=eb, 
                       legend.position="none",
                       title=element_text(size=25))
      
      
      kmp1 <- ggplot(d, aes(x, y)) +
            geom_raster(fill=d$km) +
            map +
            coord_fixed(ratio=1.2)
      
      kmp2 <- ggplot(d, aes(bio1, bio12)) +
            geom_point(color=d$km, size=3) +
            scatter +
            labs(title="k-means")
      
      
      skp1 <- ggplot(d, aes(x, y)) +
            geom_raster(fill=d$sk) +
            map +
            coord_fixed(ratio=1.2)
      
      skp2 <- ggplot(d, aes(bio1, bio12)) +
            geom_point(color=d$sk, size=3) +
            scatter +
            labs(title="SKATER")
      
      p <- arrangeGrob(kmp2, skp2, kmp1, skp1, nrow=2)
      png(paste0("clusters_", k, ".png"), width=900, height=900)
      grid.draw(p)
      dev.off()
      
}


