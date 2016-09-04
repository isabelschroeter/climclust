


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
b <- getData("GADM", country="USA", level=1)
b <- b[b$NAME_1=="California",]
b <- spTransform(b, crs(s))

# munge
r <- mask(crop(s, extent(b)), b)
s <- mask(crop(s, extent(b)), b)
r <- subset(r, c("bio1", "bio12"))
d <- as.data.frame(rasterToPoints(r))
d$bio12 <- log(d$bio12)

climate <- scale(d[,c("bio1", "bio12")])

# continuous colors
d[,c("red", "green", "blue")] <- prcomp(climate)$x %>%
      colorwheel2d() %>%
      col2rgb() %>%
      t()



# various clustering algorithms
for(method in c("kmeans", 
                "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid",
                "skater_gabriel", "skater_30nn")){
      
      # various numbers of clusters
      for(k in c(2:5, 10, 15, 30)){
            
            # CLUSTERING
            
            if(method=="kmeans"){
                  cluster <- kmeans(climate, k)$cluster
            } else if(method=="skater_gabriel"){
                  
                  # construct connectivity graph 
                  gabn <- gabrielneigh(as.matrix(d[,c("x", "y")]), nnmult=4)
                  nb <- graph2nb(gabn, sym=T)
                  #plot(nb, coords, col="red", pch=16, cex=.5)
                  
                  # weight graph edges by climate dissimiliarity
                  costs <- nbcosts(nb, climate)
                  nbw <- nb2listw(nb, costs)
                  
                  # grow minimum spanning tree within connectivity graph
                  mst <- mstree(nbw, ini=1)
                  
                  # partition tree into clusters
                  sk <- skater(mst[,1:2], climate, 
                               method="euclidean", 
                               ncuts=k-1)
                  cluster <- sk$groups
            } else if(method=="skater_30nn"){
                  
                  # construct connectivity graph 
                  knn <- knearneigh(as.matrix(d[,c("x", "y")]), k=30)
                  nb <- knn2nb(knn)
                  #gabn <- gabrielneigh(as.matrix(d[,c("x", "y")]), nnmult=4)
                  #nb <- graph2nb(gabn, sym=T)
                  #plot(nb, coords, col="red", pch=16, cex=.5)
                  
                  # weight graph edges by climate dissimiliarity
                  costs <- nbcosts(nb, climate)
                  nbw <- nb2listw(nb, costs)
                  
                  # grow minimum spanning tree within connectivity graph
                  mst <- mstree(nbw, ini=1)
                  
                  # partition tree into clusters
                  sk <- skater(mst[,1:2], climate, 
                               method="euclidean", 
                               ncuts=k-1)
                  cluster <- sk$groups
            } else {
                  cluster <- cutree(hclust(dist(scale(d[,c("bio1", "bio12")])), 
                                           method=method), 
                                    k=k)
                  
            }
            
            
            
            # PLOTS
            
            # cluster colors
            d <- ungroup(d) %>%
                  mutate(clust=cluster) %>%
                  group_by(clust) %>%
                  mutate(r=mean(red), g=mean(green), b=mean(blue)) %>%
                  mutate(color=rgb(r, g, b, maxColorValue=255))
            
            # styles
            eb <- element_blank()
            map <- theme(panel.background=eb, panel.grid=eb,
                         axis.ticks=eb, axis.text=eb, axis.title=eb,
                         legend.position="none")
            scatter <- theme(panel.background=eb, panel.grid=eb,
                             axis.ticks=eb, axis.text=eb, 
                             legend.position="none",
                             title=element_text(size=25))
            
            # builds
            p1 <- ggplot(d, aes(x, y)) +
                  geom_raster(fill=d$color) +
                  map +
                  coord_fixed(ratio=1.2)
            
            p2 <- ggplot(d, aes(bio1, bio12)) +
                  geom_point(color=d$color, size=3) +
                  scatter +
                  labs(title=paste(method, k))
            
            # export
            p <- arrangeGrob(p1, p2, nrow=1)
            png(paste0("clusters_", method, "_", k, ".png"), width=900, height=500)
            grid.draw(p)
            dev.off()
            
      }
}

