
library(raster)
library(dplyr)
library(tidyr)
library(fastcluster)
library(FNN)
library(colormap)
library(ggplot2)

setwd("C:/Lab_projects/2016_climate_classification")

# load bcm layers
files <- list.files("C:/Lab_projects/2016_Phylomodelling/Data/Climate/BCM_normals/Normals_30years", 
                    pattern="HST", full.names=T)
r <- lapply(files[2:5], readRDS) %>%
        do.call("stack", .)

# log-transform ppt
r[[4]] <- log(r[[4]])

# convert raster to matrix
v <- scale(na.omit(cbind(coordinates(r), values(r))))

# subsample pixels for speed
px <- sample(nrow(v), 100000)

# find sampled pixel most similar to each non-sampled pixel
nn <- get.knnx(v[px,3:6], v[,3:6], k=1)

# clustering
tree <- hclust.vector(v[px,3:6], method="ward")

for(k in c(5, 10, 15, 20)){
        
        # cut tree into specified number of clusters
        clust <- cutree(tree, k)
        
        # transfer cluster identities to non-sampled pixels
        cluster <- clust[nn$nn.index]
        
        # get some color values
        clrs <- distant_colors(k)[cluster]
        
        # build plot
        p <- ggplot(as.data.frame(v), aes(x, y)) + 
                geom_raster(fill=clrs) +
                ggmap::theme_nothing() +
                coord_fixed()
        
        # save plot
        png(paste0("clusters_", k, ".png"), width=5, height=6, units="in", res=1000)
        plot(p)
        dev.off()
}



