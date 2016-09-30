
library(raster)
library(dplyr)
library(tidyr)
library(fastcluster)
library(FNN)
library(colormap)
library(ggplot2)

setwd("C:/Lab_projects/2016_climate_classification/climclust")

# load bcm layers
files <- list.files("C:/Lab_projects/2016_Phylomodelling/Data/Climate/BCM_normals/Normals_30years", 
                    pattern="HST", full.names=T)
r <- lapply(files[2:5], readRDS) %>%
        do.call("stack", .)

# log-transform ppt
r[[4]] <- log(r[[4]])

# convert raster to matrix
v <- na.omit(cbind(coordinates(r), scale(values(r))))

# subsample pixels for speed
px <- sample(nrow(v), 10000) # change this to 100k for production run

# find sampled pixel most similar to each non-sampled pixel
nn <- get.knnx(v[px,3:6], v[,3:6], k=1)

# pca for colorspace
pc <- prcomp(v[,3:6])$x[,1:3]
col3d <- colors3d(pc) %>%
        col2rgb() %>%
        t()

# clustering
tree <- hclust.vector(v[px,3:6], method="ward")

for(k in c(5, 10, 15, 20)){
        
        # cut tree into specified number of clusters
        clust <- cutree(tree, k)
        
        # transfer cluster identities to non-sampled pixels
        cluster <- clust[nn$nn.index]
        
        # back-convert to raster format and export
        kr <- r[[1]]
        kr[!is.na(values(kr))] <- cluster
        writeRaster(kr, paste0("hclust/clusters_", k, ".tiff"))
        
        
        #### visualize w distant colors
        
        # get some color values
        clrs <- distant_colors(k)[cluster]
        
        # build plot
        p <- ggplot(as.data.frame(v), aes(x, y)) + 
                geom_raster(fill=clrs) +
                ggmap::theme_nothing() +
                coord_fixed()
        
        # save plot
        png(paste0("hclust/clusters_", k, ".png"), width=5, height=6, units="in", res=1000)
        plot(p)
        dev.off()
        
        # save colors
        saveRDS(clrs, paste0("hclust/clusters_", k, "_colors.rds"))
        
        
        
        #### visualize w hierarchical colors
        
        # get hierarchical color values
        hclrs <- as.data.frame(cbind(cluster, col3d)) %>%
                group_by(cluster) %>%
                mutate_each(funs(mean)) %>%
                mutate(hex=rgb(red, green, blue, maxColorValue=255))
        
        # build plot
        p <- ggplot(as.data.frame(v), aes(x, y)) + 
                geom_raster(fill=hclrs$hex) +
                ggmap::theme_nothing() +
                coord_fixed()
        
        # save plot
        png(paste0("hclust/clusters_", k, "_hcolors.png"), width=5, height=6, units="in", res=1000)
        plot(p)
        dev.off()
        
        # save colors
        saveRDS(hclrs, paste0("hclust/clusters_", k, "_hcolors.rds"))
}



