
library(raster)
library(dplyr)
library(tidyr)
library(fastcluster)
library(FNN)
library(colormap)
library(ggplot2)
library(rgdal)

setwd("C:/Lab_projects/2016_climate_classification/climclust/webinar_figures")


##### data setup ######


# load bcm layers
files <- list.files("C:/Lab_projects/2016_Phylomodelling/Data/Climate/BCM_normals/Normals_30years", 
                    pattern="HST", full.names=T)
r <- lapply(files[2:5], readRDS) %>%
        do.call("stack", .)

# log-transform ppt
r[[4]] <- log(r[[4]])

# convert raster to matrix
v <- na.omit(cbind(coordinates(r), scale(values(r))))
colnames(v) <- c("x", "y", "cwd", "djf", "jja", "ppt")

# subsample pixels for speed
px <- sample(nrow(v), 200000) # change this to 100k for production run

# find sampled pixel most similar to each non-sampled pixel
nn <- get.knnx(v[px,3:6], v[,3:6], k=1)

# pca for colorspace
pc <- prcomp(v[,3:6])$x[,1:3]
col3d <- colors3d(pc) %>%
        col2rgb() %>%
        t()

# continuous color plot
p <- ggplot(as.data.frame(v), aes(x, y)) + 
        geom_raster(fill=rgb(col3d, maxColorValue=255)) +
        ggmap::theme_nothing() +
        coord_fixed()
png(paste0("continuous.png"), width=5, height=6, units="in", res=1000)
plot(p)
dev.off()

# perservation ranch bounary shapefile
pr <- readOGR("preservation_ranch", "PreservationRanch_boundary")
prd <- broom::tidy(pr)

# coastal conservancy shapefile
cc <- readOGR("Acquisitions", "projects_scc_2016_07_13_10_40_28")
ccd <- broom::tidy(cc)


######   build state-level hclust tree   ##########

tree <- hclust.vector(v[px,3:6], method="ward")



###### figure 1 #######

# histogram of percent land area per type, for state vs coastal conservancy, at k=20

clust <- cutree(tree, 20)
cluster <- clust[nn$nn.index]

ccr <- rasterize(cc, r[[1]]) %>% reclassify(c(-1, Inf, 1))
kr <- r[[1]]
kr[!is.na(values(kr))] <- cluster
kr <- stack(kr, ccr)
names(kr) <- c("cluster", "conservancy")

cd <- as.data.frame(values(kr)) %>%
        filter(!is.na(cluster))
ccdd <- filter(cd, !is.na(conservancy))
cd$conservancy <- 0
cd <- rbind(cd, ccdd)


cdh <- group_by(cd, conservancy, cluster) %>%
        summarize(n=n()) %>%
        group_by(conservancy) %>%
        mutate(p=n/sum(n))

cdo <- cd %>%
        group_by(cluster) %>%
        summarize(jja=mean(jja)) %>%
        arrange(jja)

cdh$cluster <- factor(cdh$cluster, levels=cdo$cluster)
cdh <- arrange(cdh, cluster, conservancy)



p <- ggplot(cdh, aes(cluster, p, group=conservancy, fill=factor(conservancy, labels=c("state", "conservancy")))) +
        geom_bar(stat="identity", position="dodge", width=.75) +
        scale_fill_manual(values=c("black", "red")) +
        theme_minimal() +
        labs(y="proportion of of total land within domain",
             fill="domain") +
        theme(legend.position=c(.5,.9))
ggsave("histogram.png", p, width=9, height=6, units="in")






###### figure 2 #######

# statewide and preservation ranch cluster maps

for(k in c(20, 50, 100, 1000)){
        
        clust <- cutree(tree, k)
        cluster <- clust[nn$nn.index]
        kr <- r[[1]]
        kr[!is.na(values(kr))] <- cluster
        
        palette <- distant_colors(k)
        clrs <- palette[cluster]
        hclrs <- as.data.frame(cbind(cluster, col3d)) %>%
                group_by(cluster) %>%
                mutate_each(funs(mean)) %>%
                mutate(hex=rgb(red, green, blue, maxColorValue=255))
        
        
        kd <-  kr %>% 
                rasterToPoints() %>%
                as.data.frame() %>%
                mutate(color=palette[layer.1])
        p <- ggplot(kd, aes(x, y)) +
                geom_raster(fill=kd$color) +
                geom_polygon(data=prd, aes(long, lat, group=group), fill=NA, color="black") +
                ggmap::theme_nothing() +
                coord_fixed()
        ggsave(paste0("statewide_", k, ".png"), p, width=6, height=9, units="in")
        
        prkd <- crop(kr, pr) %>% 
                rasterToPoints() %>%
                as.data.frame() %>%
                mutate(color=palette[layer.1])
        p <- ggplot(prkd, aes(x, y)) +
                geom_raster(fill=prkd$color) +
                geom_polygon(data=prd, aes(long, lat, group=group), fill=NA, color="black") +
                ggmap::theme_nothing() +
                coord_fixed()
        ggsave(paste0("pr_", k, ".png"), p, width=6, height=6, units="in")
        
}




