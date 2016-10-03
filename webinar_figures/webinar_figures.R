
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
names(r) <- c("cwd", "djf", "jja", "ppt")

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

# coastal conservancy acquisitions shapefile
cc <- readOGR("Acquisitions", "projects_scc_2016_07_13_10_40_28")
ccd <- broom::tidy(cc)

# coastal jusrisdiction shapefile
cj <- readOGR("SCCJurisdiction2015", "SCCJurisdiction2015_Dissolve")
cj <- spTransform(cj, crs(cc))
cj <- crop(cj, r)
cjd <- broom::tidy(cj)




######   build state-level hclust tree   ##########

tree <- hclust.vector(v[px,3:6], method="ward")



###### figure 1 #######

# histogram of percent land area per type, for state vs coastal conservancy, at k=20

# cut tree into clusters and transfer to rasters
clust <- cutree(tree, 20)
cluster <- clust[nn$nn.index]
kr <- r[[1]]
kr[!is.na(values(kr))] <- cluster

# rasterize shapefiles and stack with clusters
ccr <- rasterize(cc, r[[1]]) %>% reclassify(c(-1, Inf, 1))
cjr <- rasterize(cj, r[[1]]) %>% reclassify(c(-1, Inf, 1))
kr <- stack(kr, ccr, cjr)
names(kr) <- c("cluster", "conservancy", "coastal")
kr <- stack(kr, r)

# create conservancy vs all partitions, by double-adding conservancy lands
cd1 <- as.data.frame(rasterToPoints(kr)) %>%
        filter(!is.na(cluster))
ccdd <- filter(cd1, !is.na(conservancy))
cd1$conservancy <- 0
cd <- rbind(cd1, ccdd) 

cdh <- group_by(cd, conservancy, cluster) %>%
        filter(!is.na(coastal)) %>%
        summarize(n=n(), coastal=length(na.omit(coastal))) %>%
        group_by(conservancy) %>%
        mutate(p=n/sum(n))# %>%
        #filter(coastal > 0) # exclude climate types entirely outside the coastal region
coastal_types <- unique(cdh$cluster[cdh$coastal!=0])

cdo <- cd %>%
        group_by(cluster) %>%
        summarize(jja=mean(jja)) %>%
        arrange(jja)

cdh$cluster <- factor(cdh$cluster, levels=cdo$cluster)
cdh <- arrange(cdh, cluster, conservancy)

# expand
cdh <- expand.grid(cluster=unique(cdh$cluster),
            conservancy=unique(cdh$conservancy)) %>%
        left_join(cdh)
cdh$cluster <- factor(cdh$cluster, levels=cdo$cluster)

# reference map
p <- ggplot() +
        geom_raster(data=as.data.frame(rasterToPoints(r[[1]])),
                    aes(x,y), fill="gray85") +
        geom_polygon(data=cjd, aes(long, lat, group=group), 
                     fill="darkseagreen", color=NA) +
        geom_polygon(data=ccd, aes(long, lat, group=group), 
                     fill="darkgreen", color="darkgreen") +
        ggmap::theme_nothing() +
        coord_fixed() +
        xlim(extent(r)[c(1,2)]) +
        ylim(extent(r)[c(3,4)]) +
        annotate(geom="text", label=c("SCC Jurisdiction", "SCC Acquisitions"),
                 x=150000, y=c(250000, 200000), color=c("darkseagreen", "darkgreen"),
                 size=6, hjust=0, fontface="bold")
ggsave("reference_map.png", p, width=6, height=9, units="in")


# histogram
p <- ggplot(filter(cdh, cluster %in% coastal_types), aes(cluster, p, group=conservancy, 
                     fill=factor(conservancy, labels=c("state", "conservancy")))) +
        geom_bar(stat="identity", position="dodge", width=.9) +
        scale_fill_manual(values=c("gray", "darkgreen")) +
        theme_minimal() +
        scale_y_continuous(breaks=seq(0, 1, .1)) +
        labs(y="proportion of of total land within domain",
             fill="domain",
             x="climate type (coastal types only, sorted by ascending JJA)") +
        theme(legend.position=c(.5,.9))
ggsave("histogram.png", p, width=9, height=6, units="in")
ggsave("histogram_tall.png", p, width=9, height=16, units="in")


# coastal cluster map
#clrs <- distant_colors(length(unique(cdh$cluster)))
clrs <- distant_colors(length(unique(cd$cluster)))
eb <- element_blank()
p <- ggplot(cd) +
        geom_raster(aes(x, y, fill=factor(cluster, levels=cdo$cluster))) +
        geom_polygon(data=cjd, aes(long, lat, group=group), 
                     fill=NA, color="black") +
        theme(panel.background=eb, panel.grid=eb,
              axis.text=eb, axis.title=eb, axis.ticks=eb) +
        scale_fill_manual(values=clrs) +
        labs(fill="climate\ntype")
ggsave("coastal_cluster_map.png", p, width=6, height=6, units="in")


# histogram colored to match map
p <- ggplot() + 
        geom_bar(data=cdh, 
                 aes(cluster, p, fill=cluster,
                     group=conservancy),
                 stat="identity", position="dodge", width=.9,
                 color=NA) +
        geom_bar(data=cdh, 
                 aes(cluster, p,
                     alpha=factor(conservancy, labels=c("state", "conservancy")), 
                     group=conservancy),
                 stat="identity", position="dodge", width=.9, 
                 fill="black", color=NA) +
        scale_fill_manual(values=clrs[unique(cd$cluster) %in% coastal_types], 
                          guide=F) +
        scale_alpha_manual(values=c(0, 1)) +
        theme_minimal() +
        scale_y_continuous(breaks=seq(0, 1, .1)) +
        labs(y="proportion of of total land within domain",
             alpha="domain",
             x="climate type (coastal types only, sorted by ascending JJA)") +
        theme(legend.position=c(.5,.9))
ggsave("histogram_colored.png", p, width=9, height=6, units="in")
ggsave("histogram_colored_tall.png", p, width=9, height=16, units="in")




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




