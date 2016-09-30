
# bin bcm data


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
r <- lapply(files[c(2,5,6,7)], readRDS) %>%
        do.call("stack", .)

# mean of tmin and tmax
temp <- mean(r[[3]], r[[4]])

# log-transform ppt
ppt <- log(r[[2]])


# convert raster to matrix
r <- stack(temp, ppt, r[[1]])
names(r) <- c("tmean", "ppt", "cwd")
v <- na.omit(cbind(coordinates(r), values(r)))


f <- as.data.frame(v)

cutoffs <- function(x, n, stat){
        if(stat=="reg") return(seq(min(x, na.rm=T), max(x, na.rm=T), length.out=n+1))
        if(stat=="rank") quantile(x, prob=seq(0, 1, length.out=n+1))
}



# regularly spaced bins

reg <- function(x, n) as.integer(cut(x, seq(min(x, na.rm=T), max(x, na.rm=T), length.out=n+1)))


l <- apply(as.matrix(f[3:5]), 2, cutoffs, n=3, stat="reg") %>%
        as.data.frame()


b <- apply(v[,3:5], 2, reg, n=3)

colors <- colors3d(b, order=4, inversion=4)

p <- ggplot(f, aes(x, y)) +
        geom_raster(fill=colors) +
        ggmap::theme_nothing() +
        coord_fixed()

# save plot
png(paste0("bins_regular.png"), width=5, height=6, units="in", res=1000)
plot(p)
dev.off()


# rank-spaced bins

rank <- function(x, n) as.integer(cut(x, quantile(x, prob=seq(0, 1, length.out=n+1))))

b <- apply(v[,3:5], 2, rank, n=3)


l <- apply(as.matrix(f[3:5]), 2, cutoffs, n=3, stat="rank") %>%
        as.data.frame()

colors <- colors3d(b, order=4, inversion=4)

p <- ggplot(f, aes(x, y)) +
        geom_raster(fill=colors) +
        ggmap::theme_nothing() +
        coord_fixed()

# save plot
png(paste0("bins_rank.png"), width=5, height=6, units="in", res=1000)
plot(p)
dev.off()
