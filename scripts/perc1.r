############################################## #
## R-Script Perception
##                  xxx
## =========================================== =
## Project: CAA Hamburg 2016
## Authors: Oliver Nakoinz & Franziska Faupel
## Version: 01
## Date of last changes: 24.01.2016
## Data: SHKR
## Author of data: SHKR
## Purpose: presentation
## Content: 1. xxxxx, 2. xxx
##          3. xxx, 4. xxx
## Description: applies some basic techniques of
##             boundary analysis
## Licence data: -
## Licence Script: GPL 
##    (http://www.gnu.org/licenses/gpl-3.0.html)
############################################## #

# 0. Preparation ===============================
#  working direktory
wd <- "/home/fon/daten/analyse/perception_IA/proj_perIA1"
setwd(wd)
# load("7ws/ws01.rws")

# 1. Functions ===============================
# cost functions
# Tobler1993 velocity
tobler1993a <- function(s){6 * exp(-3.5 * abs(s + 0.05))}      # km/h
tobler1993b <- function(s){0.36 * exp(-3.5 * abs(s + 0.05))}   # m/min
# Minetti2002 Metabolische Kosten in J/(kg*m) für laufen
minetti2002w <- function(s){(280.5 * s^5 - 58.7 * s^4 - 76.8 * s^3 + 51.9 * s^2 + 19.6  * s + 2.5)}
minetti2002wi <- function(s){1/(280.5 * s^5 - 58.7 * s^4 - 76.8 * s^3 + 51.9 * s^2 + 19.6 * s + 2.5)}
# Herzog2012 (after data from Minetti) walking. metabolic costs J/(kg*m)
herzog2012_w <- function(s){(1337.8 * s^6 + 278.19 * s^5 - 517.39 * s^4 - 78.199 * s^3  + 93.419 * s^2 + 19.825 * s + 1.64)}
herzog2012_wi <- function(s){1/(1337.8 * s^6 + 278.19 * s^5 - 517.39 * s^4 - 78.199 *   s^3 + 93.419 * s^2 + 19.825 * s + 1.64)}

pdf("6pictures/c9_cfunc.pdf", height=3, width=6, bg = "white") 
par( mai = c(1, 1, 0.1, 0.1))
matplot(seq(-0.5,0.5,0.01), herzog2012_wi(seq(-0.5,0.5,0.01)), type="l", xlab=expression(paste("slope (", Delta, "h/", Delta, "d)")), ylab="metabolic costs (kg m/J)")
lines(seq(-0.5,0.5,0.01), minetti2002wi(seq(-0.5,0.5,0.01)), type="l", lty=2)
dev.off() 

# auxilliary function
hdiff <- function(x){x[2]-x[1]}



trans.pol <- function(a, b=c(0,0)){
    x   <- a[1]
    y   <- a[2]
    xt  <- b[1]
    yt  <- b[2]
    r <- (((x-xt)^2)+((y-yt)^2))^0.5
    if ((x-xt) >= 0 & (y-yt) >= 0)  phi <- atan((y-yt)/(x-xt))
    if ((x-xt) < 0  & (y-yt) >= 0)  phi <- atan((y-yt)/(x-xt))  + pi
    if ((x-xt) < 0  & (y-yt) < 0)   phi <- atan((y-yt)/(x-xt))  - pi
    if ((x-xt) >= 0 & (y-yt) < 0)   phi <- atan((y-yt)/(x-xt))  + 2 * pi
    return(c(r, phi))
}

trans.cartes <- function(a, b=c(0,0)){
    r   <- a[1]
    phi <- a[2]
    xt  <- b[1]
    yt  <- b[2]
    x   <- r*cos(phi) + xt
    y   <- r*sin(phi) + yt
    return(c(x, y))
}

# 2. Data ===============================
file_fs    <- "fs_ag" 
file_v    <- "viewpoints" 
file_e    <- "ecken" 
file_h    <- "huegel" 
file_cdist    <- "2data/cdist38.csv" 
#file_srtm_hs <- "2geodata/srtm_gk3_hs.asc" 
file_srtm_a <- "3geodata/srtm_gk3_a2.tif" 
file_grid_v <- "3geodata/srtm_gk3_v.asc" 

crs1 <- "+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel +towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7 +units=m +no_defs"
#crs2 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

install.packages('sp') 
library(sp)      
#library(proj4)           
library(rgdal)         
#library(spatstat)     
#library(gdata) 

avs <- paste(wd,"/2data",sep="")
spdf_fs <- readOGR(avs, p4s=NULL, file_fs) 
spdf_fs@data <- cbind(x=coordinates(spdf_fs)[,1], y=coordinates(spdf_fs)[,2],u=0,v=0,spdf_fs@data)
spdf_v <- readOGR(avs, p4s=NULL, file_v) 
spdf_e <- readOGR(avs, p4s=NULL, file_e) 
spdf_h <- readOGR(avs, p4s=NULL, file_h) 

spdf_cdist <- read.table(file_cdist, sep=';', header=TRUE)
spdf_cdist <- cbind(spdf_cdist, sdist=0, u=0, v=0)
spdf_cdist <- rbind(spdf_cdist, c(356, 3530750, 5328900,0.0,0,0,0))

sgdf_srtm <- readGDAL(file_srtm_a) 
sgdf_grid_v <- read.asciigrid(file_grid_v) 

proj4string(sgdf_srtm)  <- CRS(as.character(crs1)) 
proj4string(sgdf_grid_v) <- CRS(as.character(crs1)) 

library(raster)
srtm_slope  <- terrain(raster(sgdf_srtm), opt='slope')
srtm_aspect <- terrain(raster(sgdf_srtm), opt='aspect')
srtm_shade  <- hillShade(srtm_slope, srtm_aspect, 40, 270)
#top.colors = colorRampPalette(c("#81BB7C", "#E5CE98", "#B89E83", "#9E5D4C"), bias=1.2) 
top.colors = colorRampPalette(c("#3BA632", "#E5CE98",  "#5C382E"), bias=1.2) 

plot(raster(sgdf_srtm), col = top.colors(20))
plot(srtm_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)

# unveränderliche Punkte
pt_stat <- cbind(x1=c(coordinates(spdf_e)[,1], 3530750), y1=c(coordinates(spdf_e)[,2], 5328900), x2=c(coordinates(spdf_e)[,1], 3530750), y2=c(coordinates(spdf_e)[,2], 5328900))

# 2. cultural distances ===============================
for (i in 1:length(spdf_cdist[,1]))   {
    x1 <- spdf_cdist[i,2]
    y1 <- spdf_cdist[i,3]
    x2 <- 3530750
    y2 <- 5328900
   spdf_cdist[i,5] <- sqrt( (x1-x2)^2 + (y1-y2)^2   ) 
}
spdf_cdist[,4] <- spdf_cdist[,4] / max(spdf_cdist[,4])
coordinates(spdf_cdist )=~x+y  


# 3. cost distance ===============================
library(gdistance)
vp <- c(3530750, 5328900)
nz= 8
# transitional object 
hd <- transition(raster(sgdf_srtm), hdiff, nz, symm=TRUE)    
slope <- geoCorrection(hd,scl=FALSE)    
adj <- adjacent(x=raster(sgdf_srtm), cells=1:ncell(raster(sgdf_srtm)), direction=nz) 
cost <- slope       
cost[adj] <- herzog2012_wi(slope[adj]) 
conduct <- geoCorrection(cost, scl=FALSE) 

cost_surface1 <- accCost(conduct, vp)


image(cost_surface1)
#cl <- contour(cost_surface1, add=T, method = "edge", levels = c(40000), drawlabels = F)    

sgdf_cost <- as(cost_surface1, 'SpatialGridDataFrame')
im <- as.image.SpatialGridDataFrame(sgdf_cost )
cl <- contourLines(im, levels=40000)
library(maptools)
SLDF <- ContourLines2SLDF(cl)
writeOGR(SLDF, ".", "costdist_cl", driver="ESRI Shapefile")

spdf_cd_pt <- readOGR(avs, p4s=NULL, "cd_pt") 

pdf("5pictures/cd_line.pdf", height=4.3, width=6) 
par(mar = c(0, 0, 0, 0))
    plot(raster(sgdf_srtm), col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE)
    plot(srtm_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
    #points(spdf_h, cex=0.3, col=grey, pch=16)
    points(spdf_fs, cex=1, col="#FF004E", pch=16)
    lines(SLDF)
dev.off() 

png("5pictures/cd_lin.png", height=1787*4, width=2135*4, units = "px") 
par(mar = c(0, 0, 0, 0))
    plot(raster(sgdf_srtm), col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE, margin=F)
    plot(srtm_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
    #points(spdf_h, cex=0.3, col=grey, pch=16)
    points(spdf_fs, cex=20, col="#FF004E", pch=16)
    lines(SLDF,lwd=20)
dev.off() 












# 4. model 1: huegel ====================================
# Start
vp <- c(3530750, 5328900)








# 5. model 2: kultur ====================================
pt_dyn1 <- coordinates(spdf_cdist)

vp <- c(3530750, 5328900)

# 6.1 transform to polar coordinates ===========
tp2 <- function(a,b) trans.pol(a,vp)
pt_dyn2 <- apply(pt_dyn1, 1, tp2)

# 6.2 manipulate distances and angles ===========
pt_d <- t(pt_dyn2)[,1]

cdist <- spdf_cdist@data[,2] *0.5 + 0.5
pt_d[] <- cdist * spdf_cdist@data[,3] 
pt_w <- t(pt_dyn2)[,2]
pt_w[355] <- 0
pt_dyn3 <- cbind(pt_d, pt_w)

# 6.3 transform back ===========
tc2 <- function(a,b)  trans.cartes (a,vp)
pt_dyn4 <- t(apply(pt_dyn3, 1, tc2))
pt_dyn1 
pt_dyn <- cbind(x1=pt_dyn1[,1], y1=pt_dyn1[,2], x2=pt_dyn4[,1], y2=pt_dyn4[,2])
pt <- rbind (pt_dyn)

# produce a spatialpointsdataframe
#pt_cd1 <- SpatialPoints(pt[,1:2], proj4string = CRS(as.character(crs1)))
#pt_cd2 <- SpatialPoints(pt[,3:4], proj4string = CRS(as.character(crs1)))
pt_cdist1 <- SpatialPointsDataFrame(pt[,1:2], data.frame(pt[,1:2]), proj4string = CRS(as.character(crs1)))
pt_cdist2 <- SpatialPointsDataFrame(pt[,3:4], data.frame(pt[,3:4]), proj4string = CRS(as.character(crs1)))

pdf("5pictures/cdist_pt.pdf", height=4.3, width=6) 
par(mar = c(0, 0, 0, 0))
plot(raster(sgdf_srtm), col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE)
plot(srtm_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
plot(pt_cdist1, add=T)
points(pt_cdist2)
dev.off() 

png("5pictures/cdist_pt.png", height=1787*4, width=2135*4, units = "px") 
par(mar = c(0, 0, 0, 0))
plot(raster(sgdf_srtm), col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE, margin=F)
plot(srtm_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
plot(pt_cdist1, add=T, lwd= 20 ,cex=15)
points(pt_cdist2, lwd= 20, cex=15)
dev.off() 


# 6.4 rubber shed ===========
sgdf <- sgdf_srtm 
bb <- sgdf@bbox
sgdf@data[,1]  <- 0
srtm <- data.frame(x=coordinates(sgdf_srtm)[,1], y=coordinates(sgdf_srtm)[,2], z=sgdf_srtm@data[,1], u=0, v=0)
library(spatstat)
cm_pp <- ppp(pt[,1], pt[,2], window=owin(c(bb[1,1],bb[1,2]),c(bb[2,1],bb[2,2])))
cm_del <- delaunay(cm_pp)
cm_pp2 <- ppp(pt[,3], pt[,4], window=owin(c(bb[1,1],bb[1,2]),c(bb[2,1],bb[2,2])))
cm_del2 <- delaunay(cm_pp2)
trans.x <- function(x) u <- b0*x[1] +b1*x[2] + b2
trans.y <- function(x) u <- b3*x[1] +b4*x[2] + b5
for (i in 1:cm_del$n)   {
    x1 <- tiles(cm_del)[[i]]$bdry[[1]][[1]][1]
    y1 <- tiles(cm_del)[[i]]$bdry[[1]][[2]][1]
    x2 <- tiles(cm_del)[[i]]$bdry[[1]][[1]][2]
    y2 <- tiles(cm_del)[[i]]$bdry[[1]][[2]][2]
    x3 <- tiles(cm_del)[[i]]$bdry[[1]][[1]][3]
    y3 <- tiles(cm_del)[[i]]$bdry[[1]][[2]][3]
    u1 <- pt_cdist2@data[which(pt_cdist1@data[,1] == x1), 1][1]
    v1 <- pt_cdist2@data[which(pt_cdist1@data[,2] == y1), 2][1]
    u2 <- pt_cdist2@data[which(pt_cdist1@data[,1] == x2), 1][1]
    v2 <- pt_cdist2@data[which(pt_cdist1@data[,2] == y2), 2][1]
    u3 <- pt_cdist2@data[which(pt_cdist1@data[,1] == x3), 1][1]
    v3 <- pt_cdist2@data[which(pt_cdist1@data[,2] == y3), 2][1]
    cmdf <- cbind(x=c(x1,x2,x3), y=c(y1,y2,y3), u=c(u1,u2,u3), v=c(v1,v2,v3))
    n <- length(cmdf[,1])
    p <- sum((cmdf[,1] - mean(cmdf[,1]))  * (cmdf[,2] - mean(cmdf[,2])))  /  sum((cmdf[,2] - mean(cmdf[,2]))^2)
    q <- sum((cmdf[,1] - mean(cmdf[,1]))  * (cmdf[,2] - mean(cmdf[,2])))  /  sum((cmdf[,1] - mean(cmdf[,1]))^2)
    b1 <- (sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,3] - cmdf[,1])) - q  * sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,3] - cmdf[,1]))) /  (sum((cmdf[,2] - mean(cmdf[,2]))^2)  - q * sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,2] - mean(cmdf[,2]))))            
    b3 <- (sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,4] - cmdf[,2])) - p  * sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,4] - cmdf[,2]))) /  (sum((cmdf[,1] - mean(cmdf[,1]))^2)  - p * sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,1] - mean(cmdf[,1]))))    
    b4 <- 1 - p * b3 + sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,4] - cmdf[,2]))  /  sum((cmdf[,2] - mean(cmdf[,2]))^2)   
    b0 <- 1 - q * b1 + sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,3] - cmdf[,1]))  /  sum((cmdf[,1] - mean(cmdf[,1]))^2)   
    b2 <- (1/n) * sum(cmdf[,3]-cmdf[,1]) +  mean(cmdf[,1]) - b0 * mean(cmdf[,1]) - b1 * mean(cmdf[,2])
    b5 <- (1/n) * sum(cmdf[,4]-cmdf[,2]) +  mean(cmdf[,2]) - b4 * mean(cmdf[,2]) - b3 * mean(cmdf[,1])
    cm_ind <- which(inside.owin(srtm[,1], srtm[,2], tiles(cm_del)[[i]]) == T)
    srtm[cm_ind,4] <- apply(srtm[cm_ind,1:2], 1, trans.x)
    srtm[cm_ind,5] <- apply(srtm[cm_ind,1:2], 1, trans.y)
    
    cm_ind <- which(inside.owin(spdf_fs@coords[,1] , spdf_fs@coords[,2] , tiles(cm_del)[[i]]) == T)
    if(length(cm_ind) > 0){
        spdf_fs@data[cm_ind,3] <- apply(spdf_fs@data[cm_ind,1:2], 1, trans.x)
        spdf_fs@data[cm_ind,4] <- apply(spdf_fs@data[cm_ind,1:2], 1, trans.y)
    }
}
coordinates(srtm)  = ~u+v
proj4string(srtm)  <- CRS(as.character(crs1)) 
cm_raster2 <- rasterize(srtm, raster(sgdf), field='z', update=TRUE, proj4string = CRS(as.character(crs1)))
cm_raster2[which(getValues(cm_raster2) == 0)] <- NA
cm_raster2 <- focal(cm_raster2, w= matrix(rep(1,25), nrow=5),  na.rm=TRUE, mean, NAonly=TRUE) 
#cm_raster2[is.na(cm_raster2)] <- 0

sgdf_srtm_cdist <- cm_raster2
srtm_cdist_slope  <- terrain(sgdf_srtm_cdist, opt='slope')
srtm_cdist_aspect <- terrain(sgdf_srtm_cdist, opt='aspect')
srtm_cdist_shade  <- hillShade(srtm_cdist_slope, srtm_cdist_aspect, 40, 270)
spdf_fs_cdist <- data.frame(cbind(u=spdf_fs@data[,3], v=spdf_fs@data[,4]))
coordinates(spdf_fs_cdist)  = ~u+v

pdf("5pictures/mod_cdist.pdf", height=4.3, width=6) 
par(mar = c(0, 0, 0, 0))
plot(sgdf_srtm_cdist, col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE)
plot(srtm_cdist_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
points(spdf_fs_cdist, cex=1, col="#FF004E", pch=16)
dev.off() 

png("5pictures/mod_cdist.png", height=1787*4, width=2135*4, units = "px") 
par(mar = c(0, 0, 0, 0))
plot(sgdf_srtm_cdist, col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE)
plot(srtm_cdist_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
points(spdf_fs_cdist, cex=20, col="#FF004E", pch=16)
dev.off() 




#  linien
LinesList <- list()         # leere LinesList Liste erstellen
pt2 <- data.frame(pt)
pt2 <- cbind(pt2, name="xx")

for(i in seq(along=pt[,1])) {         
            m <- matrix(data = c(pt[i,1],pt[i,3],pt[i,2],pt[i,4]), nrow=2, ncol=2)  
            L <- Line(m); LL <- list(L)
            name  <- paste("edge", "_", i, sep="")
            LLL <- Lines(LL, ID = name)             
            LinesList[length(LinesList)+1] <- LLL   
            pt2[i,5] <- name         
        }

plot(sdf)

#deldf <- data.frame(deldf_i,deldf_x1,deldf_y1,deldf_k,deldf_x2,deldf_y2,deldf_name)
#sldf2 <- data.frame(pt2)
sl <- SpatialLines(LinesList, proj4string = CRS(as.character(crs1)))
sdf <- SpatialLinesDataFrame(sl, pt2, match.ID = FALSE)


pdf("5pictures/cdist_pt2.pdf", height=4.3, width=6) 
par(mar = c(0, 0, 0, 0))
plot(raster(sgdf_srtm), col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE)
plot(srtm_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
lines(sdf, col="#6D5123")
plot(pt_cdist1, add=T)
points(pt_cdist2)
points(spdf_fs_cdist, cex=1, col="#FF004E", pch=16)
dev.off() 

png("5pictures/cdist_pt2.png", height=1787*4, width=2135*4, units = "px") 
par(mar = c(0, 0, 0, 0))
plot(raster(sgdf_srtm), col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE, margin=F)
plot(srtm_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
lines(sdf, lwd= 20, col="#6D5123")
plot(pt_cdist1, add=T, lwd= 20 ,cex=15)
points(pt_cdist2, lwd= 20, cex=15)
points(spdf_fs_cdist, cex=20, col="#FF004E", pch=16)
dev.off() 


# 6. model 3: lcp ====================================
pt_dyn1 <- coordinates(spdf_cd_pt)

vp <- c(3530750, 5328900)
#cm_denspt <- data.frame(coordinates(spdf_cd_pt), 10000000 * over(cent_meg, sgdf_meg_dens), r=0, phi=0, r2=0, phi2=0, x2=0, y2=0)

# 6.1 transform to polar coordinates ===========
tp2 <- function(a,b) trans.pol(a,vp)
pt_dyn2 <- apply(pt_dyn1, 1, tp2)

# 6.2 manipulate distances and angles ===========
pt_d <- t(pt_dyn2)[,1]
pt_d[] <- 52000
pt_w <- t(pt_dyn2)[,2]
pt_dyn3 <- cbind(pt_d, pt_w)

# 6.3 transform back ===========
tc2 <- function(a,b)  trans.cartes (a,vp)
pt_dyn4 <- t(apply(pt_dyn3, 1, tc2))
pt_dyn1 
pt_dyn <- cbind(x1=pt_dyn1[,1], y1=pt_dyn1[,2], x2=pt_dyn4[,1], y2=pt_dyn4[,2])
pt <- rbind (pt_dyn, pt_stat)

# produce a spatialpointsdataframe
#pt_cd1 <- SpatialPoints(pt[,1:2], proj4string = CRS(as.character(crs1)))
#pt_cd2 <- SpatialPoints(pt[,3:4], proj4string = CRS(as.character(crs1)))
pt_cd1 <- SpatialPointsDataFrame(pt[,1:2], data.frame(pt[,1:2]), proj4string = CRS(as.character(crs1)))
pt_cd2 <- SpatialPointsDataFrame(pt[,3:4], data.frame(pt[,3:4]), proj4string = CRS(as.character(crs1)))



pdf("5pictures/cd_pt.pdf", height=4.3, width=6) 
par(mar = c(0, 0, 0, 0))
    plot(raster(sgdf_srtm), col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE)
    plot(srtm_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
    plot(pt_cd1, add=T)
    points(pt_cd2)
dev.off() 

png("5pictures/cd_pt.png", height=1787*4, width=2135*4, units = "px") 
par(mar = c(0, 0, 0, 0))
    plot(raster(sgdf_srtm), col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE, margin=F)
    plot(srtm_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
    plot(pt_cd1, add=T, lwd= 20 ,cex=15)
    points(pt_cd2, lwd= 20, cex=15)
dev.off() 


# 6.4 rubber shed ===========
sgdf <- sgdf_srtm 
bb <- sgdf@bbox
sgdf@data[,1]  <- 0
srtm <- data.frame(x=coordinates(sgdf_srtm)[,1], y=coordinates(sgdf_srtm)[,2], z=sgdf_srtm@data[,1], u=0, v=0)
library(spatstat)
cm_pp <- ppp(pt[,1], pt[,2], window=owin(c(bb[1,1],bb[1,2]),c(bb[2,1],bb[2,2])))
cm_del <- delaunay(cm_pp)
cm_pp2 <- ppp(pt[,3], pt[,4], window=owin(c(bb[1,1],bb[1,2]),c(bb[2,1],bb[2,2])))
cm_del2 <- delaunay(cm_pp2)
trans.x <- function(x) u <- b0*x[1] +b1*x[2] + b2
trans.y <- function(x) u <- b3*x[1] +b4*x[2] + b5
for (i in 1:cm_del$n)   {
    x1 <- tiles(cm_del)[[i]]$bdry[[1]][[1]][1]
    y1 <- tiles(cm_del)[[i]]$bdry[[1]][[2]][1]
    x2 <- tiles(cm_del)[[i]]$bdry[[1]][[1]][2]
    y2 <- tiles(cm_del)[[i]]$bdry[[1]][[2]][2]
    x3 <- tiles(cm_del)[[i]]$bdry[[1]][[1]][3]
    y3 <- tiles(cm_del)[[i]]$bdry[[1]][[2]][3]
    u1 <- pt_cd2@data[which(pt_cd1@data[,1] == x1), 1][1]
    v1 <- pt_cd2@data[which(pt_cd1@data[,2] == y1), 2][1]
    u2 <- pt_cd2@data[which(pt_cd1@data[,1] == x2), 1][1]
    v2 <- pt_cd2@data[which(pt_cd1@data[,2] == y2), 2][1]
    u3 <- pt_cd2@data[which(pt_cd1@data[,1] == x3), 1][1]
    v3 <- pt_cd2@data[which(pt_cd1@data[,2] == y3), 2][1]
    cmdf <- cbind(x=c(x1,x2,x3), y=c(y1,y2,y3), u=c(u1,u2,u3), v=c(v1,v2,v3))
    n <- length(cmdf[,1])
    p <- sum((cmdf[,1] - mean(cmdf[,1]))  * (cmdf[,2] - mean(cmdf[,2])))  /  sum((cmdf[,2] - mean(cmdf[,2]))^2)
    q <- sum((cmdf[,1] - mean(cmdf[,1]))  * (cmdf[,2] - mean(cmdf[,2])))  /  sum((cmdf[,1] - mean(cmdf[,1]))^2)
    b1 <- (sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,3] - cmdf[,1])) - q  * sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,3] - cmdf[,1]))) /  (sum((cmdf[,2] - mean(cmdf[,2]))^2)  - q * sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,2] - mean(cmdf[,2]))))            
    b3 <- (sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,4] - cmdf[,2])) - p  * sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,4] - cmdf[,2]))) /  (sum((cmdf[,1] - mean(cmdf[,1]))^2)  - p * sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,1] - mean(cmdf[,1]))))    
    b4 <- 1 - p * b3 + sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,4] - cmdf[,2]))  /  sum((cmdf[,2] - mean(cmdf[,2]))^2)   
    b0 <- 1 - q * b1 + sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,3] - cmdf[,1]))  /  sum((cmdf[,1] - mean(cmdf[,1]))^2)   
    b2 <- (1/n) * sum(cmdf[,3]-cmdf[,1]) +  mean(cmdf[,1]) - b0 * mean(cmdf[,1]) - b1 * mean(cmdf[,2])
    b5 <- (1/n) * sum(cmdf[,4]-cmdf[,2]) +  mean(cmdf[,2]) - b4 * mean(cmdf[,2]) - b3 * mean(cmdf[,1])
    cm_ind <- which(inside.owin(srtm[,1], srtm[,2], tiles(cm_del)[[i]]) == T)
    srtm[cm_ind,4] <- apply(srtm[cm_ind,1:2], 1, trans.x)
    srtm[cm_ind,5] <- apply(srtm[cm_ind,1:2], 1, trans.y)
    
    cm_ind <- which(inside.owin(spdf_fs@coords[,1] , spdf_fs@coords[,2] , tiles(cm_del)[[i]]) == T)
    if(length(cm_ind) > 0){
        spdf_fs@data[cm_ind,3] <- apply(spdf_fs@data[cm_ind,1:2], 1, trans.x)
        spdf_fs@data[cm_ind,4] <- apply(spdf_fs@data[cm_ind,1:2], 1, trans.y)
    }
}
coordinates(srtm)  = ~u+v
proj4string(srtm)  <- CRS(as.character(crs1)) 
cm_raster2 <- rasterize(srtm, raster(sgdf), field='z', update=TRUE, proj4string = CRS(as.character(crs1)))
cm_raster2[which(getValues(cm_raster2) == 0)] <- NA
cm_raster2 <- focal(cm_raster2, w= matrix(rep(1,25), nrow=5),  na.rm=TRUE, mean, NAonly=TRUE) 
#cm_raster2[is.na(cm_raster2)] <- 0

sgdf_srtm_cd <- cm_raster2
srtm_cd_slope  <- terrain(sgdf_srtm_cd, opt='slope')
srtm_cd_aspect <- terrain(sgdf_srtm_cd, opt='aspect')
srtm_cd_shade  <- hillShade(srtm_cd_slope, srtm_cd_aspect, 40, 270)
spdf_fs_cd <- data.frame(cbind(u=spdf_fs@data[,3], v=spdf_fs@data[,4]))
coordinates(spdf_fs_cd)  = ~u+v

pdf("5pictures/mod_cd.pdf", height=4.3, width=6) 
par(mar = c(0, 0, 0, 0))
plot(sgdf_srtm_cd, col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE)
plot(srtm_cd_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
points(spdf_fs_cd, cex=1, col="#FF004E", pch=16)
dev.off() 

png("5pictures/mod_cd.png", height=1787*4, width=2135*4, units = "px") 
par(mar = c(0, 0, 0, 0))
plot(sgdf_srtm_cd, col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE)
plot(srtm_cd_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
points(spdf_fs_cd, cex=20, col="#FF004E", pch=16)
dev.off() 


# 7. model 4: sicht  ====================================
pt_dyn1 <- coordinates(spdf_v)

vp <- c(3530750, 5328900)
#cm_denspt <- data.frame(coordinates(spdf_cd_pt), 10000000 * over(cent_meg, sgdf_meg_dens), r=0, phi=0, r2=0, phi2=0, x2=0, y2=0)

# 7.1 transform to polar coordinates ===========
tp2 <- function(a,b) trans.pol(a,vp)
pt_dyn2 <- apply(pt_dyn1, 1, tp2)

# 7.2 manipulate distances and angles ===========
pt_d <- t(pt_dyn2)[,1]
pt_d[] <- 30000
pt_w <- t(pt_dyn2)[,2]
pt_dyn3 <- cbind(pt_d, pt_w)

# 7.3 transform back ===========
tc2 <- function(a,b)  trans.cartes (a,vp)
pt_dyn4 <- t(apply(pt_dyn3, 1, tc2))
pt_dyn1 
pt_dyn <- cbind(x1=pt_dyn1[,1], y1=pt_dyn1[,2], x2=pt_dyn4[,1], y2=pt_dyn4[,2])
pt <- rbind (pt_dyn, pt_stat)

# produce a spatialpointsdataframe
#pt_cd1 <- SpatialPoints(pt[,1:2], proj4string = CRS(as.character(crs1)))
#pt_cd2 <- SpatialPoints(pt[,3:4], proj4string = CRS(as.character(crs1)))
pt_v1 <- SpatialPointsDataFrame(pt[,1:2], data.frame(pt[,1:2]), proj4string = CRS(as.character(crs1)))
pt_v2 <- SpatialPointsDataFrame(pt[,3:4], data.frame(pt[,3:4]), proj4string = CRS(as.character(crs1)))


pdf("5pictures/v_pt.pdf", height=4.3, width=6) 
par(mar = c(0, 0, 0, 0))
plot(raster(sgdf_srtm), col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE)
plot(srtm_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
plot(pt_v1, add=T)
points(pt_v2)
dev.off() 

png("5pictures/v_pt.png", height=1787*4, width=2135*4, units = "px") 
par(mar = c(0, 0, 0, 0))
plot(raster(sgdf_srtm), col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE, margin=F)
plot(srtm_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
plot(pt_cd1, add=T, lwd= 20 ,cex=15)
points(pt_cd2, lwd= 20, cex=15)
dev.off() 

# 7.4 rubber shed ===========
sgdf <- sgdf_srtm 
bb <- sgdf@bbox
sgdf@data[,1]  <- 0
srtm <- data.frame(x=coordinates(sgdf_srtm)[,1], y=coordinates(sgdf_srtm)[,2], z=sgdf_srtm@data[,1], u=0, v=0)
library(spatstat)
cm_pp <- ppp(pt[,1], pt[,2], window=owin(c(bb[1,1],bb[1,2]),c(bb[2,1],bb[2,2])))
cm_del <- delaunay(cm_pp)
cm_pp2 <- ppp(pt[,3], pt[,4], window=owin(c(bb[1,1],bb[1,2]),c(bb[2,1],bb[2,2])))
cm_del2 <- delaunay(cm_pp2)
trans.x <- function(x) u <- b0*x[1] +b1*x[2] + b2
trans.y <- function(x) u <- b3*x[1] +b4*x[2] + b5
for (i in 1:cm_del$n)   {
    x1 <- tiles(cm_del)[[i]]$bdry[[1]][[1]][1]
    y1 <- tiles(cm_del)[[i]]$bdry[[1]][[2]][1]
    x2 <- tiles(cm_del)[[i]]$bdry[[1]][[1]][2]
    y2 <- tiles(cm_del)[[i]]$bdry[[1]][[2]][2]
    x3 <- tiles(cm_del)[[i]]$bdry[[1]][[1]][3]
    y3 <- tiles(cm_del)[[i]]$bdry[[1]][[2]][3]
    u1 <- pt_v2@data[which(pt_v1@data[,1] == x1), 1][1]
    v1 <- pt_v2@data[which(pt_v1@data[,2] == y1), 2][1]
    u2 <- pt_v2@data[which(pt_v1@data[,1] == x2), 1][1]
    v2 <- pt_v2@data[which(pt_v1@data[,2] == y2), 2][1]
    u3 <- pt_v2@data[which(pt_v1@data[,1] == x3), 1][1]
    v3 <- pt_v2@data[which(pt_v1@data[,2] == y3), 2][1]
    cmdf <- cbind(x=c(x1,x2,x3), y=c(y1,y2,y3), u=c(u1,u2,u3), v=c(v1,v2,v3))
    n <- length(cmdf[,1])
    p <- sum((cmdf[,1] - mean(cmdf[,1]))  * (cmdf[,2] - mean(cmdf[,2])))  /  sum((cmdf[,2] - mean(cmdf[,2]))^2)
    q <- sum((cmdf[,1] - mean(cmdf[,1]))  * (cmdf[,2] - mean(cmdf[,2])))  /  sum((cmdf[,1] - mean(cmdf[,1]))^2)
    b1 <- (sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,3] - cmdf[,1])) - q  * sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,3] - cmdf[,1]))) /  (sum((cmdf[,2] - mean(cmdf[,2]))^2)  - q * sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,2] - mean(cmdf[,2]))))            
    b3 <- (sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,4] - cmdf[,2])) - p  * sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,4] - cmdf[,2]))) /  (sum((cmdf[,1] - mean(cmdf[,1]))^2)  - p * sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,1] - mean(cmdf[,1]))))    
    b4 <- 1 - p * b3 + sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,4] - cmdf[,2]))  /  sum((cmdf[,2] - mean(cmdf[,2]))^2)   
    b0 <- 1 - q * b1 + sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,3] - cmdf[,1]))  /  sum((cmdf[,1] - mean(cmdf[,1]))^2)   
    b2 <- (1/n) * sum(cmdf[,3]-cmdf[,1]) +  mean(cmdf[,1]) - b0 * mean(cmdf[,1]) - b1 * mean(cmdf[,2])
    b5 <- (1/n) * sum(cmdf[,4]-cmdf[,2]) +  mean(cmdf[,2]) - b4 * mean(cmdf[,2]) - b3 * mean(cmdf[,1])
    cm_ind <- which(inside.owin(srtm[,1], srtm[,2], tiles(cm_del)[[i]]) == T)
    srtm[cm_ind,4] <- apply(srtm[cm_ind,1:2], 1, trans.x)
    srtm[cm_ind,5] <- apply(srtm[cm_ind,1:2], 1, trans.y)
    
    cm_ind <- which(inside.owin(spdf_fs@coords[,1] , spdf_fs@coords[,2] , tiles(cm_del)[[i]]) == T)
    if(length(cm_ind) > 0){
        spdf_fs@data[cm_ind,3] <- apply(spdf_fs@data[cm_ind,1:2], 1, trans.x)
        spdf_fs@data[cm_ind,4] <- apply(spdf_fs@data[cm_ind,1:2], 1, trans.y)
    }
}
coordinates(srtm)  = ~u+v
proj4string(srtm)  <- CRS(as.character(crs1)) 
cm_raster2 <- rasterize(srtm, raster(sgdf), field='z', update=TRUE, proj4string = CRS(as.character(crs1)))
cm_raster2[which(getValues(cm_raster2) == 0)] <- NA
cm_raster2 <- focal(cm_raster2, w= matrix(rep(1,361), nrow=19),  na.rm=TRUE, mean, NAonly=TRUE) 
cm_raster2 <- focal(cm_raster2, w= matrix(rep(1,361), nrow=19),  na.rm=TRUE, mean, NAonly=TRUE) 
cm_raster2 <- focal(cm_raster2, w= matrix(rep(1,361), nrow=19),  na.rm=TRUE, mean, NAonly=TRUE) 

sgdf_srtm_v <- cm_raster2
srtm_v_slope  <- terrain(sgdf_srtm_v, opt='slope')
srtm_v_aspect <- terrain(sgdf_srtm_v, opt='aspect')
srtm_v_shade  <- hillShade(srtm_v_slope, srtm_cd_aspect, 40, 270)
spdf_fs_v <- data.frame(cbind(u=spdf_fs@data[,3], v=spdf_fs@data[,4]))
coordinates(spdf_fs_v)  = ~u+v

pdf("5pictures/mod_v.pdf", height=4.3, width=6) 
par(mar = c(0, 0, 0, 0))
plot(sgdf_srtm_v, col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE)
plot(srtm_v_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
points(spdf_fs_v, cex=1, col="#FF004E", pch=16)
dev.off() 

png("5pictures/mod_v.png", height=1787*4, width=2135*4, units = "px") 
par(mar = c(0, 0, 0, 0))
plot(sgdf_srtm_v, col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE)
plot(srtm_v_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
points(spdf_fs_v, cex=20, col="#FF004E", pch=16)
dev.off() 



# 9. export ====================================





pdf("5pictures/mod_xx.pdf", height=4.3, width=6) 
par(mar = c(0, 0, 0, 0))
    plot(raster(sgdf_srtm), col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE)
    plot(srtm_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
    #points(spdf_h, cex=0.3, col=grey, pch=16)
    points(spdf_fs, cex=1, col="#FF004E", pch=16)
dev.off() 


png("5pictures/mod_xx.png", height=1787*4, width=2135*4, units = "px") 
par(mar = c(0, 0, 0, 0))
plot(raster(sgdf_srtm), col = top.colors(20), axes=FALSE, box=FALSE, legend=FALSE)
plot(srtm_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
#points(spdf_h, cex=0.3, col=grey, pch=16)
points(spdf_fs, cex=20, col="#FF004E", pch=16)
dev.off() 


 



save.image("7ws/ws01.rws")




#for FILE in ./*.pdf; do
#pdfcrop "${FILE}" "${FILE}"
#done



