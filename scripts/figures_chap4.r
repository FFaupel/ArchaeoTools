
arbverz <- "/home/fon/daten/analyse/modproj_qaam"  # variable with the path
setwd(arbverz)                                 # set the working directory

file_meg    <- "1data/meg_dw.csv"          # megalithic tombs
file_tum    <- "1data/tum_dw.csv"          # tumuli
file_coast1 <- "2geodata/coast_gk3.shp"    # Shapefile with folder and extension
file_coast  <- "coast_gk3"                 # Shapefile without folder and extension
file_srtm   <- "2geodata/dw_gk3_50_ag.asc" # srtm DTM
file_vil    <- "1data/villages.xls"        # foundation of villages in dw

crs1 <- "+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=WGS84 +units=m +no_defs" # projected gk3 (command in two lines)
crs2 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"  # unprojected WGS84

spdf_meg <- read.table(file_meg, sep=';', header=TRUE)
spdf_tum <- read.table(file_tum, sep=';', header=TRUE)
df_meg <- read.table(file_meg, sep=';', header=TRUE)
df_tum <- read.table(file_tum, sep=';', header=TRUE)

library("sp")
library("proj4")

coordinates(spdf_meg)=~x+y  
coordinates(spdf_tum)=~x+y  

avs <- paste(arbverz,"/2geodata",sep="")
#coast <- readOGR(arbverz, p4s=NULL, file_coast1) 

library("gdata")
df_vil_wgs84 <- read.xls(file_vil, 1)
spdf_vil_wgs84 <- df_vil_wgs84
coordinates(spdf_vil_wgs84)=~x+y  

sgdf_srtm <- read.asciigrid(file_srtm) 

df_vil_coord <- project(cbind(spdf_vil_wgs84@coords[,1], spdf_vil_wgs84@coords[,2]), 
                        crs1)
df_vil_k <- cbind(x=df_vil_coord[,1],y=df_vil_coord[,2])
df_vil <- data.frame(id=df_vil_wgs84[,1], village=as.character(df_vil_wgs84[,2]), 
                     AD=df_vil_wgs84[,3])
spdf_vil <- SpatialPointsDataFrame(df_vil_k, as.data.frame(df_vil), 
                                   proj4string=CRS(as.character(crs1)))

proj4string(spdf_meg)  <- CRS(as.character(crs1)) 
proj4string(spdf_tum)  <- CRS(as.character(crs1)) 
proj4string(sgdf_srtm) <- CRS(as.character(crs1)) 




cb <- c(1200,1250,1300,1350,1400)
count <- 1:4
for (i in 1:4) {  
    higher <- which(df_vil[,3] > cb[i])
    lower  <- which(df_vil[,3] < cb[i+1])
    hl <- intersect(higher,lower)
    length(hl)
    count[i] <- length(hl)
}
years <- c("1200-1250","1250-1300","1300-1350","1350-1400")
v.df <- data.frame(years,count)
ggplot(v.df, aes(years,count)) +
    geom_bar(stat="identity") +
    ggtitle("Village foundations in different periods")

vil_fd <- df_vil[,3]
hist(df_vil[,3],col="black",6)

library(ggplot2)

pdf("/home/fon/daten/pub/springer/ssmod_v06/figures/ss2013_v04-EM-Vil-1_GGPLOT",paper = "special",width = 8, height = 4)
ggplot(df_vil, aes(AD)) +
    geom_histogram(binwidth = 20) +
    xlab("Time A.D.") +
    ylab("Frequency") +
    ggtitle("Histogram of village foundations in different periods")
dev.off()

library("KernSmooth")
ks_vil <- bkde(df_vil[,3], kernel="normal", 
   bandwidth=5, gridsize=201, range.x = 
   c(1200,1400))
plot(ks_vil,col="black", pch=16)
lines(ks_vil)

pdf("/home/fon/daten/pub/springer/ssmod_v06/figures/ss2013_v04-EM-Vil-2_GGPLOT",paper = "special",width = 8, height = 4)
ks_vil <- data.frame(ks_vil)
ggplot(ks_vil, aes(ks_vil$x,ks_vil$y)) +
    xlim(1240,1360) +
    xlab("Time A.D.") + 
    ylab("density of village foundations") +
    theme(axis.title.y = element_text(vjust=.3)) + 
    geom_point() + 
    geom_line()
dev.off()

plot(df_vil[,3],df_vil[,1],col="black", pch=16)
lines(df_vil[,3],df_vil[,1])

pdf("/home/fon/daten/pub/springer/ssmod_v06/figures/ss2013_v04-EM-Vil-3_GGPLOT",paper = "special",width = 8, height = 4)
ggplot(df_vil, aes(AD,id)) +
    xlab("Time A.D.") + 
    geom_point() + 
    geom_line()
dev.off()

interval <- c(df_vil[,3],df_vil[13,3]) - c(df_vil[1,3],df_vil[,3])
plot(interval[2:12],col="black", pch=16)
lines(interval[2:12])

pdf("/home/fon/daten/pub/springer/ssmod_v06/figures/ss2013_v04-EM-Vil-4_GGPLOT",paper = "special",width = 8, height = 4)
qplot(seq_along(interval[2:12]), interval[2:12]) +
    xlab("Index") +
    ylab("interval") +
    geom_point() + 
    geom_line()
dev.off()

file_meg    <- "1data/meg_dw.csv" 
df_meg <- read.table(file_meg, sep=';', header=TRUE)
spdf_meg <- read.table(file_meg, sep=';', header=TRUE)
library(sp)
coordinates(spdf_meg)=~x+y
crs1 <- "+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=WGS84 +units=m +no_defs"
proj4string(spdf_meg)  <- CRS(as.character(crs1)) 
file_srtm   <- "2geodata/dw_gk3_50_ag.asc"
sgdf_srtm <- read.asciigrid(file_srtm) 
bb = bbox(sgdf_srtm)    
library(spatstat)
win <- owin(xrange=c(bb[1,1],bb[1,2]), yrange=c(bb[2,1],bb[2,2]), unitname="m")
spdf_meg <-remove.duplicates(spdf_meg, zero=0,remove.second=TRUE)
ppp_meg <- ppp(spdf_meg@coords[,1], spdf_meg@coords[,2], window=win)

anzahl <- ppp_meg$n
dx <- (ppp_meg$window$xrange[2] - 
       ppp_meg$window$xrange[1]) / 1000
dy <- (ppp_meg$window$yrange[2] - 
       ppp_meg$window$yrange[1]) / 1000
dichte_1 <- anzahl / (dx*dy)
dichte_1

ch <- convexhull.xy(ppp_meg$x, ppp_meg$y) 
fl <- area.owin(ch)      
dichte_2 <- anzahl / fl  
dichte_3 <- ppp_meg$n / area.owin(convexhull.xy(ppp_meg$x, ppp_meg$y))
f <- (dichte_2 * 1000000) / dichte_1 
f

rw      <- 3000    
xmin    <- ppp_meg$window$xrange[1] - rw/2
xmax    <- ppp_meg$window$xrange[2] + rw/2
ymin    <- ppp_meg$window$yrange[1] - rw/2
ymax    <- ppp_meg$window$yrange[2] + rw/2
zeilen  <- round((ymax-ymin)/rw, 0) + 1
spalten <- round((xmax-xmin)/rw, 0) + 1
z <- cbind(1:(spalten*zeilen))
df <- data.frame(z)
gt <- GridTopology(c(ppp_meg$window$xrange[1] - rw/2,ppp_meg$window$yrange[1] - rw/2), c(rw,rw), c(spalten,zeilen))
sgdf <- SpatialGridDataFrame(gt, df, proj4string = CRS(as.character(crs1)))
gt
for (i in seq(along=coordinates(gt)[,1])){
    x <- coordinates(gt)[i,1] - rw/2
    y <- coordinates(gt)[i,2] - rw/2
    xi <- which(ppp_meg$x>x & ppp_meg$x<x+rw)  
    yi <- which(ppp_meg$y>y & ppp_meg$y<y+rw)  
    pz <- length(intersect(xi,yi))   
    sgdf@data$z[i]<- pz / (rw/1000)^2 
}
image(sgdf, col = gray.colors(25, start = 0.97, end = 0.4))      
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)

pdf("/home/fon/daten/pub/springer/ssmod_v06/figures/ss2013_v04-density_chess_board_EM-Meg-3-GGPLOT.pdf",paper = "special",width = 7, height = 4)
df_sgdf <- data.frame(sgdf)
colnames(df_sgdf) <- c("density","x","y")
ggplot(data=df_sgdf) +
    coord_equal() + labs(x=NULL, y=NULL) + 
    geom_tile(aes(x, y,fill=density)) +
    scale_fill_gradient(low="#ffffff", high="#999999") +
    geom_point(data = df_meg, aes(x = x,y = y))
dev.off()

sgdf_kde <- sgdf
sd      <- 3000
for (i in seq(along=coordinates(gt)[,1])){
    x <- coordinates(gt)[i,1]
     y <- coordinates(gt)[i,2]
     g2 <- 0
     for (j in seq(along=ppp_meg$x)){  
         distanz <- sqrt((ppp_meg$x[j] - x)^2 + 
                         (ppp_meg$y[j] - y)^2)
         g1 <- dnorm(distanz, mean=0, sd=sd)
         g2 <-g2 + g1
     }
     sgdf_kde@data$z[i]<- g2     
 }
image(sgdf_kde, col = gray.colors(25, start = 
                    0.97, end = 0.4))   
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4) 

pdf("/home/fon/daten/pub/springer/ssmod_v06/figures/ss2013_v04-kde_EM-Meg-4-GGPLOT.pdf",paper = "special",width = 7, height = 4)
df_sgdf_kde <- data.frame(sgdf_kde)
colnames(df_sgdf_kde) <- c("density","x","y")
ggplot(data=df_sgdf_kde) +
    coord_equal() + labs(x=NULL, y=NULL) + 
    geom_tile(aes(x, y,fill=density)) +
    scale_fill_gradient(low="#ffffff", high="#999999") +
    geom_point(data = df_meg, aes(x = x,y = y))
dev.off()


rw <- 1000  
sd <- 2000
dens_p <- density(ppp_meg, sd, edge=TRUE, 
   at="points")   
dens_r5 <- density(ppp_meg, sd, eps=rw, 
   edge=TRUE, at="pixels") 
plot(dens_r5, col = gray.colors(25, start = 0.97, 
   end = 0.4))
contour(dens_r5, add=T)    
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)  
 

pdf("/home/fon/daten/pub/springer/ssmod_v06/figures/ss2013_v04-kde_EM-Meg-5-GGPLOT.pdf",paper = "special",width = 7, height = 4)
df_dens_r5 <- data.frame(dens_r5)
colnames(df_dens_r5) <- c("x","y","density")
ggplot(data=df_dens_r5) +
    coord_equal() + labs(x=NULL, y=NULL) + 
    geom_tile(aes(x, y,fill=density)) +
    stat_contour(aes(x, y,z=density),colour="#000000") +
    geom_point(data = df_meg, aes(x = x,y = y)) +
    scale_fill_gradient(low="#ffffff", high="#999999")
dev.off()


sdev <- 3*mean(nndist(ppp_meg))  
dens_r6 <- density(ppp_meg, sdev, eps=rw, 
                   edge=TRUE, at="pixels")
plot(dens_r6, col = gray.colors(25, start = 
                  0.97, end = 0.4))
contour(dens_r6, add=T)    
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)  

pdf("/home/fon/daten/pub/springer/ssmod_v06/figures/ss2013_v04-kde_EM-Meg-6-GGPLOT.pdf",paper = "special",width = 7, height = 4)
df_dens_r6 <- data.frame(dens_r6)
colnames(df_dens_r6) <- c("x","y","density")
ggplot(data=df_dens_r6) +
    coord_equal() + labs(x=NULL, y=NULL) + 
    geom_tile(aes(x, y,fill=density)) +
    stat_contour(aes(x, y,z=density),colour="#000000") +
    geom_point(data = df_meg, aes(x = x,y = y)) +
    scale_fill_gradient(low="#ffffff", high="#999999")
dev.off()


dens_r <- density(ppp_meg, bw = "nrd", eps=rw, 
   edge=TRUE, at="pixels")
plot(dens_r, col = gray.colors(25, start = 0.97, 
   end = 0.4))
contour(dens_r, add=T) 
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)  

pdf("/home/fon/daten/pub/springer/ssmod_v06/figures/ss2013_v04-kde_EM-Meg-7-GGPLOT.pdf",paper = "special",width = 7, height = 4)
df_dens_r <- data.frame(dens_r)
colnames(df_dens_r) <- c("x","y","density")
ggplot(data=df_dens_r) +
    coord_equal() + labs(x=NULL, y=NULL) + 
    geom_tile(aes(x, y,fill=density)) +
    stat_contour(aes(x, y,z=density),colour="#000000") +
    geom_point(data = df_meg, aes(x = x,y = y)) +
    scale_fill_gradient(low="#ffffff", high="#999999")
dev.off()



#gt   <- getGridTopology(sgdf_srtm)


#sgdf <- SpatialGridDataFrame(gt, df, proj4string = CRS(as.character(crs1)))
#sgdf <- sgdf_srtm



rw <- 1000
fs <- cbind(x=spdf_meg@coords[,1],y=spdf_meg@coords[,2])
zeilen  <- round((bbox(spdf_meg)[2,2]-bbox(spdf_meg)[2,1])/rw, 0) + 2
spalten <- round((bbox(spdf_meg)[1,2]-bbox(spdf_meg)[1,1])/rw, 0) + 2
z    <- cbind(1:(spalten*zeilen))
df   <- data.frame(cbind(1:((round((bbox(spdf_meg)[2,2]-bbox(spdf_meg)[2,1])/rw, 0) + 2)*(round((bbox(spdf_meg)[1,2]-bbox(spdf_meg)[1,1])/rw, 0) + 2))))
gt   <- GridTopology(c(bbox(spdf_meg)[1,1] - rw/2,bbox(spdf_meg)[2,1] - rw/2),c(rw,rw), c(spalten,zeilen))
sgdf <- SpatialGridDataFrame(gt, df, proj4string = CRS(as.character(crs1)))

library(tripack)
fsv   <- voronoi.mosaic(spdf_meg$x, spdf_meg$y, duplicate = 'remove') 
rad <- fsv$radius
fsvsp <- SpatialPointsDataFrame(cbind(fsv$x, fsv$y), as.data.frame(rad), proj4string=CRS(as.character(crs1))) 

library(gstat)
g <- gstat(formula=fsvsp@data$rad ~ 1, data=fsvsp, nmin = 5, maxdist = 10000, nmax = 15) #8
vt <- variogram(g)
v.fit <- fit.variogram(vt, vgm(1, "Gau", 10000, 1), fit.sills = TRUE, fit.ranges = TRUE,fit.method = 1)
g <- gstat(g, id="var1", model=v.fit )
k <- predict(g, model=v.fit, newdata=sgdf)



image(k, col = gray.colors(25, start = 0.4, end = 0.97))
contour(k, add=T)    
points(spdf_meg, pch=16, cex=0.4)



pdf("/home/fon/daten/pub/springer/ssmod_v06/figures/em-meg-8-GGPLOT.pdf",paper = "special",width = 7, height = 4)
df_k <- data.frame(k)
colnames(df_k) <- c("radius","var","x","y")
ggplot(data=df_k) +
    coord_equal() + labs(x=NULL, y=NULL) + 
    geom_tile(aes(x, y,fill=radius)) +
    stat_contour(aes(x, y,z=radius),colour="#000000") +
    geom_point(data = df_meg, aes(x = x,y = y)) +
    scale_fill_gradient(low="#999999", high="#ffffff")
dev.off()


 library('raster')

# KRIGIN CODE v2
pdf("/home/fon/daten/pub/springer/ssmod_v06/figures/em-meg-8-v2.pdf",paper = "special",width = 7, height = 5.33)
    plot(raster(k), col = gray.colors(25, start = 0.4, end = 0.97), cex.axis = .9)
    contour(k, add=T)    
    points(ppp_meg$x, ppp_meg$y, pch=20)
dev.off()
# KRIGIN CODE v3
pdf("/home/fon/daten/pub/springer/ssmod_v06/figures/em-meg-8-v2.pdf",paper = "special",width = 7, height = 5.33)
    plot(raster(k), col = gray.colors(25, start = 0.4, end = 0.97),  cex.axis = .9, legend.shrink=0.25, legend.width=1,    legend.args=list(text="          Radius (m)",line=.5))
    contour(k, add=T)    
    points(ppp_meg$x, ppp_meg$y, pch=20)
dev.off()


