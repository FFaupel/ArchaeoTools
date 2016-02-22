# Please improve the script. 
#   add comments
#   change variable names
#   apply style guid
#   enhance code
#   finde artificial and natural errors
#   add additions from the exercises


################################################################################
## Didactic R-Script for GSHDL Modelling Summer School
## =============================================================================
## Project: GSHDL Summer School
## Author: Oliver Nakoinz
## Version: 01
## Date of last changes: 30.07.2013
## Data: srtm, monuments
## Author of data: gshdl
## Purpose: Didactic
## Inhalt: 1. preparation, 2. data import, ...
## Licence data: -
## Licence Script: GPL
################################################################################

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

library("sp")
library("proj4")

coordinates(spdf_meg)=~x+y  
coordinates(spdf_tum)=~x+y  

avs <- paste(arbverz,"/2geodata",sep="")
coast <- readOGR(arbverz, p4s=NULL, file_coast1) 

library("gdata")
df_vil_wgs84 <- read.xls(file_vil, 1)
spdf_vil_wgs84 <- df_vil_wgs84
coordinates(spdf_vil_wgs84)=~x+y  

sgdf_srtm <- read.asciigrid(file_srtm) 

spdf_vil_wgs84@coords[,1]        
spdf_vil_wgs84@coords[,"x"]  
spdf_vil_wgs84@coords[,2]  

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

save.image("4ws/ws01.rws")

hist(df_vil[,3],col="black",6)

library("KernSmooth")
ks_vil <- bkde(df_vil[,3], kernel="normal", bandwidth=5, gridsize=201, range.x = c(1200,1400))
plot(ks_vil,col="black", pch=16)
lines(ks_vil)

################ time series

ts_vil <- ts(ks_vil$y, start=c(1200), end=c(1400), frequency = 1) 
plot(ts_vil)



pdf("/home/fon/daten/pub/springer/ssmod_v06/figures/4dens_ts1.pdf", 4.65, 3, title="Time Series",    pointsize=11)
plot(ts_vil,ylab="density")
dev.off()


pdf("/home/fon/daten/pub/springer/ssmod_v06/figures/4dens_ts2.pdf", 4.65, 3, title="Time Series",     pointsize=11)
acf(ts_vil, lag.max = 100)
dev.off()


acf(ts_vil)
acf(ts_vil, lag.max = 100)


ks_vil[,2]

spectrum(ts_vil)


library("TTR")
ts_vil_s10 <- SMA(ts_vil,n=10)
ts_vil_s25 <- SMA(ts_vil,n=25)
ts_vil_s50 <- SMA(ts_vil,n=50)
plot(ts_vil_s10)
lines(ts_vil_s25, lty=2)
lines(ts_vil_s50, lty=3)
plot(ts_vil-ts_vil_s50)


> ts_vil_dec <- decompose(ts_vil)
Fehler in decompose(ts_vil) : 
  Zeitreihe hat keine oder weniger als 2 Perioden

#######################################################xxxxxxxxxxxxxxxxxxxxx
























#######################################################xxxxxxxxxxxxxxxxxxxxx


plot(df_vil[,3],df_vil[,1],col="black", pch=16)
lines(df_vil[,3],df_vil[,1])  

plot(interval[2:12],col="black", pch=16)
lines(interval[2:12]) 

plot(df_vil[,3],df_vil[,1],col="black", pch=16)
lines(c(df_vil[1,3],df_vil[13,3]),c(df_vil[1,1],df_vil[13,1]) ) 
years <- 1259:1350
y <- 1 + 0.003* (years-1259)^2
lines (years,y)

E <-  df_vil[,1]
T1 <- 1 + (12/91)*(df_vil[,3] - 1259)
T2 <- 1 + 0.003* (df_vil[,3]-1259)^2
sum(((E-T1)^2)/T1)
sum(((E-T2)^2)/T2)

yr <- df_vil[,3] - 1259  # first we deal with c
lm_vil <- lm(E ~ yr)
coef(lm_vil)

plot(df_vil[,3]- 1259,df_vil[,1],col="black", pch=16)
abline(lm_vil)

flm_vil <- fitted(lm_vil)
sum(((E-flm_vil)^2)/flm_vil)

lm_vil2 <- lm(E ~ I(yr^2))
coef(lm_vil2)
flm_vil2 <- fitted(lm_vil2)
sum(((E-flm_vil2)^2)/flm_vil2)

lm_vil3 <- lm(E ~ poly(yr, 4, raw=TRUE))
coef(lm_vil3)
flm_vil3 <- fitted(lm_vil3)
sum(((E-flm_vil3)^2)/flm_vil3)

plot(df_vil[,3]- 1259,df_vil[,1],col="black", pch=16)
lines(yr,flm_vil3)

chi <- 1:20
for(i in seq(1::20)){
    lm_vil4 <- lm(E ~ poly(yr, i, raw=TRUE))
    flm_vil4 <- fitted(lm_vil4)
    chi[i] <- sum(((E-flm_vil4)^2)/flm_vil4)
}
chi
plot(chi,col="black", pch=16)
lines(chi)

anzahl <- ppp_meg$n
dx <- (ppp_meg$window$xrange[2] - ppp_meg$window$xrange[1]) / 1000
dy <- (ppp_meg$window$yrange[2] - ppp_meg$window$yrange[1]) / 1000
dichte_1 <- anzahl / (dx*dy)
dichte_1

library(spatstat)  
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
gt <- GridTopology(c(ppp_meg$window$xrange[1] - rw/2,ppp_meg$window$yrange[1] - 
                         rw/2), c(rw,rw), c(spalten,zeilen))
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

rw      <- 3000    
xmin    <- ppp_meg$window$xrange[1] - rw/2
xmax    <- ppp_meg$window$xrange[2] + rw/2
ymin    <- ppp_meg$window$yrange[1] - rw/2
ymax    <- ppp_meg$window$yrange[2] + rw/2
zeilen  <- round((ymax-ymin)/rw, 0) + 1
spalten <- round((xmax-xmin)/rw, 0) + 1
z <- cbind(1:(spalten*zeilen))
df <- data.frame(z)
gt <- GridTopology(c(ppp_meg$window$xrange[1] - rw/2,ppp_meg$window$yrange[1] - 
                         rw/2), c(rw,rw), c(spalten,zeilen))
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
image(sgdff, col = gray.colors(25, start = 0.97, end = 0.4))      
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)    

sgdf_kde <- sgdf
sd      <- 3000
for (i in seq(along=coordinates(gt)[,1])){
    x <- coordinates(gt)[i,1]
    y <- coordinates(gt)[i,2]
    g2 <- 0
    for (j in seq(along=ppp_meg$x)){  
        distanz <- sqrt((ppp_meg$x[j] - x)^2 + (ppp_meg$y[j] - y)^2)
        g1 <- dnorm(distanz, mean=0, sd=sd)
        g2 <-g2 + g1
    }
    sgdf_kde@data$z[i]<- g2     
}
image(sgdf_kde, col = gray.colors(25, start = 0.97, end = 0.4))   
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4) 

sgdf_kde <- sgdf
sd      <- 3000
for (i in seq(along=coordinates(gt)[,1])){
    x <- coordinates(gt)[i,1]
    y <- coordinates(gt)[i,2]
    g2 <- 0
    for (j in seq(along=ppp_meg$x)){  
        distanz <- sqrt((ppp_meg$x[j] - x)^2 + (ppp_meg$y[j] - y)^2)
        g1 <- dnorm(distanz, mean=0, sd=sd)
        g2 <-g2 + g1
    }
    sgdf_kde@data$z[i]<- g2     
}
image(sgdf_kde, col = gray.colors(25, start = 0.97, end = 0.4))   
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)   

library(spatstat)
rw <- 1000  
sd <- 2000
dens_p <- density(ppp_meg, sd, edge=TRUE, at="points")   
dens_r5 <- density(ppp_meg, sd, eps=rw, edge=TRUE, at="pixels") 
plot(dens_r5, col = gray.colors(25, start = 0.97, end = 0.4))
contour(dens_r5, add=T)    
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)  

library(spatstat)
rw <- 1000  
sd <- 2000
dens_p <- density(ppp_meg, sd, edge=TRUE, at="points")   
dens_r5 <- density(ppp_meg, sd, eps=rw, edge=TRUE, at="pixels") 
plot(dens_r5, col = gray.colors(25, start = 0.97, end = 0.4))
contour(dens_r5, add=T)    
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)  

sdev <- 3*mean(nndist(ppp_meg))  
dens_r6 <- density(ppp_meg, sdev, eps=rw, edge=TRUE, at="pixels")
plot(dens_r6, col = gray.colors(25, start = 0.97, end = 0.4))
contour(dens_r6, add=T)    
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)  

library(spatstat)
sdev <- 3*mean(nndist(ppp_meg))  
dens_r6 <- density(ppp_meg, sdev, eps=rw, edge=TRUE, at="pixels")
plot(dens_r6, col = gray.colors(25, start = 0.97, end = 0.4))
contour(dens_r6, add=T)    
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)  

dens_r <- density(ppp_meg, bw = "nrd", eps=rw, edge=TRUE, at="pixels")
plot(dens_r, col = gray.colors(25, start = 0.97, end = 0.4))
contour(dens_r, add=T) 
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4) 

dens_r <- density(ppp_meg, bw = "nrd", eps=rw, edge=TRUE, at="pixels")
plot(dens_r, col = gray.colors(25, start = 0.97, end = 0.4))
contour(dens_r, add=T)    
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)  

library(tripack); library(gstat)
rw <- 1000
fs <- cbind(x=spdf_meg@coords[,1],y=spdf_meg@coords[,2])
zeilen  <- round((bbox(spdf_meg)[2,2]-bbox(spdf_meg)[2,1])/rw, 0) + 2
spalten <- round((bbox(spdf_meg)[1,2]-bbox(spdf_meg)[1,1])/rw, 0) + 2
z    <- cbind(1:(spalten*zeilen))
df   <- data.frame(cbind(1:((round((bbox(spdf_meg)[2,2]-bbox(spdf_meg)[2,1])
                                   /rw, 0) + 2)*(round((bbox(spdf_meg)[1,2]-bbox(spdf_meg)[1,1])/rw, 0) + 2))))
gt   <- GridTopology(c(ppp_meg$window$xrange[1] - rw/2,ppp_meg$window$yrange[1] - 
                           rw/2), c(rw,rw), c(spalten,zeilen))
sgdf <- SpatialGridDataFrame(gt, df, proj4string = CRS(as.character(crs1)))

fsv   <- voronoi.mosaic(ppp_meg$x, ppp_meg$y, duplicate = 'remove') 
fsvsp <- SpatialPointsDataFrame(cbind(fsv$x, fsv$y), as.data.frame(fsv$x), 
                                proj4string=CRS(as.character(crs1)))
fspv  <- ppp(fsvsp@coords[,1], fsvsp@coords[,2], window=win) 
fs_vd <- cbind(fspv$x,fspv$y,nncross(fspv,ppp_meg)$dist)            
fs_vd_spdf <- SpatialPointsDataFrame(cbind(fs_vd[,1],fs_vd[,2]), 
                                     as.data.frame(fs_vd[,3]), proj4string=CRS(as.character(crs1)))

vt    <- variogram(fs_vd_spdf@data$fs_vd ~ 1, fs_vd_spdf)
v.fit <- fit.variogram(vt, vgm(1, "Gau", 10000, 1), fit.sills = TRUE, 
                       fit.ranges = TRUE,fit.method = 1)
k     <- krige(fs_vd_spdf@data$fs_vd ~ 1, fs_vd_spdf, sgdf, v.fit, 
               nmin = 3, maxdist = 10000, nmax = 8)
image(k, col = gray.colors(25, start = 0.4, end = 0.97))
contour(k, add=T)    
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)  

diff6_5 <- dens_r5
diff6_5$v <- dens_r6$v - dens_r5$v

save.image("4ws/ws02.rws")

library(gstat)
meg_idw     <- idw(fs_vd_spdf@data$fs_vd ~ 1, fs_vd_spdf, sgdf)
image(meg_idw, col = gray.colors(25, start = 0.4, end = 0.97))
contour(k, add=T)    
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)  

library(maptools)
rw <- 500  
sd <- 2000
dens_r <- density(ppp_meg, sd, eps=rw, edge=TRUE, at="pixels") 
#plot(dens_r, col = gray.colors(25, start = 0.97, end = 0.4))
df_dens_r <- as.SpatialGridDataFrame.im(dens_r)
meg_kde_samppoints <- spsample(df_dens_r, 1000, type="random")
meg_kde_samp <- overlay(x=df_dens_r, y=meg_kde_samppoints)

vt2    <- variogram(meg_kde_samp@data$v ~ 1, meg_kde_samp)
v.fit2 <- fit.variogram(vt2, vgm(1, "Gau", 5000, 1), fit.sills = TRUE, 
                        fit.ranges = TRUE,fit.method = 7)
plot(vt2,v.fit2)

k2     <- krige(meg_kde_samp@data$v ~ 1, meg_kde_samp, df_dens_r, v.fit2, 
                nmin = 3, maxdist = 10000, nmax = 8)
image(k2, col = gray.colors(25, start = 0.97, end = 0.4))
contour(k2, add=T)    
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)  

library(spatstat)
sdev <- 2*mean(nndist(ppp_meg)+mean(nndist(ppp_tum)))
bb   <- bbox(sgdf_srtm)    
win  <- owin(xrange=c(bb[1,1],bb[1,2]), yrange=c(bb[2,1],bb[2,2]), unitname="m")
meg_dens <- density(ppp_meg, kernel="gaussian", bw=sdev, dimyx=c(36,56), w=win, 
                    edge=TRUE, at="pixels")
tum_dens <- density(ppp_tum, kernel="gaussian", bw=sdev, dimyx=c(36,56), w=win, 
                    edge=TRUE, at="pixels")

library(maptools)
sgdf_meg_dens     <- as.SpatialGridDataFrame.im(meg_dens)
meg_dens_samppt   <- spsample(sgdf_meg_dens, 500, type="random")
meg_dens_samppt2  <- spsample(sgdf_meg_dens, 500, type="random")
meg_dens_samp     <- overlay(x=sgdf_meg_dens, y=meg_dens_samppt)
meg_dens_samp2    <- overlay(x=sgdf_meg_dens, y=meg_dens_samppt2)

sgdf_tum_dens     <- as.SpatialGridDataFrame.im(tum_dens)
tum_dens_samppt   <- spsample(sgdf_tum_dens, 500, type="random")
tum_dens_samp     <- overlay(x=sgdf_tum_dens, y=tum_dens_samppt)

dens_samp         <- meg_dens_samp
names(dens_samp)[names(dens_samp) == 'v'] <- 'meg'
dens_samp@data    <- cbind(dens_samp@data,tum_dens_samp@data$v)
names(dens_samp)[names(dens_samp) == 'tum_dens_samp@data$v'] <- 'tum'
dens_samp <- dens_samp[-which(is.na(dens_samp@data$meg)),]

library(maptools)
cor.test(dens_samp@data$meg, dens_samp@data$tum, method="p")
ks.test(dens_samp@data$meg, dens_samp@data$tum)
ks.test(dens_samp@data$meg, meg_dens_samp2@data$v)

meg_env_g <- envelope(ppp_meg, fun = Gest, nrank = 2, nsim = 99)   
plot(meg_env_g)

meg_env_f <- envelope(ppp_meg, fun = Fest, nrank = 2, nsim = 99)   
plot(meg_env_f)

meg_env_k <- envelope(ppp_meg, fun = Kest, nrank = 2, nsim = 99)   
plot(meg_env_k)

library(spatstat)
modspec <- list(cif="straush",par=list(beta=2,gamma=0.2,r=0.7,hc=0.3), 
                w=c(bb[1,1],bb[1,2],bb[2,1],bb[2,2]))
modsim <- rmh(model=modspec,start=list(n.start=200), control=list(nrep=10,nverb=5))
plot(modsim)

meg_env_t <- envelope(ppp_meg, fun = Tstat, nrank = 2, nsim = 20)   
plot(meg_env_t)

library(spatstat)
ch_meg <- convexhull(ppp_meg)
ch_tum <- convexhull(ppp_tum)
plot(ch_tum, border="grey", main="")
plot(ch_meg, add=TRUE)
points(ppp_tum$x, ppp_tum$y, pch=17, cex=0.6, col="grey")  
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)

###it takes some time, so it is commented
#buf_tum <- dilation(ppp_tum, 1000, polygonal=TRUE, tight=F) 
#buf_meg <- dilation(ppp_meg, 1000, polygonal=TRUE, tight=F) 
#plot(buf_tum, border="grey", main="")
#plot(buf_meg, add=TRUE)
#points(ppp_tum$x, ppp_tum$y, pch=17, cex=0.6, col="grey")  
#points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)  

library("classInt")
nb_meg <- classIntervals(dens_samp@data$meg, style = "fisher", dataPrecision = NULL)  
nb_tum <- classIntervals(dens_samp@data$tum, style = "fisher", dataPrecision = NULL)  
contour(sgdf_meg_dens, add=F, method = "edge", levels = nb_meg$brks, drawlabels = F)    
contour(sgdf_tum_dens, add=T, method = "edge", levels = nb_tum$brks, drawlabels = FF, col="grey")  
points(ppp_tum$x, ppp_tum$y, pch=17, cex=0.6, col="grey")  
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)

library("cluster")
dens_samp_clus <- pam(dens_samp@data, 5)
dens_samp@data <- cbind(dens_samp@data,dens_samp_clus$clustering)
names(dens_samp)[names(dens_samp) == 'dens_samp_clus$clustering'] <- 'clus'
plot(dens_samp, pch=dens_samp@data$clus)

dpd <- ppp_meg
dpd$d <- density(ppp_meg, 1000, edge=TRUE, at="points")          
k <- 1                                                         
s <- -20000000                                                 
dpd$id <- seq(along=dpd$x); dpd$v_id <- 0; dpd$v_x <- 0; dpd$v_y <- 0; 
dpd$v_d <- 0; dpd$dist <- 0; dpd$e <- 0 
dpd.d <- dist(cbind(dpd$x,dpd$y),upper = T)               
dpd.m <- as.matrix(dpd.d)                                   
maxd <- max(dpd$d)
maxm <- max(dpd.m)
dpd$d <- dpd$d * maxm / maxd                           
dpd.e <- dpd.m                               

for (i in seq(along=dpd$d))   {                     
    dpd.e[i,] <- (dpd$d - dpd$d[i]) - (k * dpd.m[i,])         
    w <-max(dpd.e[i,])                                        
    wi <- which(dpd.e[i,] == w))                             
    if (i == wi | w < s)    {dpd$v_id[i] <- 0}              
    else                    {dpd$v_id[i] <- wi}
}
dpd$v_id      

for (i in seq(along=dpd$d))   {                     
    dpd.e[i,] <- (dpd$d - dpd$d[i]) - (k * dpd.m[i,])         
    w <-max(dpd.e[i,])                                        
    wi <- which(dpd.e[i,] == w)                                
    if (i == wi | w < s)    {dpd$v_id[i] <- 0}              
    else                    {dpd$v_id[i] <- wi}
}
dpd$v_id                                               
for (i in seq(along=dpd$d))   {                           
    if (dpd$v_id[i] > 0)  {
        dpd$v_x[i] <- dpd$x[dpd$v_id[i]]
        dpd$v_y[i] <- dpd$y[dpd$v_id[i]]
        dpd$v_d[i] <- dpd$d[dpd$v_id[i]]
        dpd$dist[i]<- dpd.m[i,dpd$v_id[i]]
        dpd$e[i]   <- dpd.e[i,dpd$v_id[i]]
    }
    else {
        dpd$v_x[i] <- dpd$x[i]
        dpd$v_y[i] <- dpd$y[i]
        dpd$v_d[i] <- dpd$d[i]
    }
}
dpd

LinesList <- list()                                          
for(i in seq(along=dpd$id)) {                   
    m <- matrix(data = c(dpd$x[i],dpd$v_x[i],dpd$y[i],dpd$v_y[i]), nrow=2, ncol=2)   
    L <- Line(m)
    LL <- list(L)
    name  <- paste("Zuordnung_", dpd$id,"_", dpd$v_id, sep="") 
    LLL <- Lines(LL, ID = name[i])                        
    LinesList[length(LinesList)+1] <- LLL  
}
sl <- SpatialLines(LinesList, proj4string = CRS(crs1))

sdf <- SpatialLinesDataFrame(sl, as.data.frame(dpd$id), match.ID = FALSE)
plot(sdf)
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4) 

ras <- sgdf_tum_dens
r <- 5000      
ras@data$v[which(is.na(ras@data$v))] <- 0
m <- max(ras@data$v)
s <- m / 10                          
indmax  <- c()
indplan <- c()

for (i in seq(along=ras@data$v))  {      
    x <- coordinates(ras)[i,1]
    y <- coordinates(ras)[i,2]
    z <- ras@data$v[i]
    indx <- which((coordinates(ras)[,1] > x - r) & (coordinates(ras)[,1] < x + r))
    indy <- which((coordinates(ras)[,2] > y - r) & (coordinates(ras)[,2] < y + r))
    indxy <- intersect(indx,indy)
    if (max(ras@data[indxy,1]) == z & z > s) {indmax[length(indmax)+1]   <- i}  
    if (sd(ras@data[indxy,1]) == 0)  {indplan[length(indplan)+1] <- i}
    rm(indx)
    rm(indy)
    rm(indxy)
}

mn <- length(indmax)      
mx <- coordinates(ras)[indmax,1]
my <- coordinates(ras)[indmax,2]
mx2 <- coordinates(ras)[indplan,1]
my2 <- coordinates(ras)[indplan,2]
mz <- ras@data[indmax,1]
maxima <- data.frame(cbind(mx,my,mz))

image(ras, col = gray.colors(25, start = 0.97, end = 0.4))
points(mx2,my2,pch=16, col="gray")              # points in the plane
points(mx,my,pch=16, col="black")               # lokal maxima

library(deldir)
try <- deldir(maxima[,1],maxima[,2],plot=TRUE,wl='te')

cent <- data.frame(cbind(id=seq(1:length(maxima[,1])),x=maxima[,1],y=maxima[,2],meg=0, tum=0))
cent2 <- data.frame(cbind(id=seq(1:length(maxima[,1])),x=maxima[,1]+1000,y=maxima[,2] +1000,meg=0,tum=0))
coordinates(cent)=~x+y
proj4string(cent)  <- CRS(as.character(crs1))
coordinates(cent2)=~x+y
proj4string(cent2)  <- CRS(as.character(crs1))
cent_meg  <- overlay(x=sgdf_meg_dens, y=cent2)
cent_tum  <- overlay(x=sgdf_tum_dens, y=cent)
cent@data$meg <- cent_meg@data$v
cent@data$tum <- cent_tum@data$v
dens_samp@data <- cbind(dens_samp@data,cent=0)

for(i in seq(along=(dens_samp@data$cent))) {
    d1 <- sqrt((dens_samp@data$meg[i] - cent@data$meg[1])^2 + (dens_samp@data$tum[i] 
                                                               - cent@data$tum[1])^2)  
    d2 <- sqrt((dens_samp@data$meg[i] - cent@data$meg[2])^2 + (dens_samp@data$tum[i] 
                                                               - cent@data$tum[2])^2) 
    d3 <- sqrt((dens_samp@data$meg[i] - cent@data$meg[3])^2 + (dens_samp@data$tum[i] 
                                                               - cent@data$tum[3])^2) 
    d4 <- sqrt((dens_samp@data$meg[i] - cent@data$meg[4])^2 + (dens_samp@data$tum[i] 
                                                               - cent@data$tum[4])^2) 
    d <- c(d1,d2,d3,d4)
    mindist <- min(d1,d2,d3,d4)
    id  <- which(d == mindist)
    dens_samp@data$cent[i] <- id
}
plot(dens_samp, pch=dens_samp@data$cent)

library(raster)
library(gdistance)
nz <- 8   # neigbourhood number: 4,8,16
projection(ras) <- crs1
ras <- focal(ras, w=matrix(1/9,nrow=3,ncol=3), NAonly=TRUE)
plot(ras, col = gray.colors(25, start = 0.97, end = 0.4))
starts <- cbind(cent@coords[,1],cent@coords[,2])

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

# auxilliary function
hdiff <- function(x){x[2]-x[1]}

# transitional object 
hd <- transition(ras,hdiff,nz,symm=FALSE)    
slope <- geoCorrection(hd,scl=FALSE)             # transitional object with slope
adj <- adjacent(x=ras, cells=1:ncell(ras), direction=nz) 
cost <- slope       
cost[adj] <- herzog2012_wi(slope[adj]) 
conduct <- geoCorrection(cost, scl=FALSE) # conductivity=cost/dist; time=1/conductivity

cost_surface1 <- accCost(conduct, starts[1,])
cost_surface2 <- accCost(conduct, starts[2,])
cost_surface3 <- accCost(conduct, starts[3,])
cost_surface4 <- accCost(conduct, starts[4,])

#  voronoi in economic space
csm1 <- cost_surface1 == min(cost_surface1,cost_surface2,cost_surface3,cost_surface4)
csm2 <- cost_surface2 == min(cost_surface1,cost_surface2,cost_surface3,cost_surface4)
df <- data.frame(id=c(0,1), v=c(0,2))
csm2 <- subs(x=csm2, y=df)
csm3 <- cost_surface3 == min(cost_surface1,cost_surface2,cost_surface3,cost_surface4)
df <- data.frame(id=c(0,1), v=c(0,3))
csm3 <- subs(csm3, df)
csm4 <- cost_surface4 == min(cost_surface11,cost_surface2,cost_surface3,cost_surface4)
df <- data.frame(id=c(0,1), v=c(0,4))
csm4 <- subs(csm4, df)
csm <- csm1 + csm2 + csm3 + csm4
plot(csm, col = gray.colors(25, start = 0.97, end = 0.4))
points(starts, pch=16, cex=0.8)  

#  weighted voronoi in economic space
cost_surface1b <- cost_surface1 * 2
csm1 <- cost_surface1b == min(cost_surface1b,cost_surface2,cost_surface3,cost_surface4)
csm2 <- cost_surface2 == min(cost_surface1b,cost_surface2,cost_surface3,cost_surface4)
df <- data.frame(id=c(0,1), v=c(0,2))
csm2 <- subs(x=csm2, y=df)
csm3 <- cost_surface3 == min(cost_surface1b,cost_surface2,cost_surface3,cost_surface4)
df <- data.frame(id=c(0,1), v=c(0,3))
csm3 <- subs(csm3, df)
csm4 <- cost_surface4 == min(cost_surface1b,cost_surface2,cost_surface3,cost_surface4)
df <- data.frame(id=c(0,1), v=c(0,4))
csm4 <- subs(csm4, df)
csm <- csm1 + csm2 + csm3 + csm4
plot(csm, col = gray.colors(25, start = 0.97, end = 0.4))
points(starts, pch=16, cex=0.8)  

# start points
library(deldir)
try <- deldir(starts[,1], starts[,2], plot=TRUE,wl='tr')         
LinesList <- list()      
id <- seq(along=try$delsgs[,1])
for(i in seq(along=try$delsgs[,1])) {    
    m <- matrix(data = c(try$delsgs[i,1],try$delsgs[i,3],try$delsgs[i,2],
                         try$delsgs[i,4]), nrow=2, ncol=2)   
    L <- Line(m)
    LL <- list(L)
    name  <- paste("Kante", "_", try$delsgs[i,5],"_", try$delsgs[i,6], sep="")
    LLL <- Lines(LL, ID = name)               
    LinesList[length(LinesList)+1] <- LLL     
}
deldf <- data.frame(try$delsgs[,5], try$delsgs[,1], try$delsgs[,2], try$delsgs[,6], 
                    try$delsgs[,3], try$delsgs[,4], paste("Kante", "_", try$delsgs[,5],"_"), 
                    try$delsgs[,6], sep="")
cols <- c("a","b","c","d","e","f","g","h")
colnames(deldf) <- cols
sl <- SpatialLines(LinesList, proj4string = CRS("+proj=tmerc +lat_0=0 +lon_0=9 +k=1 
                                                +    +x_0=3500000 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
starts2 <- SpatialLinesDataFrame(sl, deldf, match.ID = FALSE)

# least cost path
for(i in 1:length(starts2[,1])){
    i1 <- i*2-1
    i2 <- i*2
    s <- c(starts2@data$b[i],starts2@data$c[i])
    z <- c(starts2@data$e[i],starts2@data$f[i])
    sz <- shortestPath(conduct, s, z, output="SpatialLines")
    zs <- shortestPath(conduct, z, s, output="SpatialLines")
    sz@lines[[1]]@ID <- as.character(i1)    
    zs@lines[[1]]@ID <- as.character(i2)           
    if(i==1){sdf <-rbind(sz,zs)}
    if(i>1){sdf <- rbind(sdf,sz,zs, makeUniqueIDs = TRUE)} 
    if(i==1){df <- cbind(c(1,2), c("sz","zs"), c(starts2@data$g[i]))}
    if(i>1){df <- cbind(c(df[,1],i1,i2), c(df[,2],"sz","zs"), c(df[,3],
                                                                starts2@data$g[i],starts2@data$g[i]))}
}
lcp_df <- as.data.frame(df)
lcp_sldf <- SpatialLinesDataFrame(sdf,lcp_df, match.ID = FALSE) 














co <- coordinates(cent)  
coords=as.matrix(coordinates(cent))
ids <- row.names(as.data.frame(cent))
wts <- fs[,1]; wts[] <- 1
fs_nb_del <- tri2nb(co, row.names=ids)       
del <- nb2lines(fs_nb_del, wts=wts, coords=coords, proj4string = 
                    CRS(as.character(crs1)))
fs_nb_soi <- graph2nb(soi.graph(fs_nb_del, co), row.names=ids)
soi <- nb2lines(fs_nb_soi, wts=wts, coords=coords, proj4string = 
                    CRS(as.character(crs1)))
fs_nb_gabriel <- graph2nb(gabrielneigh(co), row.names=ids)  
gabriel <- nb2lines(fs_nb_gabriel, wts=wts, coords=coords, proj4string = 
                        CRS(as.character(crs1)))
fs_nb_relativ <- graph2nb(relativeneigh(co), row.names=ids) 
relativ <- nb2lines(fs_nb_relativ, wts=wts, coords=coords, proj4string = 
                        CRS(as.character(crs1)))

par(mfrow=c(2,2))
image(sgdf_srtm, col = gray.colors(25, start = 0.97, end = 0.4))  
lines(del)
image(sgdf_srtm, col = gray.colors(25, start = 0.97, end = 0.4))  
lines(soi)
image(sgdf_srtm, col = gray.colors(25, start = 0.97, end = 0.4))  
lines(gabriel)
image(sgdf_srtm, col = gray.colors(25, start = 0.97, end = 0.4))  
lines(relativ)
par(mfrow=c(1,1))

library(rgdal)
fsd <- tri.mesh(spdf_meg, duplicate = 'remove')   # ermittelt die delauny trinangulation
fsnn <- neighbours(fsd)     # Vektor mit Listen der Nachbarn 
LinesList <- list()         # leere LinesList Liste erstellen
sldf <- c();deldf_i <- c();deldf_x1 <- c();deldf_y1 <- c();deldf_k <- c();deldf_x2 <- 
    c();deldf_y2 <- c();deldf_name <- c()
for(i in seq(along=fsg@coords[,1])) {         
    pid1 <- i                                  
    x1 <- fsg@coords[i,1]; y1 <- fsg@coords[i,2]
    for(k in seq(along=(fsnn[i][[1]]))) {   
        pid2 <- fsnn[[i]][k]    
        if (pid2 > pid1) {     
            x2 <- fsg@coords[pid2,1]
            y2 <- fsg@coords[pid2,2] 
            m <- matrix(data = c(x1,x2,y1,y2), nrow=2, ncol=2)  
            L <- Line(m); LL <- list(L)
            name  <- paste("Kante", "_", pid1,"_", pid2, sep="")
            LLL <- Lines(LL, ID = name)             
            LinesList[length(LinesList)+1] <- LLL   
            sldf[length(sldf)+1] <- name         
            j <- length(deldf_i) + 1
            deldf_i[j]    <- i
            deldf_x1[j]   <- x1; deldf_y1[j]   <- y1
            deldf_k[j]    <- pid2
            deldf_x2[j]   <- x2; deldf_y2[j]   <- y2
            deldf_name[j] <- name
        }
    }
}
deldf <- data.frame(deldf_i,deldf_x1,deldf_y1,deldf_k,deldf_x2,deldf_y2,deldf_name)
sldf2 <- data.frame(sldf)
sl <- SpatialLines(LinesList, proj4string = CRS(as.character(crs1)))
sdf <- SpatialLinesDataFrame(sl, sldf2, match.ID = FALSE)
writeOGR(sdf, arbverz, "delaunay_meg", driver="ESRI Shapefile", overwrite_layer=TRUE)
write(rbind(deldf_i,deldf_x1,deldf_y1,deldf_k,deldf_x2,deldf_y2,deldf_name), file 
      = "modproj/5result/delaunay_meg.csv", ncolumns = 7, sep = ";")

nb_k1 <- knn2nb(knearneigh(spdf_meg, k=1),row.names=ids) 
nb_k2 <- knn2nb(knearneigh(spdf_meg, k=2),row.names=ids)   
nb_k3 <- knn2nb(knearneigh(spdf_meg, k=3),row.names=ids)    
nb_k4 <- knn2nb(knearneigh(spdf_meg, k=4),row.names=ids)  

par(mfrow=c(2,2))
plot.nb(nb_k1, fs, pch= 16, col="black", points=TRUE, add=FALSE, arrows=FALSE, cex=0.4)
title("k = 1")
plot.nb(nb_k2, fs,  pch= 16, col="black", points=TRUE, add=FALSE, arrows=FALSE, cex=0.4)
title("k = 2")
plot.nb(nb_k3, fs,  pch= 16, col="black", points=TRUE, add=FALSE, arrows=FALSE, cex=0.4)
title("k = 3")
plot.nb(nb_k4, fs,  pch= 16, col="black", points=TRUE, add=FALSE, arrows=FALSE, cex=0.4)
title("k = 4")
par(mfrow=c(1,1))











































library(igraph)
library(spdep)
co <- coordinates(spdf_meg)  
coords=as.matrix(coordinates(spdf_meg))
ids <- row.names(as.data.frame(spdf_meg))
fs_nb_del <- tri2nb(fs, row.names=ids) 
m <- nb2mat(fs_nb_del)             
g <- graph.adjacency(m, mode="lower", weighted=T)
g <- set.vertex.attribute(g, "x", index=V(g), coordinates(spdf_meg)[,1])   
g <- set.vertex.attribute(g, "y", index=V(g), coordinates(spdf_meg)[,2]) 

dpd        <- deldf
distanzen  <- sqrt((dpd[2] - dpd[5])^2 + (dpd[3] - dpd[6])^2 )
distanzen2 <- distanzen^2
dpd$dist   <- distanzen
dpd$dist2  <- distanzen2
g          <- set.edge.attribute(g, "distanz2", index=E(g), distanzen2)
E(g)$weight <- distanzen

c.degree <- degree(g, v=V(g))
c.closness <- closeness(g, v=V(g))
c.betweenness <- betweenness(g, v=V(g))
c.bonnacich.power <- bonpow(g, nodes=V(g))

ctab <- data.frame(cbind(id = V(g)+1, x = fs[,1], y = fs[,2], degree=c.degree, closness=c.closness, betweenness=c.betweenness, bpower=c.bonnacich.power))
coordinates(ctab)=~x+y  
proj4string(ctab)  <- CRS(as.character(crs1)) 
image(sgdf_srtm, col = gray.colors(25, start = 0.97, end = 0.4))   
points(ctab, pch=16, cex=sqrt(ctab$betweenness)/50)

spoints <- data.frame(cbind(x=(seq(1:44)*500) + 3552000, y=rep(6031270, 44)))
coordinates(spoints)=~x+y
proj4string(spoints)  <- CRS(as.character(crs1)) 

i_kde <- overlay(x=sgdf_meg_dens, y=spoints)
i_kde_tum <- overlay(x=sgdf_tum_dens, y=spoints)
names(i_kde)[names(i_kde) == 'v'] <- 'meg'
i_kde@data    <- cbind(i_kde@data,i_kde_tum@data$v)
names(i_kde)[names(i_kde) == 'i_kde_tum@data$v'] <- 'tum'

distance <- sqrt((i_kde$meg[1] - i_kde$meg)^2 + (i_kde$tum[1] - i_kde$tum)^2)

plot(distance, col="black", pch=16)
lines(distance,lty=2,col="black")






