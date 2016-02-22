# version 2 ist noch nicht getestet                                             ##### testen


################################################################################
## wvoro
## =====
## R-Skript weighted voronoi analysis
## =============================================================================
## Projekt: congress Cambridge 2013
## Skriptbearbeiter: Oliver Nakoinz
## Version: 02
## letzte Bearbeitung: 13.09.2013
## Daten: SRTM DGM; Fürstensitze
## Datenbearbeiter: O. Nakoinz
## Zweck: Veranschaulichung des gewichteten Voronoi-Graphen
## Inhalt:  0. Vorbereitung, 1. Daten laden etc., 2. vworo, 
##          3. Export
## Lizenz Daten: -
## Lizenz Skript: GPL (http://www.gnu.de/documents/gpl.en.html)
################################################################################


# @@ RGEDIT LANDMARK @@: 0. Vorbereitung =======================================
# 0.1 Variablen setzen ---------------------------------------------------------
arbverz <- "/home/fon/daten/aktuell/analyse/bw_kulturdist/analyse_2012" 
setwd(arbverz)
crs1 <- "+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=WGS84 +units=m +no_defs" # gk3
#grau <- gray.colors(12, 0.55, 0.95) 
#grau <- gray.colors(12, 0.95, 0.55) 
grau <- gray.colors(12, 0.95, 0.40) 
nz <- 8   # Nachbarschaftszahl für LCP: 4,8,16
breite  <- 17  # breite der abb: in cm pz 2 spalten= 17 cm; pz 1 spalten= 8.3 cm
breitei <- 6.6929134   # breite der abb: in cm pz 2 spalten= 6.6929134 inch; pz 1 spalten= 2.7559055 inch

breite2  <- 8.3  # breite der abb: in cm pz 2 spalten= 17 cm; pz 1 spalten= 8.3 cm
breitei2 <- 3.2677   # breite der abb: in cm pz 2 spalten= 5.511811 inch; pz 1 spalten= 2.7559055 inch

#  Gewichte                                                                     ##### Gewichte am besten mit den Zentren einlesen
#Würzburg, Bad Dürkheim, Leutenheim, Asperg, Ipf, Hundersingen, Villingen, Breisach
k <- -5000 # Gewichtungsfaktor der additiven Gewichte
daw <- k * c(1,2,1,1,1,3,1,3)  # Anzahl der zentralen Funktionen nach Posluschny, leicht verändert # dynamische additives Gewicht
dmw <- 1/(c(1,1,1,3,1,3,2,3)^1) # Anzahl von Wagen/Pferd-Fundstellen in 25 km Radius um die FS       # dynamische multiplikatives Gewicht

##Graphikparameter
par(
    ps=10,
    family='sans',
    mar=c(5, 4, 4, 2)   # numerical vector indicating margin size c(bottom, left, top, right) in lines. default = c(5, 4, 4, 2) + 0.1
    )   
# 0.2 Pakete -------------------------------------------------------------------
library(sp)               
library(proj4)              
library(rgdal)              
library(maptools)
library(raster)
library(gdistance)
library(gstat)
library(colorRamps)
library(RColorBrewer)
library(spdep)

# 0.3 Funktionen ---------------------------------------------------------------
edist <- function(a,b){sqrt(sum((a-b) ^ 2))}
herzog2012_wi <- function(s){1/(1337.8 * s^6 + 278.19 * s^5 - 517.39 * s^4 - 78.199 * s^3 + 93.419 * s^2 + 19.825 * s + 1.64)}
hdiff <- function(x){x[2]-x[1]}
minid <- function(x){which(x == min(x))[1]    }

# @@ RGEDIT LANDMARK @@: 1. Daten laden etc. ===================================
time1 <- Sys.time()  # Startzeit
# 1.1 Raster laden -------------------------------------------------------------
ras <- raster("srtm_250_gk3.asc")
projection(ras) <- crs1
ras <- focal(ras, w=matrix(1/9,nrow=3,ncol=3), NAonly=TRUE)

# 1.2 Shapes Raster laden ------------------------------------------------------
riv_l <- readOGR(arbverz, p4s=NULL, 'fluss_gr')
riv_s <- readOGR(arbverz, p4s=NULL, 'fluss_kl')
lake <- readOGR(arbverz, p4s=NULL, 'See')
fs1 <- readOGR(arbverz, p4s=NULL, 'fs_ag')
#fsid <- c(3,8,9,10,11,12)
fsid <- c(2,3,6,7,8,9,10,11)
fs <- coordinates(fs1)[fsid,]
starts <- cbind(fs[,1],fs[,2])

# 1.4 Karte zeichnen -----------------------------------------------------------
image(ras,col=grau)
lines(riv_l,col='grey70',lwd=1.3)
lines(riv_s,col='grey70',lwd=0.5)
points(fs1, cex=5,col="red")
#text(fs, cex=5,col="red")
plot(lake,col='grey80', add=T,border='grey70')

# @@ RGEDIT LANDMARK @@: 2. cost surface =======================================
# 2.1 transitional object 
hd <- transition(ras,hdiff,nz,symm=FALSE)    
slope <- geoCorrection(hd,scl=FALSE)             # transitional object with slope
adj <- adjacent(x=ras, cells=1:ncell(ras), direction=nz) 
cost <- slope       
cost[adj] <- herzog2012_wi(slope[adj]) 
conduct <- geoCorrection(cost, scl=FALSE) # conductivity=cost/dist; time=1/conductivity

save.image("wvoro2013-1.rws")

# @@ RGEDIT LANDMARK @@: 3. Schleifen ==========================================
# 2.2 akkumulierte Kosten
cost_surface <- list()
for (i in seq(along=coordinates(fs1)[,1])){
    
    cost_surface[[i]] <- accCost(conduct, starts[i,])
    }    

# 2.2 Gewichtung
cost_surfaceb <- cost_surface
for (i in seq(along=coordinates(fs1)[,1])){
    cost_surfaceb[[i]] <- cost_surface[[i]] * dmw[i] + daw[i]
    }

# 2.3 Auswahl
csmx <- cost_surface 
for (i in seq(along=coordinates(fs1)[,1])){
    csmx[[i]] <- cost_surfaceb[[i]] == min(cost_surfaceb[[]])                   ##### prüfen ob das geht
    df <- data.frame(id=c(0,1), v=c(0,i))
    csmx[[i]] <- subs(x=csmx[[i]], y=df)
    }

# 2.4 Zusammenführen
csm <- sum(csmx[[]])                                                            ##### prüfen ob das geht


plot(csm, col = gray.colors(8, start = 0.97, end = 0.4))
lines(csm_po)
points(starts, pch=16, cex=0.8)  

# @@ RGEDIT LANDMARK @@: 4. Export =============================================
#csm_pi <- as(csm, "SpatialPixels") 
#csm_po <- as(csm_pi, "SpatialPolygons") 

writeRaster(csm, '/home/fon/daten/bilder/qA/grenzen/wvoro_fs_2013-06.asc', 'ascii', overwrite=TRUE)

# Reimport der Grenzen ### Raster in Vektor verwandeln (QGIS)
v1 <- readOGR(arbverz, p4s=NULL, 'voro1')
v2 <- readOGR(arbverz, p4s=NULL, 'voro2')
v3 <- readOGR(arbverz, p4s=NULL, 'voro3')

##oder 
# csm_pi <- as(csm, "SpatialPixels)
# csm_po <- as(csm_pi, "SpatialPolygons") 


# Plot
image(ras,col=grau)
#lines(riv_l,col='grey70',lwd=1.3)
#lines(riv_s,col='grey70',lwd=0.5)
points(fs1, pch=16,cex=2,col="black")
plot(lake,col='grey80', add=T,border='grey70')
plot(v2,border="orange",lwd=2.3,add=T)
plot(v1,border="red",lwd=1.5,add=T)
plot(v3,border="blue",lwd=1,add=T)



png("/home/fon/daten/bilder/qA/grenzen/wvoro_fs_2013-x.png", height=1.2*breite, width=breite, units="cm", res=300, bg = "white") 
    image(ras,col=grau)
    lines(riv_l,col='grey70',lwd=1.3)
    lines(riv_s,col='grey70',lwd=0.5)
    plot(lake,col='grey80', add=T,border='grey70')
    plot(v2,border="orange",lwd=4,add=T,lty=3)
    plot(v1,border="red",lwd=3,add=T,lty=3)
    plot(v3,border="blue",lwd=2,add=T,lty=3)
    points(fs1, pch=16,cex=2,col="black")
dev.off()    

save.image("wvoro2013-2.rws")




time3 <- Sys.time()  # Endzeit
zi2 <- (time3-time2)  # Dauer der Berechnung
zi2
zi



