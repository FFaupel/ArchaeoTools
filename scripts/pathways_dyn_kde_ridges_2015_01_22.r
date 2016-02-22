################################################################################
## R-Script for pathway reconstruction based on monuments
## =============================================================================
## Project: GSHDL Modelling Summer School 2013 'Pathways for Megalithictombs'
## Author: Oliver Nakoinz with contributions from 
##         Arne Windler, Luise Lorenz, Rene Ohlrau, Hendaya Serrano Gil
## Version: 07
## Date of last changes: 22.01.2015
## Data: Megalithic tombs from Schleswig-Holstein
## Author of data: Archäologisches Landesamt Schleswig-Holstein / SPP 1400
## Purpose: empirical pathway model
## Content: 1. preparation, 2. data, 3. grid, 4. static & dynamic KDE, 5. pathway detection, 
##         6. correction, 7. relative neighbour graph, 8. export
## Description: First a static KDE, with high sd is calculated. Then a dynamic 
##       or adaptive KDE, where sd depends on the static KDE is calculated and a 
##       ridge detection applied. Dynamik KDE and ridge detection is iterativ. 
##       the grid-cells of the ridges were added to the points and the next interation is
##       calculated. Finally the relative-neighbour-graph is calculated to connect the fragments
##       of pathways. 
##       The kernel is a cosius with an added treeangular roof in the centre and 
##       1/dist in the remote areas. 
## Licence data: -
## Licence Script: GPL
################################################################################

time1 <- Sys.time()  # start time
Sys.time()
# @@ RGEDIT LANDMARK @@: 1. PREPARATION
########################################################## START parameter input
################################################################################
################################################################################
### 1.1 set working directory and defining file-names
arbverz      <- "/home/fon/daten/analyse/wege_sh2014_meg"    # variable with the path
setwd(arbverz)                                    # set the working directory
avs <- paste(arbverz,"/2data/",sep="")            # working-directory data (for shape-files)
avs4 <- paste(arbverz,"/4result/",sep="")  
avs3 <- paste(arbverz,"/3geodata/",sep="")  
file_meg <- "sh_mega_dg"
#file_meg <- "rpt01a"
file_meg_ex <- "sh_mega_dg2"                          # name for file with export-points for iteration
iter <- 2                                       # steps of iteration

#file_ag  <-  "ag"                                  # Shp-Datei Ausdehunung Arbeitsgebiet
file_srtm  <- "3geodata/sh_90_gk3.tif" # srtm DTM

### 1.2 define values
rw     <- 500   # Rasterweite definieren
tresh  <- 0.05  # Pfade  in Flächen mit sehr niedrigem Dichtewert (KDE), also Rechenartefakte werden unterhalb dieses Schwellenwertes (treshold) entfernt (genullt) 
f_sd1  <- 4     # Faktor für die Größe des ersten Kernels, der die grundlegende Struktur des dynamischen Kernels bildet
f1     <- 0.2      # Faktor für die Untergrenze des dynamischen Kernels (Rasterweite) f1*mean(nn)  ## 0.2
f2     <- 0.4      # Faktor für die Obergrenze des dynamischen Kernels f2*mean(nn)
f3     <- 0.5      # Minimalinentität des Kernels
f4     <- 1        # Maximalinentität des Kernels
s      <- -0.3     # Kernelparameter: Steigung von Punkt 1
de     <- 0.7      # Höhe der zusätzlichen Kernelspitze
sw     <- 12       # Bildbreite in cm
mwin   <- 9        # Mowing-window-size for ridge detection (4,9,16)
#sup    <- 1        # Faktor für die Überlagerung von statischer Dichte und dynamischer Dichte (d = sup*d_statisch + d_dynamisch) 
### area of interest (ag-shapefile, comment if not needed): 
#xmin    <- 765270 
#xmax    <- 813280
#ymin    <- 6013270
#ymax    <- 6071990
#xp <- 750     # Kernelparameter: x-wert   von Punkt 1

### 1.3 define projections
crs1 <- "+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel +towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7 +units=m +no_defs"
crs3 <- "+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0" # Ziel CRS (UTM-Zone 32 N)
crs2 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"                              # Quell CRS (unprojected WGS84) 
################################################################################
################################################################################
########################################################## ENDE  parameter input

### 1.3 load packages
library(rgdal) 
library(sp)
library(spatstat)
library(raster)
library(maptools)
library(spdep)
library(raster)

### 1.4 functions
edist  <- function(x){sqrt((x[1] - x[3])^2 + ((x[2] - x[4])^2))}  # euklidische Distanz
gau1   <- function(x, sd){dnorm(edist(x), mean=0, sd=sd)}         # euklidische Distanzen gewichtet mit normalverteilung des statischen kernels

kernel.par <- function(xp,s,int){
    #Teil 1 der Funktion
    x1 <- asin(-s)     # x-wert mit Steigung s der cosinus-komponente
    y1 <- cos(x1)
    #Teil 2 der Funktion
    x2 <- (-1/s)^0.5     # x-wert mit Steigung s der 1/x-komponente   
    y2 <- 1/x2
    # Parameter zum Anpassen der Komponenten
    a <- y2-y1
    b <- x2-x1
    c <- x1/xp
    #   int: scaling factor for intensity (result*int)
    #
    # Rückgabewerte zusammenfassen
    return(c(xp,a,b,c,int))
    }

kernel1 <- function(x,kp){
    # compound kernel using cos(x) and 1/x for weighting
    #   the distance between points
    # Args:
    #   x: distance between points 
    #   kp: vector made with kernel.par
    # Return:
    #   y: weighted distance
    xp <- kp[1]*kp[4]
    a  <- kp[2]
    b  <- kp[3]
    c  <- kp[4]
    int  <- kp[5]
    x <- x*c
    if (x <= xp) {y <- cos(x) + a + (de-de*x/xp)}     # cosinus mit einem aufgesetzten Dreieck der Höhe de
    if (x > xp)  {y <- 1/(x + b)}
    y <- y*int
    return(y)
    }

kernel1d <- function(x,kp){   # mit Distanzberechnung
    d <- edist(x)
    xp <- kp[1]*kp[4]
    a  <- kp[2]
    b  <- kp[3]
    c  <- kp[4]
    int  <- kp[5]
    x <- d*c
    if (x <= xp) {y <- cos(x) + a + (de-de*x/xp)}     # cosinus mit einem aufgesetzten Dreieck der Höhe de
    if (d > xp)  {y <- 1/(x + b)}
    y <- y*int
    return(y)
    }

# @@ RGEDIT LANDMARK @@: 2. DATA
### 2.1 load data
###### 2.1.1 area of interest
#ag     <- readOGR(avs3, p4s=NULL, file_srtm)
#sgdf_srtm <- read.asciigrid(file_srtm) 
sgdf_srtm <- readGDAL(file_srtm)
xmin    <- sgdf_srtm@bbox[1,1] 
xmax    <- sgdf_srtm@bbox[1,2]
ymin    <- sgdf_srtm@bbox[2,1]
ymax    <- sgdf_srtm@bbox[2,2]

###### 2.1.2 points/monuments
######### 2.1.2a csv
#meg <- read.table(file_meg, header=TRUE, sep=";", dec=".") # Tabelle mit Megalithgräbern auf Ruegen
#coordinates(meg) <- c("koo_rechts", "koo_hoch")            # Spalten für Koordinaten werden definiert
######### 2.1.2b shp
meg      <- readOGR(avs, p4s=NULL, file_meg)
meg <- project(cbind(meg@coords[,1], meg@coords[,2]), crs1)
meg <- SpatialPoints(meg)
projection(meg)  <- crs1
### 2.2 change projection of the coordinate system and coordinates
#proj4string(meg)  <- CRS(as.character(crs2))                              # set original coordinate system
#meg_coord <- project(cbind("x"=meg@coords[,1], "y"=meg@coords[,2]), crs1) # change projection of the data
#df_meg <- data.frame(meg_coord)                                           # create a data frame from the matrix 
#coordinates(df_meg) <- c("X1", "X2")                                      # define coordinates again
#proj4string(df_meg)  <- CRS(as.character(crs1))                           # define projection again

### 2.3 SRTM-Daten laden
# lsgdf_srtm <- readGDAL(file_srtm)

# @@ RGEDIT LANDMARK @@: 3. GRID
### 3.1 Rahmeneckpunkte werden genutzt, das Seitenverhältnis des Ergebnisses zu definieren
sv <- (xmax-xmin)/(ymax-ymin)

### 3.2 Rahmenausdehnung wird definiert und Rahmen wird erstellt
zeilen  <- round((ymax-ymin)/rw, 0) + 1                                     #Anzahl der Zeilen
spalten <- round((xmax-xmin)/rw, 0) + 1                                     #Anzahl der Spalten
v <- cbind(1:(spalten*zeilen))                                              #Vektor für alle Gridzellen
df <- data.frame(v)                                                         #Konvertiert den Vektor zu einem Dataframe
gt <- GridTopology(c(xmin, ymin), c(rw, rw), c(spalten, zeilen))            #Gridtopology wird erstellt (c(linke untere ecke), c(Zellgröße), c(Gridausdehnung) )                       
sgdf <- SpatialGridDataFrame(gt, df, proj4string = CRS(as.character(crs1))) # Grid wird projeziert
#image(sgdf)
#points(meg)

# @@ RGEDIT LANDMARK @@: 4. DYNAMIC KERNEL DENSITY ESTIMATION
    ### für Flächen mit hoher Dichte ist ein Kernel mit einer kleineren rw vorteilhaft, 
    #### für Flächen mit niedriger Dichte ein Kernel mit größerer rw
    ### Berechnung eines Grundlagenkernels um Aussagen über Flächen mit hohen und niedrigen Dichtewert zu treffen 
    #### bei der nächsten Dichteberechnung Berechnung Anpassung der Kernelgröße an die 
    ### Dichtewerte des Grundlagenkernels, Dichtewert des Grundlagenkernels wird als Faktor verwendet

### 4.1 preparation
###### define boundingbox 
bb   <- bbox(sgdf)
win  <- owin(xrange=c(bb[1,1],bb[1,2]), yrange=c(bb[2,1],bb[2,2]), unitname="m")
#ppp_meg <- ppp(df_meg@coords[,1], df_meg@coords[,2], window=win)
ppp_meg <- ppp(meg@coords[,1], meg@coords[,2], window=win)

###### Variablen in Abhängigkeit zur mittleren nächsten Nachbardistanz setzen
nn   <- mean(nndist(ppp_meg))
sd1  <- f_sd1*nn

### 4.2 calculate static-KDE
if(exists("stest") == FALSE) {stest <- 0}
if(stest != sd1) {
    base_kde <- sgdf
    for (i in seq(along=base_kde@data$v)){
        pdist <- cbind(coordinates(base_kde)[i,1],coordinates(base_kde)[i,2],ppp_meg$x[],ppp_meg$y[])     ### x,y,x1,y1
        base_kde@data$v[i]<- sum(apply(pdist,1,gau1,sd=sd1))
        }
    image(base_kde, col = gray.colors(25, start = 0.97, end = 0.4))   
    points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)   
    }
stest <- sd1

for(i  in 1:iter) {       ######################################################## (main loop)
### 4.3 Scaling factor for dynamic KDE 
###### 4.3.1 bandwidth
    d2 <- min(base_kde@data$v)   # d1 minimal density
    d1 <- max(base_kde@data$v)   # d2 maximal density
    a  <- -(d2*f1-d1*f2)/(d1-d2) # a coefficient 1
    b  <- (f1-f2)/(d1-d2)        # b coefficient 2
    factor <- function(x){a+b*x} # factor(base_kde@data$v[i])
    fsd2 <- factor(base_kde@data$v)   # Vektor with sd values for dynamic kernel
    sd2 <- fsd2*nn               # sd = faktor * mean distance to nearest neighbour
###### 4.3.2 intensity (maximal y-value)
    a  <- -(d1*f3-d2*f4)/(d2-d1) # a coefficient 1
    b  <- (f3-f4)/(d2-d1)        # b coefficient 2
    factor_i <- function(x){a+b*x} # factor(base_kde@data$v[i])
    fint <- factor_i(base_kde@data$v)   # Vektor with sd values for dynamic kernel
    
    
### 4.4 dynamic KDE  (in the loop)
    dyn_kde <- sgdf
    date()
    for (i in seq(along=dyn_kde@data$v)){
        sdi <- sd2[i]
        finti <- fint[i]
        kerpar <- kernel.par(sdi,s,finti)
        pdist <- cbind(coordinates(dyn_kde)[i,1],coordinates(dyn_kde)[i,2],ppp_meg$x[],ppp_meg$y[])     ### x,y,x1,y1
        #dyn_kde@data$v[i]<- sum(apply(pdist,1,gau1,sd=sdi))
        dyn_kde@data$v[i]<- sum(apply(pdist,1,kernel1d,kp=kerpar))
    #kernel1d(pdist[22,],kp=kerpar)
    }
    date()
    ##############dyn_kde@data$v <- dyn_kde@data$v + sup * base_kde@data$v    ### groben und dynamischen Kernel kombinieren
    #image(dyn_kde, col = gray.colors(25, start = 0.97, end = 0.4))   
    #points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)   
    
# @@ RGEDIT LANDMARK @@: 5. PATHWAY DETECTION (in the loop)
### 5.1 Ridges Peucker-Douglas-Algorithmus (in the loop)
    ras_ridges <- sgdf
    ras        <- dyn_kde
    ras_ridges@data$v <- 1        # initial value of 1 for ridges
    ras@data$v[is.na(ras@data$v)] <- 10000000000
    
    if (mwin==4){
        for (i in 1:(length(ras@data$v)-spalten))  { 
            ind <- c(i,i+1,i+spalten,i+spalten+1)
            ind_min1 <- which(ras@data$v[ind]==min(ras@data$v[ind]))
            ind_min2 <-ind[ind_min1]
            ras_ridges@data$v[ind_min2] <- 0
            }
        image(ras_ridges, col = gray.colors(25, start = 0.97, end = 0.5))   
        points(ppp_meg, pch=16, cex=0.4)  
        }
    
    if (mwin==9){
        for (i in 1:(length(ras@data$v)-(2*spalten)-2))  { 
            ind <- c(i,i+1,i+2,i+spalten,i+spalten+1,i+spalten+2,i+(2*spalten),i+(2*spalten)+1,i+(2*spalten)+2)
            ind_min1 <- which(ras@data$v[ind]==min(ras@data$v[ind]))
            ind_min2 <-ind[ind_min1]
            ras_ridges@data$v[ind_min2] <- 0
            }
        image(ras_ridges, col = gray.colors(25, start = 0.97, end = 0.5))   
        points(ppp_meg, pch=16, cex=0.4)  
        }
    
    if (mwin==16){
        for (i in 1:(length(ras@data$v)-(3*spalten)-3))  { 
            ind <- c(i,i+1,i+2,i+3,i+spalten,i+spalten+1,i+spalten+2,i+spalten+3,i+(2*spalten),i+(2*spalten)+1,i+(2*spalten)+2,i+(2*spalten)+3,i+(3*spalten),i+(3*spalten)+1,i+(3*spalten)+2,i+(3*spalten)+3)
            ind_min1 <- which(ras@data$v[ind]==min(ras@data$v[ind]))
            ind_min2 <-ind[ind_min1]
            ras_ridges@data$v[ind_min2] <- 0
            }
        image(ras_ridges, col = gray.colors(25, start = 0.97, end = 0.5))   
        points(ppp_meg, pch=16, cex=0.4)  
        }
    
# @@ RGEDIT LANDMARK @@: 6. CORRECTION (in the loop)
        ### Korrektur zu niedriger Dichtewerte und Rechenartefakte (Rauschen) die zu künstlichen Pfaden im berechneten Wegesystem führen
        ### Ridges im Bereich niedriger Dichten werden entfernt
### 6.1 remove noise
    ras_ridges_corr <- ras_ridges
    tresh_dens <- tresh*max(ras@data$v)
    ras_ridges_corr@data$v[which (ras@data$v < tresh_dens)] <- 0
    
### 6.2 sample and combine ridge and monuments points
    coordsnew <- rbind(coordinates(ras_ridges_corr)[which(ras_ridges_corr@data$v==1),], coordinates(meg))
    ddd <- data.frame(id=1:length(coordsnew[,1]))
    pointges = SpatialPointsDataFrame(coordsnew, ddd)
    ppp_meg <- as.ppp(pointges)
} ############################################################# end of main loop

### 6.3 Plot results
#image(ras_ridges_corr, col = gray.colors(25, start = 0.97, end = 0.5)) 
#points(ppp_meg, pch=16, cex=0.4)  

#par(mfrow=c(1,2))
    #image(dyn_kde, col = gray.colors(25, start = 0.97, end = 0.4))  
    #points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)   
    #image(ras_ridges, col = gray.colors(25, start = 0.97, end = 0.5))   
    #points(ppp_meg, pch=16, cex=0.4)   
#par(mfrow=c(1,1))

# @@ RGEDIT LANDMARK @@: 7. Relative Neighbour graph of ridge and points
### 7.1 sample points without monuments
coordsnew <- rbind(coordinates(ras_ridges_corr)[which(ras_ridges_corr@data$v==1),])

### 7.2 calculate rn graph
fs_nb_relativ <- graph2nb(relativeneigh(coordsnew) ) 
wts <- seq(1:length(fs_nb_relativ)); wts[] <- 1
nbsldf_relativ <- nb2lines(fs_nb_relativ, wts, coordsnew, proj4string=CRS(as.character(crs1)))
#plot.nb(fs_nb_relativ, coordsnew, points=FALSE, add=FALSE, arrows=FALSE)

### 7.3 convert rn graph to raster
rras <- raster(dyn_kde)
rn_grid <-  rasterize(nbsldf_relativ, rras)
rn_grid@data@values[which(is.na(rn_grid@data@values)==FALSE)] <- 1
#rn_grid@data@values[which(is.na(rn_grid@data@values)==TRUE)] <- 0
#image(rn_grid)


# @@ RGEDIT LANDMARK @@: 8. EXPORT 
### 8.1 Points
coordsnew <- rbind(coordinates(ras_ridges_corr)[which(ras_ridges_corr@data$v==1),])
ddd <- data.frame(id=1:length(coordsnew[,1]))
pointnew = SpatialPointsDataFrame(coordsnew, ddd)
folder <- paste(arbverz, "/4result", sep = "") 
writeOGR(pointnew, folder, file_meg_ex, "ESRI Shapefile", overwrite=TRUE)

filename <- paste(file_meg_ex, "a", sep = "") 
coordsnew <- rbind(coordinates(ras_ridges_corr)[which(ras_ridges_corr@data$v==1),], coordinates(meg))
ddd <- data.frame(id=1:length(coordsnew[,1]))
pointges = SpatialPointsDataFrame(coordsnew, ddd)
folder <- paste(arbverz, "/4result", sep = "") 
writeOGR(pointnew, folder, filename, "ESRI Shapefile", overwrite=TRUE)

### 8.2 Lines
filename <- paste(arbverz, "/4result/ridges_rn_", as.character(rw), "_", as.character(s), "_", as.character(f1), "_", as.character(f2), "_", as.character(f3), "_", as.character(f4), "_", as.character(f_sd1), "_it", as.character(iter), "_de", as.character(de), sep = "") 
writeLinesShape(nbsldf_relativ, filename)

### 8.3 Esri-Asc-grid
filename <- paste(arbverz, "/4result/ridges_", as.character(rw), "_", as.character(s), "_", as.character(f1), "_", as.character(f2), "_", as.character(f3), "_", as.character(f4), "_", as.character(f_sd1), "_it", as.character(iter), "_de", as.character(de), ".asc", sep = "") 
writeAsciiGrid(ras_ridges_corr, filename, attr = 1, na.value = -999999, dec = ".")

filename <- paste(arbverz, "/4result/kde_", as.character(rw), "_", as.character(s), "_", as.character(f1), "_", as.character(f2), "_", as.character(f3), "_", as.character(f4), "_", as.character(f_sd1), "_it", as.character(iter), "_de", as.character(de), ".asc", sep = "") 
writeAsciiGrid(dyn_kde, filename, attr = 1, na.value = -999999, dec = ".")

filename <- paste(arbverz, "/4result/ridges_rn_", as.character(rw), "_", as.character(s), "_", as.character(f1), "_", as.character(f2), "_", as.character(f3), "_", as.character(f4), "_", as.character(f_sd1), "_it", as.character(iter), "_de", as.character(de), ".asc", sep = "") 
writeRaster(rn_grid, filename, "ascii", overwrite=TRUE)

### 8.4 png-picture
filename <- paste(arbverz, "/5pictures/ridges_", as.character(rw), "_", as.character(s), "_", as.character(f1), "_", as.character(f2), "_", as.character(f3), "_", as.character(f4), "_", as.character(f_sd1), "_it", as.character(iter), "_de", as.character(de), ".png", sep = "") 
png(filename, height=sw/sv, width=sw, units="cm", res=300, bg = "white") 
    image(ras_ridges_corr, col = gray.colors(25, start = 0.97, end = 0.5))   
    #points(ppp_meg, pch=16, cex=0.4)  
    dev.off()

filename <- paste(arbverz, "/5pictures/ridges_rn_", as.character(rw), "_", as.character(s), "_", as.character(f1), "_", as.character(f2), "_", as.character(f3), "_", as.character(f4), "_", as.character(f_sd1), "_it", as.character(iter), "_de", as.character(de), ".png", sep = "") 
png(filename, height=sw/sv, width=sw, units="cm", res=500, bg = "white") 
    image(rn_grid, col = gray.colors(25, start = 0.97, end = 0.5))   
    points(meg, pch=16, cex=0.03)  
    dev.off()

### 8.5 Analysemetadaten
time2 <- Sys.time()  # end of calculation
zi <- (time2-time1)  # time of calculation

filename <- paste(arbverz, "/4result/ridges_", as.character(rw), "_", as.character(s), "_", as.character(f1), "_", as.character(f2), "_", as.character(f3), "_", as.character(f4), "_", as.character(f_sd1), "_it", as.character(iter), ".txt", sep = "")  # ridges_rw_tresh_fsd1_fsd2min_fsd2max_sup.txt

categ <- c("file      ",
         "rw        ",
         "xmin      ",
         "xmax      ",
         "ymin      ",
         "ymax      ",
         "gridcells ",
         "points    ",
         "tresh     ",
         "f_sd1     ",
         "f1        ",
         "f2        ",
         "f3        ",
         "f4        ",
         "nn        ",
         "s         ",
         "de        ",
         "time      ",
         "time unit "
        )
val <- c(file_meg,
         rw,
         xmin,
         xmax,
         ymin,
         ymax,
         length(ras_ridges@data$v),
         ppp_meg$n,
         tresh,
         f_sd1,
         f1,
         f2,
         f3,
         f4,
         nn,
         s,
         de,
         zi,
         as.character(attributes(zi)[1])
        )

write(rbind(categ,val), file = filename, ncolumns = 2, sep = "\t")

save.image("7ws/ws3.rws")  
#load("7ws/ws3.rws") 

Sys.time()
zi
