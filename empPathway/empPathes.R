#################################################################################
## R-skript for Testanalysis Modelling Empiric Pathways in the Upper Rhine Valley
## ==============================================================================
## Project: IMODEZ
##
####
######   ALL SITES, BUT NOTHING LATER THAN Lt B 
####
##
## Author:          Oliver Nakoinz, with contributions from Franziska Faupel
## Version:         02
## Data:            ArkeoGIS, SHKR
## Author of data:  see ArkeoGIS, SHKR
## Purpose:         empirical pathway model
## Content:         1. SETTING AREA OF INTEREST AND GRID DATA,  2. DYNAMIC KERNEL DENSITY ESTIMATION, 3. PATHWAY DETECTION (in the loop), 
##                  4. CORRECTION (in the loop), 5. RELATIVE NEIGHBOUR GRAPH OF RIDGE AND POINTS, 6. EXPORT
##Description:  EMPIRICAL PATHWAY MODEL
##                  First a static KDE, with high sd is calculated. Then a dynamic 
##                  or adaptive KDE, where sd depends on the static KDE is calculated and a 
##                 ridge detection applied. Dynamik KDE and ridge detection is iterativ. 
##                  the grid-cells of the ridges were added to the points and the next interation is
##                  calculated. Finally the relative-neighbour-graph is calculated to connect the fragments
##                  of pathways. The kernel is a cosius with an added treeangular roof in the centre and 
##                  1/dist in the remote areas.
##
## Date of last changes:  25.11.2015
## Licence data:          -
## Licence Script:        
## Calculation time estimated: 50h
## Intended use: IMODEZ
##
## R Version needed: R.3.1.1
################################################################################

################################################################################
##################################################
################################
####################
#########
######
##
#          GETTING STARTED
##
######
#########
####################
################################
##################################################
################################################################################

time1 <- Sys.time()  # start time
Sys.time()
date <- Sys.time()

# Prepare packages
#--------------------------------------------------------------------------------------------------------------
#install.packages("stats"); install.packages('sp'); install.packages('proj4'); install.packages('rgdal'); install.packages('spatstat');install.packages('RSQLite');install.packages('gdata'); install.packages('KernSmooth'); install.packages('maptools');install.packages('deldir'); install.packages('graphics'); install.packages('gdistance'); install.packages('reaster');install.packages('spdep');
library(stats)
library(sp)
library(proj4)
library(rgdal)
library(spatstat)
library(gdata)
library(KernSmooth)
library(maptools)
library(plyr)
library(raster)
library(deldir)
library(graphics)
library(gdistance)
library(spdep)

# Working Directory
#---------------------------------------------------------------------------------------------------------------
arbverz <- "C:\\Diss\\Pathway modelling\\Pathways Alsace"
setwd(arbverz) 


# Functions
#----------------------------------------------------------------------------------------------------------------
# projection
crs1 <- "+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel +towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7 +units=m +no_defs"
crs3 <- "+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0" # Ziel CRS (UTM-Zone 32 N)
crs2 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"                              # Quell CRS (unprojected WGS84) 

# distance functions
edist  <- function(x){sqrt((x[1] - x[3])^2 + ((x[2] - x[4])^2))}  # euklidische Distanz
gau1   <- function(x, sd){dnorm(edist(x), mean=0, sd=sd)}         # euklidische Distanzen gewichtet mit normalverteilung des statischen kernels

# KDE
# 1. KDE Function ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
kernel.par <- function(xp,s,int){
  # Teil 1 der Funktion
  x1 <- asin(-s)     # x-wert mit Steigung s der cosinus-komponente
  y1 <- cos(x1)
  # Teil 2 der Funktion
  x2 <- (-1/s)^0.5     # x-wert mit Steigung s der 1/x-komponente   
  y2 <- 1/x2
  # Parameter zum Anpassen der Komponenten
  a <- y2-y1
  b <- x2-x1
  c <- x1/xp
  #   int: scaling factor for intensity (result*int)  # Rückgabewerte zusammenfassen
  return(c(xp,a,b,c,int))
}

# 2. KDE Function ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
  if (x <= xp) {y <- cos(x) + a + (de-de*x/xp)}     # cosinus mit einem aufgesetzten Dreieck der HÃÂ¶he de
  if (x > xp)  {y <- 1/(x + b)}
  y <- y*int
  return(y)
}

# 3. KDE Function ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
kernel1d <- function(x,kp){   # mit Distanzberechnung
  d <- edist(x)
  xp <- kp[1]*kp[4]
  a  <- kp[2]
  b  <- kp[3]
  c  <- kp[4]
  int  <- kp[5]
  x <- d*c
  if (x <= xp) {y <- cos(x) + a + (de-de*x/xp)}     # cosinus mit einem aufgesetzten Dreieck der HÃÂ¶he de
  if (d > xp)  {y <- 1/(x + b)}
  y <- y*int
  return(y)
}

# define values
#---------------------------------------------------------------------------------------------------------------
rw     <- 500   # width of raster defined in m
tresh  <- 0.05  # treshold, to delete phaths in areas with really low density values (KDE), meaning calculation artefacts 
f_sd1  <- 4     # factor defining size if the first kernel, which generate the stucture of dynamic kernel
f1     <- 0.2      # factor defining the minimum border of dynamic kernel (raster width) f1*mean(nn)  ## 0.2
f2     <- 0.4      # factor defining the maximum border of dynamic kernel f2*mean(nn)
f3     <- 0.5      # MinimalinentitÃÂ¤t des Kernels
f4     <- 1        # MaximalinentitÃÂ¤t des Kernels
s      <- -0.3     # Kernelparameter: incline starting from ponit 1
de     <- 0.7      # hight of additional kernell peak
sw     <- 12       # width of picture, cm
mwin   <- 9        # Mowing-window-size for ridge detection (4,9,16)
#sup    <- 1        # Faktor fÃÂ¼r die ÃÂberlagerung von statischer Dichte und dynamischer Dichte (d = sup*d_statisch + d_dynamisch) 
xp <- 750     # Kernelparameter: x-wert   von Punkt 1

# Import Data
#--------------------------------------------------------------------------------------------------------------
avs4 <- paste(arbverz,"/4result/",sep="")  
avs3 <- paste(arbverz,"/3geodata/",sep="")  
file_g <- "2data/g.csv"                              # latest change 3677 sites (from Bibracte: data resolving)

g <- read.table(file_g, header=TRUE, sep=';')
spdf_g <- g
coordinates(spdf_g)=~coords.x1+coords.x2
spdf_g <- SpatialPoints(spdf_g)
projection(spdf_g)  <- crs1


#file_ag  <-  "ag"                                  # Shp-Datei Ausdehunung Arbeitsgebiet
file_srtm  <- "3geodata/bw_srtm.asc" # srtm DTM
sgdf_srtm <- read.asciigrid(file_srtm) 
lsgdf_srtm <- readGDAL(file_srtm)

file_g_ex <- "2data/g_2"                          # name for file with export-points for iteration
iter <- 2                                       # steps of iteration

#============================================================================================================
save.image("7ws/ws_emp_1.rws")
#============================================================================================================

################################################################################
##################################################
################################
####################
#########
######
##
#          1. Setting Area of Interet and Grid Data
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_emp_1.rws")

# area of interest                                      # Rheinebene bis Bodensee, Schwäbische Alb bis zum nördl. Ende des Schwarz Waldes
#-----------------------------------------------------------------------------------------------
xmin    <- 3354223
xmax    <- 3587393
ymin    <- 5260128
ymax    <- 5448006
ext_ai <- extent(3354223, 3587393, 5260128, 5448006)
raster_srtm <- raster(sgdf_srtm)
sgdf_ai <- crop(raster_srtm, ext_ai, snap='near', overwrite=TRUE)
sgdf_ai <- as(sgdf_ai, 'SpatialGridDataFrame')
ras_ai <- crop(raster_srtm, ext_ai, filename="raster_ai", snap='near', overwrite=TRUE)

# Coordinates to set frame corner to define aspect ratio                      # Rahmeneckpunkte werden genutzt, das Seitenverhältnis des Ergebnisses zu definieren
sv <- (xmax-xmin)/(ymax-ymin)

# Defintion of frame expantion and defining frame                              # Rahmenausdehnung wird definiert und Rahmen wird erstellt
rows  <- round((ymax-ymin)/rw, 0) + 1 ; rows                                    #Anzahl der Zeilen
colums <- round((xmax-xmin)/rw, 0) + 1 ; colums                                     #Anzahl der Spalten
v <- cbind(1:(colums*rows))                                              #Vektor fÃÂ¼r alle Gridzellen
df <- data.frame(v)                                                         #Konvertiert den Vektor zu einem Dataframe
gt <- GridTopology(c(xmin, ymin), c(rw, rw), c(colums, rows))            #Gridtopology wird erstellt (c(linke untere ecke), c(ZellgrÃÂ¶ÃÂe), c(Gridausdehnung) )                       
sgdf <- SpatialGridDataFrame(gt, df, proj4string = CRS(as.character(crs1))) # Grid wird projeziert
image(sgdf)
points(spdf_g)


#============================================================================================================
save.image("7ws/ws_emp_2.rws")
#============================================================================================================

################################################################################
##################################################
################################
####################
#########
######
##
#          2. DYNAMIC KERNEL DENSITY ESTIMATION
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_emp_2.rws")

### areas with low densities: smal rw kernel is advantageous 
#### areas with high densities: big rw kernel is advantageous 
### calculating a test kernel, to estimate low or high densities in areas 
#### next density estimation, modulation of kernel size to density values of test kernel  
### density values of test kernel will be used as a factor

# define boundingbox 
#----------------------------------------------------------------------------------------------------------------------
bb   <- bbox(sgdf)
win  <- owin(xrange=c(bb[1,1],bb[1,2]), yrange=c(bb[2,1],bb[2,2]), unitname="m")
spdf_g <- remove.duplicates(spdf_g, zero=0, remove.second=TRUE)
ppp_g <- ppp(spdf_g@coords[,1], spdf_g@coords[,2], window=win)
g_ai <- data.frame(ppp_g)


# Variablen in Abhä¤ngigkeit zur mittleren nächsten Nachbardistanz setzen
nn   <- mean(nndist(ppp_g))
sd1  <- f_sd1*nn

# 2. calculate static-KDE
#--------------------------------------------------------------------------------------------------------------------------
if(exists("stest") == FALSE) {stest <- 0}
if(stest != sd1) {
  base_kde <- sgdf
  for (i in seq(along=base_kde@data$v)){
    pdist <- cbind(coordinates(base_kde)[i,1],coordinates(base_kde)[i,2],ppp_g$x[],ppp_g$y[])     ### x,y,x1,y1
    base_kde@data$v[i]<- sum(apply(pdist,1,gau1,sd=sd1))
  }
  image(base_kde, col = gray.colors(25, start = 0.97, end = 0.4))   
  points(ppp_g$x, ppp_g$y, pch=16, cex=0.4)   
}
stest <- sd1

for(i  in 1:iter) {       ######################################################## (main loop)
  
  # Scaling factor for dynamic KDE 
  ## bandwidth---------------------------------------------------------------------------------------------------
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
  
  
  #  dynamic KDE  (in the loop)
  ##------------------------------------------------------------------------------------------------------------------------
  dyn_kde <- sgdf
  date()
  for (i in seq(along=dyn_kde@data$v)){
    sdi <- sd2[i]
    finti <- fint[i]
    kerpar <- kernel.par(sdi,s,finti)
    pdist <- cbind(coordinates(dyn_kde)[i,1],coordinates(dyn_kde)[i,2],ppp_g$x[],ppp_g$y[])     ### x,y,x1,y1
    #dyn_kde@data$v[i]<- sum(apply(pdist,1,gau1,sd=sdi))
    dyn_kde@data$v[i]<- sum(apply(pdist,1,kernel1d,kp=kerpar))
    #kernel1d(pdist[22,],kp=kerpar)
  }
  date()
  
  ################################################################################
  ##################################################
  ################################
  ####################
  #########
  ######
  ##
  #          3. PATHWAY DETECTION (in the loop)
  ##
  ######
  #########
  ####################
  ################################
  ##################################################
  ################################################################################
  
  # 3.1 Ridges Peucker-Douglas-Algorithmus (in the loop)
  #------------------------------------------------------------------------------------------------------------------------------------------------
  ras_ridges <- sgdf
  ras        <- dyn_kde
  ras_ridges@data$v <- 1        # initial value of 1 for ridges
  ras@data$v[is.na(ras@data$v)] <- 10000000000
  
  if (mwin==4){
    for (i in 1:(length(ras@data$v)-colums))  { 
      ind <- c(i,i+1,i+colums,i+colums+1)
      ind_min1 <- which(ras@data$v[ind]==min(ras@data$v[ind]))
      ind_min2 <-ind[ind_min1]
      ras_ridges@data$v[ind_min2] <- 0
    }
    image(ras_ridges, col = gray.colors(25, start = 0.97, end = 0.5))   
    points(ppp_g, pch=16, cex=0.4)  
  }
  
  if (mwin==9){
    for (i in 1:(length(ras@data$v)-(2*colums)-2))  { 
      ind <- c(i,i+1,i+2,i+colums,i+colums+1,i+colums+2,i+(2*colums),i+(2*colums)+1,i+(2*colums)+2)
      ind_min1 <- which(ras@data$v[ind]==min(ras@data$v[ind]))
      ind_min2 <-ind[ind_min1]
      ras_ridges@data$v[ind_min2] <- 0
    }
    image(ras_ridges, col = gray.colors(25, start = 0.97, end = 0.5))   
    points(ppp_g, pch=16, cex=0.4)  
  }
  
  if (mwin==16){
    for (i in 1:(length(ras@data$v)-(3*colums)-3))  { 
      ind <- c(i,i+1,i+2,i+3,i+colums,i+colums+1,i+colums+2,i+colums+3,i+(2*colums),i+(2*colums)+1,i+(2*colums)+2,i+(2*colums)+3,i+(3*colums),i+(3*colums)+1,i+(3*colums)+2,i+(3*colums)+3)
      ind_min1 <- which(ras@data$v[ind]==min(ras@data$v[ind]))
      ind_min2 <-ind[ind_min1]
      ras_ridges@data$v[ind_min2] <- 0
    }
    image(ras_ridges, col = gray.colors(25, start = 0.97, end = 0.5))   
    points(ppp_g, pch=16, cex=0.4)  
  }
  
  ################################################################################
  ##################################################
  ################################
  ####################
  #########
  ######
  ##
  #          4. CORRECTION (in the loop)
  ##
  ######
  #########
  ####################
  ################################
  ##################################################
  ################################################################################
  
  
  ### Korrektur zu niedriger Dichtewerte und Rechenartefakte (Rauschen) die zu kÃÂ¼nstlichen Pfaden im berechneten Wegesystem fÃÂ¼hren
  ### Ridges im Bereich niedriger Dichten werden entfernt
  # 4.1 remove noise
  ##------------------------------------------------------------------------------------------------------------------------
  ras_ridges_corr <- ras_ridges
  tresh_dens <- tresh*max(ras@data$v)
  ras_ridges_corr@data$v[which (ras@data$v < tresh_dens)] <- 0
  
  # 4.2 sample and combine ridge and monuments points
  #------------------------------------------------------------------------------------------------------------------------
  coordsnew <- rbind(coordinates(ras_ridges_corr)[which(ras_ridges_corr@data$v==1),], coordinates(spdf_g)); ddd <- (id=1:length(coordsnew[,1])); coordsnew <- data.frame(coordsnew)
  coordinates(coordsnew)=~s1+s2; ddd <- data.frame(ddd)
  pointges = SpatialPointsDataFrame(coordsnew, ddd)
  ppp_g <- as.ppp(pointges)
} ############################################################# end of main loop

# 4.3 Plot results
#----------------------------------------------------------------------------------------------------------------------------
image(ras_ridges_corr, col = gray.colors(25, start = 0.97, end = 0.5)) 
points(ppp_g, pch=16, cex=0.4)  

par(mfrow=c(1,2))
image(dyn_kde, col = gray.colors(25, start = 0.97, end = 0.4))  
points(ppp_g$x, ppp_g$y, pch=16, cex=0.4)   
image(ras_ridges, col = gray.colors(25, start = 0.97, end = 0.5))   
points(ppp_g, pch=16, cex=0.4)   
par(mfrow=c(1,1))

#============================================================================================================
save.image("7ws/ws_emp_3.rws")
#============================================================================================================
################################################################################
##################################################
################################
####################
#########
######
##
#          5. Relative Neighbour graph of ridge and points
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_emp_3.rws")

# 5.1 sample points without monuments
#--------------------------------------------------------------------------------------------------------------------------
coordsnew <- rbind(coordinates(ras_ridges_corr)[which(ras_ridges_corr@data$v==1),])
coord_rid <- coordsnew

# 5.2 calculate rn graph
#--------------------------------------------------------------------------------------------------------------------------
fs_nb_relativ <- graph2nb(relativeneigh(coordsnew) ) 
wts <- seq(1:length(fs_nb_relativ)); wts[] <- 1
nbsldf_relativ <- nb2lines(fs_nb_relativ, wts, coordsnew, proj4string=CRS(as.character(crs1)))
plot.nb(fs_nb_relativ, coordsnew, points=FALSE, add=FALSE, arrows=FALSE)

# 5.3 convert rn graph to raster
#--------------------------------------------------------------------------------------------------------------------------
rras <- raster(dyn_kde)
rn_grid <-  rasterize(nbsldf_relativ, rras)
rn_grid@data@values[which(is.na(rn_grid@data@values)==FALSE)] <- 1
#rn_grid@data@values[which(is.na(rn_grid@data@values)==TRUE)] <- 0
#image(rn_grid)

# 5.4 convert to proper graphes
#--------------------------------------------------------------------------------------------------------------------------
df <- coord_rid
id <- (id=1:length(coord_rid[,1]))
df <- data.frame(df, id)
coordinates(df)=~s1+s2
df <- remove.duplicates(df, zero=1000, remove.second=TRUE)                                    # width of kde 500, so the radius needs to be 1000

ddd <- df@data$id                                                                             #extract id to reconstruct original coordinates
rid_cle <- coordsnew # loading in coordinates data
rid_cle <- cbind(rid_cle, id) # attaching to id (is needed to detact)
rid_cle <- rid_cle[rid_cle[,3] %in% ddd,] # deleting doubled points, which have been identified before in df
rid_cle <- rid_cle[,c(1:2)]                                                                 #cleaned from unnecessary id vlaues

nb_cle <- graph2nb(relativeneigh(rid_cle)) 
wts <- seq(1:length(nb_cle)); wts[] <- 1
nbl_cle <- nb2lines(nb_cle, wts, rid_cle, proj4string=CRS(as.character(crs1)))
plot(nbl_cle)

#============================================================================================================
save.image("7ws/ws_emp_4.rws")
#============================================================================================================
################################################################################
##################################################
################################
####################
#########
######
##
#          6. EXPORT
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_emp_4.rws")

# 6.1 Points
#--------------------------------------------------------------------------------------------------------------------------
coordsnew <- rbind(coordinates(ras_ridges_corr)[which(ras_ridges_corr@data$v==1),])
ddd <- data.frame(id=1:length(coordsnew[,1]))
pointnew = SpatialPointsDataFrame(coordsnew, ddd)
folder <- paste(arbverz, "/4result", sep = "") 
writeOGR(pointnew, folder, file_g_ex, layer="coord_den_ridges", driver="ESRI Shapefile", overwrite=TRUE)

filename <- paste(file_g_ex, "a", sep = "") 
coordsnew <- rbind(coordinates(ras_ridges_corr)[which(ras_ridges_corr@data$v==1),], coordinates(spdf_g));  ddd <- (id=1:length(coordsnew[,1])); coordsnew <- data.frame(coordsnew)
coordinates(coordsnew)=~s1+s2;
ddd <- data.frame(ddd)
pointges = SpatialPointsDataFrame(coordsnew, ddd)
folder <- paste(arbverz, "/4result", sep = "") 
writeOGR(pointges, folder, filename, layer="g_2a", driver="ESRI Shapefile", overwrite=TRUE)

coordsnew <- rbind(coordinates(rid_cle))
ddd <- data.frame(id=1:length(coordsnew[,1]))
pointnew = SpatialPointsDataFrame(coordsnew, ddd)
folder <- paste(arbverz, "/4result", sep = "") 
writeOGR(pointnew, folder, file_g_ex, layer="coord_den_cleaned", driver="ESRI Shapefile", overwrite=TRUE)

write.table(g_ai, file ="2data/g_ai.csv", sep= ";")

# 6.2 Lines
##--------------------------------------------------------------------------------------------------------------------------
filename <- paste(arbverz, "/4result/ridges_rn_", as.character(rw), "_", as.character(s), "_", as.character(f1), "_", as.character(f2), "_", as.character(f3), "_", as.character(f4), "_", as.character(f_sd1), "_it", as.character(iter), "_de", as.character(de), sep = "") 
writeLinesShape(nbsldf_relativ, filename)

filename <- paste(arbverz, "/4result/ridges_rn_cleaned", sep = "")
writeLinesShape(nbl_cle, filename)

# 6.3 Esri-Asc-grid
#--------------------------------------------------------------------------------------------------------------------------
filename <- paste(arbverz, "/4result/ridges_", as.character(rw), "_", as.character(s), "_", as.character(f1), "_", as.character(f2), "_", as.character(f3), "_", as.character(f4), "_", as.character(f_sd1), "_it", as.character(iter), "_de", as.character(de), ".asc", sep = "") 
writeAsciiGrid(ras_ridges_corr, filename, attr = 1, na.value = -999999, dec = ".")

filename <- paste(arbverz, "/4result/kde_", as.character(rw), "_", as.character(s), "_", as.character(f1), "_", as.character(f2), "_", as.character(f3), "_", as.character(f4), "_", as.character(f_sd1), "_it", as.character(iter), "_de", as.character(de), ".asc", sep = "") 
writeAsciiGrid(dyn_kde, filename, attr = 1, na.value = -999999, dec = ".")

filename <- paste(arbverz, "/4result/ridges_rn_", as.character(rw), "_", as.character(s), "_", as.character(f1), "_", as.character(f2), "_", as.character(f3), "_", as.character(f4), "_", as.character(f_sd1), "_it", as.character(iter), "_de", as.character(de), ".asc", sep = "") 
writeRaster(rn_grid, filename, "ascii", overwrite=TRUE)

# 6.4 png-picture
#--------------------------------------------------------------------------------------------------------------------------
filename <- paste(arbverz, "/5pictures/ridges_", as.character(rw), "_", as.character(s), "_", as.character(f1), "_", as.character(f2), "_", as.character(f3), "_", as.character(f4), "_", as.character(f_sd1), "_it", as.character(iter), "_de", as.character(de), ".png", sep = "") 
png(filename, height=sw/sv, width=sw, units="cm", res=300, bg = "white") 
image(ras_ridges_corr, col = gray.colors(25, start = 0.97, end = 0.5))   
points(spdf_g, pch=16, cex=0.4)  
dev.off()

filename <- paste(arbverz, "/5pictures/ridges_rn_", as.character(rw), "_", as.character(s), "_", as.character(f1), "_", as.character(f2), "_", as.character(f3), "_", as.character(f4), "_", as.character(f_sd1), "_it", as.character(iter), "_de", as.character(de), ".png", sep = "") 
png(filename, height=sw/sv, width=sw, units="cm", res=500, bg = "white") 
image(rn_grid, col = gray.colors(25, start = 0.97, end = 0.5))   
points(spdf_g, pch=16, cex=0.03)  
dev.off()

# 6.5 Analysemetadaten
#--------------------------------------------------------------------------------------------------------------------------
time2 <- Sys.time()  # end of calculation
zi <- (time2-time1)  # time of calculation

filename <- paste(arbverz, "/4result/ridges_", as.character(rw), "_", as.character(s), "_", as.character(f1), "_", as.character(f2), "_", as.character(f3), "_", as.character(f4), "_", as.character(f_sd1), "_it", as.character(iter), ".txt", sep = "")  # ridges_rw_tresh_fsd1_fsd2min_fsd2max_sup.txt

categ <- c("file      ", "rw        ", "xmin      ", "xmax      ", "ymin      ", "ymax      ", "gridcells ", "points    ", "tresh     ",   "f_sd1     ", "f1        ", "f2        ", "f3        ", "f4        ", "nn        ", "s         ", "de        ", "time      ", "time unit ")
val <- c(file_g, rw, xmin, xmax, ymin, ymax, length(ras_ridges@data$v),     ppp_g$n, tresh, f_sd1, f1, f2, f3, f4, nn, s, de, zi,as.character(attributes(zi)[1]))

write(rbind(categ,val), file = filename, ncolumns = 2, sep = "\t")

Sys.time()
zi
#===================================================================================================================================================
save.image("7ws/ws_emp_5.rws")                                                                 #load("7ws/ws_emp_5.rws")
#====================================================================================================================================================

# 6.6 Workspace
#--------------------------------------------------------------------------------------------------------------------------
# Remove all unnessesairy variables from last caluclation
rm(g)
rm(kernel.par)
rm(kernel1)
rm(kernel1d)
rm(edist)
rm(factor)
rm(gau1)
rm(crs1)
rm(crs2)
rm(crs3)
rm(date)
rm(time1)
rm(time2)
rm(rw)
rm(i)
rm(win)
rm(id)
rm(ext_ai)
rm(folder)
rm(filename)
rm(coordsnew)
rm(ddd)
rm(df)
rm(pdist)
rm(v)
rm(categ)
rm(fint)
rm(finti)
rm(fsd2)
rm(gt)
rm(ind)
rm(ind_min1)
rm(ind_min2)
rm(iter)
rm(kerpar)
rm(mwin)
rm(nn)
rm(pointges)
rm(pointnew)
rm(s)
rm(sd1)
rm(sd2)
rm(sdi)
rm(colums)
rm(stest)
rm(sv)
rm(sw)
rm(tresh)
rm(tresh_dens)
rm(val)
rm(wts)
rm(rows)
rm(zi)
rm(bb)
rm(d1)
rm(d2)
rm(de)
rm(f1)
rm(f2)
rm(f3)
rm(f4)
rm(f_sd1)
rm(ras)
rm(ras_ridges)
rm(rn_grid)
rm(rras)
rm(a)
rm(b)
rm(base_kde)
rm(dyn_kde)
rm(fs_nb_relativ)
rm(nbsldf_relativ)
rm(sgdf)
rm(sgdf_srtm)
rm(xp)
rm(file_g)
rm(file_g_ex)
rm(file_srtm)
rm(ppp_g)
rm(factor_i)
#===================================================================================================================================================
save.image("7ws/ws_emp_6.rws")  
#====================================================================================================================================================

savehistory(file="6report/empiric_pathes26112015.Rhistory")
