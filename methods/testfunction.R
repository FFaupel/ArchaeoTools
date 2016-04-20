library(sp)
library(raster) 
library(rgdal) 
library(spatstat) 
library(spdep) 
library(maptools)

source("methods/testfunction2.R")

testmatrix <- data.frame(
  x = abs(rnorm(100)*50), 
  y = abs(rnorm(100)*50)
)

plot(testmatrix$x, testmatrix$y)


xmin    <- 0
xmax    <- max(testmatrix$x)
ymin    <- 0
ymax    <- max(testmatrix$y)
ext_ai <- extent(xmin, xmax, ymin, ymax)
# raster_srtm <- raster(sgdf_srtm)
# sgdf_ai <- crop(raster_srtm, ext_ai, snap='near', overwrite=TRUE)
# sgdf_ai <- as(sgdf_ai, 'SpatialGridDataFrame')
# ras_ai <- crop(raster_srtm, ext_ai, filename="raster_ai", snap='near', overwrite=TRUE)

# Coordinates to set frame corner to define aspect ratio                      # Rahmeneckpunkte werden genutzt, das Seitenverhältnis des Ergebnisses zu definieren
sv <- (xmax-xmin)/(ymax-ymin)

rw     <- 10   # width of raster defined in m


# Defintion of frame expantion and defining frame                              # Rahmenausdehnung wird definiert und Rahmen wird erstellt
rows  <- round((ymax-ymin)/rw, 0) + 1 ; rows                                    #Anzahl der Zeilen
colums <- round((xmax-xmin)/rw, 0) + 1 ; colums                                     #Anzahl der Spalten
v <- cbind(1:(colums*rows))                                              #Vektor fÃÂ¼r alle Gridzellen
df <- data.frame(v)                                                         #Konvertiert den Vektor zu einem Dataframe
gt <- GridTopology(c(xmin, ymin), c(rw, rw), c(colums, rows))            #Gridtopology wird erstellt (c(linke untere ecke), c(ZellgrÃÂ¶ÃÂe), c(Gridausdehnung) )                       
#sgdf <- SpatialGridDataFrame(gt, df, proj4string = CRS(as.character(crs1))) # Grid wird projeziert
sgdf <- SpatialGridDataFrame(gt, df) # Grid wird projeziert





bb   <- bbox(sgdf)
win  <- owin(xrange=c(bb[1,1],bb[1,2]), yrange=c(bb[2,1],bb[2,2]), unitname="m")
# spdf_g <- remove.duplicates(spdf_g, zero=0, remove.second=TRUE)

g_ai <- data.frame(ppp_g)



# Variablen in Abhä¤ngigkeit zur mittleren nächsten Nachbardistanz setzen

spdf_g <- SpatialPoints(testmatrix)
ppp_g <- ppp(spdf_g@coords[,1], spdf_g@coords[,2], window=win)


#### sd1-Generator Funktion ####
# Generieren der Standardabweichung aus den Abständen der 
# Daten. 
## Parameter
# pppm: Point Pattern Object (spatstat)
# f_sd1: Faktor, um die Bandbreite der STD zu erhöhen 
## Output
# sd1: numeric value 
sd1gen <- function(pppm, f_sd1 = 4) {
  sd1  <- f_sd1*mean(nndist(pppm))
  return(sd1)
}

sd1 <- sd1gen(ppp_g)


#### Function to calculate static-KDE ####
# Calculting a static kde (kernel density estimation) and adding 
# the result to the input spatial grid dataframe
## Parameter
# sd1: STD defined by *sd1gen*
# sgdf: spatial grid dataframe of the area of interest
# df: input dataframe (coordinates)
# x: index of x-coordinates in df
# y: index of y-coordinates in df
## Output
# base_kde: spatial grid dataframe
makestatkde <- function(sd1, sgdf, df, x = 1, y = 2){
  stest <- 0
  if(stest != sd1) {
    base_kde <- sgdf
    for (i in seq(along=base_kde@data$v)){
      pdist <- cbind(coordinates(base_kde)[i,1],coordinates(base_kde)[i,2],df[,x],df[,y])     ### x,y,x1,y1
      base_kde@data$v[i]<- sum(apply(pdist,1,gau1,sd=sd1))
    }
  }
  stest <- sd1
  
  return(base_kde)
}

base_kde <- makestatkde(sd1 = sd1, sgdf = sgdf, df = testmatrix)
  






iter = 2

tresh  <- 0.05  # treshold, to delete phaths in areas with really low density values (KDE), meaning calculation artefacts 

f1     <- 0.2      # factor defining the minimum border of dynamic kernel (raster width) f1*mean(nn)  ## 0.2
f2     <- 0.4      # factor defining the maximum border of dynamic kernel f2*mean(nn)
f3     <- 0.5      # MinimalinentitÃÂ¤t des Kernels
f4     <- 1        # MaximalinentitÃÂ¤t des Kernels
s      <- -0.3     # Kernelparameter: incline starting from ponit 1
de     <- 0.7      # hight of additional kernell peak
sw     <- 12       # width of picture, cm
mwin   <- 9        # Mowing-window-size for ridge detection (4,9,16)
#sup    <- 1        # Faktor fÃÂ¼r die ÃÂberlagerung von statischer Dichte und dynamischer Dichte (d = sup*d_statisch + d_dynamisch) 
xp <- 15     # Kernelparameter: x-wert   von Punkt 1





for(i  in 1:iter) {       
  
  # Scaling factor for dynamic KDE 
  ## bandwidth---------------------------------------------------------------------------------------------------
  d2 <- min(base_kde@data$v)   # d1 minimal density
  d1 <- max(base_kde@data$v)   # d2 maximal density
  a  <- -(d2*f1-d1*f2)/(d1-d2) # a coefficient 1
  b  <- (f1-f2)/(d1-d2)        # b coefficient 2
  factor <- function(x){a+b*x} # factor(base_kde@data$v[i])
  fsd2 <- factor(base_kde@data$v)   # Vektor with sd values for dynamic kernel
  sd2 <- fsd2*mean(nndist(ppp_g))               # sd = faktor * mean distance to nearest neighbour
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
  
  #          3. PATHWAY DETECTION (in the loop)
  
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
  
  #          4. CORRECTION (in the loop)
 
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
} 