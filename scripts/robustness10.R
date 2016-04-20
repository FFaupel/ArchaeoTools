################################################################################
## R-skript: TESTING ROBUSTNESS OF RECONSTRUCTED PATHSWAY MODELLS
##        Bronze Age in Schleswig-Holstein
## =============================================================================
## 
####
######   TESTING ROBUSTNESS!!!! 10
####
##
## Author: Franziska Faupel (based on Oliver Nakoinz)
## Version: 01
## Date of last changes: 11.04.2016
## Data (Graves): Gravemounds in Schleswig-Holstein, Bronze Age
## Author of data (GRaves): GSDHL Kiel
## Data (SRTM): www.naturalearthdata.com (with contributions from F. Faupel)
## Author of data (SRTM): Tom Patterson,Nathaniel Vaughn Kelso and Franziska Faupel## Purpose: empirical pathway model, theoretical Pathway models
## Content: 1. preparation, 2. data, 3. grid, 4. static & dynamic KDE, 5. pathway detection, 
##         6. correction, 7. relative neighbour graph, 8. export
##Description - Section A: EMPIRICAL PATHWAY MODEL
##       First a static KDE, with high sd is calculated. Then a dynamic 
##       or adaptive KDE, where sd depends on the static KDE is calculated and a 
##       ridge detection applied. Dynamik KDE and ridge detection is iterativ. 
##       the grid-cells of the ridges were added to the points and the next interation is
##       calculated. Finally the relative-neighbour-graph is calculated to connect the fragments
##       of pathways. The kernel is a cosius with an added treeangular roof in the centre and 
##       1/dist in the remote areas.
## Licence data: -
## Licence Script: 
## Calculation time estimated: 20h
## Intended use: Tagung der AG Bronzezeit, TÃ¼bingen 2015 und Artikel 2016t
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


# Prepare packages
#--------------------------------------------------------------------------------------------------------------
# install.packages('sp'); install.packages('proj4'); install.packages('rgdal'); install.packages('spatstat'); install.packages('raster');install.packages('spdep');install.packages('maptools')

#sicher benÃ¶tigte packages
library(sp)
library(raster)
library(rgdal)
library(spatstat)
library(proj4)
library(spdep)
library(maptools)

date <- Sys.Date()
date

time1 <- Sys.time()  # start time
Sys.time()

# Working Directory
#---------------------------------------------------------------------------------------------------------------
arbverz <- "C:\\Diss\\Pathway modelling\\BZ_SH"
setwd(arbverz) 

# empiric modell, reconstructed in seperate skript
load("7ws/ws_rob45.rws")

file_g_ex <- "2data/rob10_2"                          # name for file with export-points for iteration
iter <- 2                                       # steps of iteration

################################################################################
##################################################
################################
####################
#########
######
##
#          1. Preparing dataset
##
######
#########
####################
################################
##################################################
################################################################################
#only 10%

g10 <- g[,c(1,2)]
g10 <- g10[sample(nrow(g10)),]
id <- (id=1:length(g[,1]))
g10 <- cbind(g10[,c(1,2)], id)
max10 <- round(661*0.1)
g10 <- g10[c(1:max10),]

spdf_g <- g10
coordinates(spdf_g)=~X+Y
spdf_g <- SpatialPoints(spdf_g)
projection(spdf_g)  <- crs4
plot(spdf_g)
#============================================================================================================
save.image("7ws/ws_rob46.rws")
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
#load("7ws/ws_rob46.rws")

### areas with low densities: smal rw kernel is advantageous 
#### areas with high densities: big rw kernel is advantageous 
### calculating a test kernel, to estimate low or high densities in areas 
#### next density estimation, modulation of kernel size to density values of test kernel  
### density values of test kernel will be used as a factor

ppp_g <- ppp(spdf_g@coords[,1], spdf_g@coords[,2], window=win)

# Variablen in AbhÃ¤Â¤ngigkeit zur mittleren nÃ¤chsten Nachbardistanz setzen
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
  
  
  ### Korrektur zu niedriger Dichtewerte und Rechenartefakte (Rauschen) die zu kÃÂÃÂ¼nstlichen Pfaden im berechneten Wegesystem fÃÂÃÂ¼hren
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
save.image("7ws/ws_rob47.rws")
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
#load("7ws/ws_rob47.rws")

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
df <- remove.duplicates(df, zero=1000, remove.second=TRUE)                                    # width of kde 100, so the radius needs to be 1000

ddd <- df@data$id                                                                             #extract id to reconstruct original coordinates
rid_cle <- coordsnew # loading in coordinates data
rid_cle <- cbind(rid_cle, id) # attaching to id (is needed to detact)
rid_cle <- rid_cle[rid_cle[,3] %in% ddd,] # deleting doubled points, which have been identified before in df
rid_cle10 <- rid_cle[,c(1:2)]                                                                 #cleaned from unnecessary id vlaues

nb_cle <- graph2nb(relativeneigh(rid_cle10)) 
wts <- seq(1:length(nb_cle)); wts[] <- 1
nbl_cle10 <- nb2lines(nb_cle, wts, rid_cle, proj4string=CRS(as.character(crs1)))
plot(nbl_cle10)

#============================================================================================================
save.image("7ws/ws_rob48.rws")
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
#load("7ws/ws_rob48.rws")

# 6.1 Points
#--------------------------------------------------------------------------------------------------------------------------
coordsnew <- rbind(coordinates(rid_cle10))
ddd <- data.frame(id=1:length(coordsnew[,1]))
pointnew = SpatialPointsDataFrame(coordsnew, ddd)
folder <- paste(arbverz, "/4result", sep = "") 
writeOGR(pointnew, folder, file_g_ex, layer="coord_den_rob10", driver="ESRI Shapefile", overwrite=TRUE)

coordsnew <- rbind(coordinates(g90))
ddd <- data.frame(id=1:length(coordsnew[,1]))
pointnew = SpatialPointsDataFrame(coordsnew, ddd)
folder <- paste(arbverz, "/4result", sep = "") 
filename <- paste(arbverz, "/g90_subsample", sep = "")
writeOGR(pointnew, folder, file_g, layer="g90", driver="ESRI Shapefile", overwrite=TRUE)

coordsnew <- rbind(coordinates(g80))
ddd <- data.frame(id=1:length(coordsnew[,1]))
pointnew = SpatialPointsDataFrame(coordsnew, ddd)
folder <- paste(arbverz, "/4result", sep = "") 
filename <- paste(arbverz, "/g80_subsample", sep = "")
writeOGR(pointnew, folder, file_g, layer="g80", driver="ESRI Shapefile", overwrite=TRUE)

coordsnew <- rbind(coordinates(g70))
ddd <- data.frame(id=1:length(coordsnew[,1]))
pointnew = SpatialPointsDataFrame(coordsnew, ddd)
folder <- paste(arbverz, "/4result", sep = "") 
filename <- paste(arbverz, "/g70_subsample", sep = "")
writeOGR(pointnew, folder, file_g, layer="g70", driver="ESRI Shapefile", overwrite=TRUE)

coordsnew <- rbind(coordinates(g60))
ddd <- data.frame(id=1:length(coordsnew[,1]))
pointnew = SpatialPointsDataFrame(coordsnew, ddd)
folder <- paste(arbverz, "/4result", sep = "") 
filename <- paste(arbverz, "/g60_subsample", sep = "")
writeOGR(pointnew, folder, file_g, layer="g60", driver="ESRI Shapefile", overwrite=TRUE)

coordsnew <- rbind(coordinates(g50))
ddd <- data.frame(id=1:length(coordsnew[,1]))
pointnew = SpatialPointsDataFrame(coordsnew, ddd)
folder <- paste(arbverz, "/4result", sep = "") 
filename <- paste(arbverz, "/g50_subsample", sep = "")
writeOGR(pointnew, folder, file_g, layer="g50", driver="ESRI Shapefile", overwrite=TRUE)

coordsnew <- rbind(coordinates(g40))
ddd <- data.frame(id=1:length(coordsnew[,1]))
pointnew = SpatialPointsDataFrame(coordsnew, ddd)
folder <- paste(arbverz, "/4result", sep = "") 
filename <- paste(arbverz, "/g40_subsample", sep = "")
writeOGR(pointnew, folder, file_g, layer="g40", driver="ESRI Shapefile", overwrite=TRUE)

coordsnew <- rbind(coordinates(g30))
ddd <- data.frame(id=1:length(coordsnew[,1]))
pointnew = SpatialPointsDataFrame(coordsnew, ddd)
folder <- paste(arbverz, "/4result", sep = "") 
filename <- paste(arbverz, "/g30_subsample", sep = "")
writeOGR(pointnew, folder, file_g, layer="g30", driver="ESRI Shapefile", overwrite=TRUE)

coordsnew <- rbind(coordinates(g20))
ddd <- data.frame(id=1:length(coordsnew[,1]))
pointnew = SpatialPointsDataFrame(coordsnew, ddd)
folder <- paste(arbverz, "/4result", sep = "") 
filename <- paste(arbverz, "/g20_subsample", sep = "")
writeOGR(pointnew, folder, file_g, layer="g20", driver="ESRI Shapefile", overwrite=TRUE)

coordsnew <- rbind(coordinates(g10))
ddd <- data.frame(id=1:length(coordsnew[,1]))
pointnew = SpatialPointsDataFrame(coordsnew, ddd)
folder <- paste(arbverz, "/4result", sep = "") 
filename <- paste(arbverz, "/g10_subsample", sep = "")
writeOGR(pointnew, folder, file_g, layer="g10", driver="ESRI Shapefile", overwrite=TRUE)

# 6.2 Lines
##--------------------------------------------------------------------------------------------------------------------------

filename <- paste(arbverz, "/4result/ridges_rn_rob10", sep = "")
writeLinesShape(nbl_cle10, filename)

# 6.3 Esri-Asc-grid
#--------------------------------------------------------------------------------------------------------------------------
filename <- paste(arbverz, "/4result/Robustness", "_", "10", ".asc", sep = "") 
writeAsciiGrid(ras_ridges_corr, filename, attr = 1, na.value = -999999, dec = ".")

# 6.5 Analysemetadaten
#--------------------------------------------------------------------------------------------------------------------------
time2 <- Sys.time()  # end of calculation
zi <- (time2-time1)  # time of calculation

filename <- paste(arbverz, "/4result/rid_rob10", as.character(rw), "_", as.character(s), "_", as.character(f1), "_", as.character(f2), "_", as.character(f3), "_", as.character(f4), "_", as.character(f_sd1), "_it", as.character(iter), ".txt", sep = "")  # ridges_rw_tresh_fsd1_fsd2min_fsd2max_sup.txt

categ <- c("file      ", "rw        ", "xmin      ", "xmax      ", "ymin      ", "ymax      ", "gridcells ", "points    ", "tresh     ",   "f_sd1     ", "f1        ", "f2        ", "f3        ", "f4        ", "nn        ", "s         ", "de        ", "time      ", "time unit ")
val <- c(file_g, rw, xmin, xmax, ymin, ymax, length(ras_ridges@data$v),     ppp_g$n, tresh, f_sd1, f1, f2, f3, f4, nn, s, de, zi,as.character(attributes(zi)[1]))

write(rbind(categ,val), file = filename, ncolumns = 2, sep = "\t")

Sys.time()
zi

#===================================================================================================================================================
save.image("7ws/ws_rob49.rws")                                                                 #load("7ws/ws_rob49.rws")
#====================================================================================================================================================

# Remove all unnessesairy variables from last caluclation
rm(coord_rid)
rm(coordsnew)
rm(ddd)
rm(pdist)
rm(rid_cle)
rm(base_kde)
rm(categ)
rm(v)
rm(a)
rm(b)
rm(d1)
rm(d2)
rm(date)
rm(df)
rm(dyn_kde)
rm(file_g_ex)
rm(fint)
rm(finti)
rm(fs_nb_relativ)
rm(fsd2)
rm(gt)
rm(i)
rm(ind)
rm(ind_min1)
rm(ind_min2)
rm(kerpar)
rm(nb_cle)
rm(nbsldf_relativ)
rm(nn)
rm(pointges)
rm(pointnew)
rm(ppp_g)
rm(ras_ridges)
rm(ras_ridges_corr)
rm(rn_grid)
rm(rras)
rm(sd1)
rm(sd2)
rm(sdi)
rm(stest)
rm(time1)
rm(time2)
rm(val)
rm(wts)
#===================================================================================================================================================
save.image("7ws/ws_rob50.rws")  
#====================================================================================================================================================

savehistory(file="6report/robustness10.Rhistory")
