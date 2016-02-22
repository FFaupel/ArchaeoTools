################################################################################
## R-skript for Modelling Theoteritc Pathways in the Upper Rhine Valley
## =============================================================================
## Project: IMODEZ
##
####
######   ALL SITES, BUT NOTHING LATER THAN Lt B 
######    THEORETICAL MODELL
####
##
## Author: Franziska Faupel, with special thanks to O. Nakoinz
## Version: 01
## Date of last changes: 25.11.2015
## Data: ArkeoGIS, SHKR
## Author of data: ArkeoGIS, SHKR
## Purpose: theoretical Pathway models
## Content: 1.
## Description  -  THEORETICAL PATHWAY MODEL
##                1. Calculating junctions from empiric modell as start and endpoints in theoretical modell 
## Licence data: -
## Licence Script: GPL
## Calculation time estimated: 50h
## Intended use: IMODEZ
##
## R Version needed: R.3.2.2 (fist Part, later R 3.1.1)
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
date <- Sys.date()
date

time1 <- Sys.time()  # start time
Sys.time()

# Prepare packages
#--------------------------------------------------------------------------------------------------------------
install.packages('stats') 
install.packages('sp')
install.packages('proj4')
install.packages('rgdal')
install.packages('spatstat')
install.packages('RSQLite')
install.packages('gdata')
install.packages('KernSmooth')
install.packages('maptools')
install.packages('deldir')
install.packages('graphics')
install.packages('gdistance')
install.packages('raster')
install.packages('spdep')

#testen welche packages ich wirklich f?r den ersten Teil brauche und den Rest hier raus l?schen
library(stats)
library(sp)
#library(proj4)
library(rgdal)
library(spatstat)
#library(RSQLite)
library(gdata)
#library(KernSmooth)
library(maptools)
#library(raster)
library(deldir)
library(graphics)
library(gdistance) 



#R 3.2.2
library(ggplot2)
library(spdep)
library(plyr)
library(raster)

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

    # COST FUNCTIONS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      #tobler1993a <- function(s){6 * exp(-3.5 * abs(s + 0.05))}      # km/h Daten in m
      tobler1993b <- function(s){0.36 * exp(-3.5 * abs(s + 0.05))}   # m/min
      #herzog2012_w <- function(s){(1337.8 * s^6 + 278.19 * s^5 - 517.39 * s^4 - 78.199 * s^3 + 93.419 * s^2 + 19.825 * s + 1.64)}
      herzog2012_wi <- function(s){1/(1337.8 * s^6 + 278.19 * s^5 - 517.39 * s^4 - 78.199 * s^3 + 93.419 * s^2 + 19.825 * s + 1.64)}
      #herzog2012_d8  <- function(s){(1 + (s / 8)^2)}
      herzog2012_d8  <- function(s){(1/(1 + (s / 8)^2))}

    # Function to calculate a Cost Surface +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      hdiff <- function(x){x[2]-x[1]}

# define values
#---------------------------------------------------------------------------------------------------------------
      p <- 2500 # buffer zone around rWalk rasters, used in loop
      #iter <- 2  
      rw <- 3000 # cell size of KDE to reconstuct crossings
      sd <- 500  # standart deviation of KDE to reconstuct crossings

# Import Data
#--------------------------------------------------------------------------------------------------------------
    avs4 <- paste(arbverz,"/4result/",sep="")  
    avs3 <- paste(arbverz,"/3geodata/",sep="")  

  # empiric modell, reconstructed in seperate skript
    load("7ws/ws_emp_6.rws")

#============================================================================================================
save.image("7ws/ws_theo_1.rws")
#============================================================================================================
################################################################################
##################################################
################################
####################
#########
######
##
#          1. Calculating junctions from empiric modell as start and endpoints in theoretical modell 
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_theo_1.rws")

  # Extracting Crossings from nbl of empiric modell
  #------------------------------------------------------------------------------------------------
      jun <- data.frame(nbl_cle) # restore neighbourhood informtaion into a data frame
      i <- jun$i # seperating points from
      j <- jun$j # separating points to
      ddd <- c(i, j) # into one vector
      num_nei <- count(ddd) # counting how many connections started from one point
      cros <- num_nei[which(num_nei$freq >=3),1] # detecting all crossings
      sing <- num_nei[which(num_nei$freq == 1),1] # detecting all single ends
      coords <- cbind(rid_cle, (id=1:length(rid_cle[,1]))) #extracting coordinates
      coords <- data.frame(coords)
      coords_cros <- coords[coords$V3 %in% cros,]#extracting coordinates from crossings
      coords_sing <- coords[coords$V3 %in% sing,]#extracting coordinates from beginnings

    #removing all crossings nearer than 7,5 km
      df <- coords_cros # prepaing df to remove to close crossings
      id <- (id=1:length(coords_cros[,1]))
      df <- data.frame(df, id)
      coordinates(df)=~s1+s2
      df <- remove.duplicates(df, zero=15000, remove.second=TRUE)  #function to remove points nearer than ...m                                 
      plot(df)
      reg_cros <- df # coordinates of regional crossings

    #removing all single points nearer than 7,5 km
      df <- coords_sing # prepaing df to remove to close beginnings
      id <- (id=1:length(coords_sing[,1]))
      df <- data.frame(df, id)
      coordinates(df)=~s1+s2
      df <- remove.duplicates(df, zero=15000, remove.second=TRUE)   #function to remove points nearer than ...m                                
      plot(df)
      sing_poi <- df # coordinates of regional beginnings

    #Combining single points and crossings
      df1 <- data.frame(sing_poi)
      ddd <- (id=1:length(sing_poi[,1])) #vector needed to produce dataframe with matching number of rows
      ddd <- paste("sing_poi", ddd, sep="_") #renaming vector as preparation to do merging with
      df1 <- cbind(df1, ddd) 
      df2 <- data.frame(reg_cros)
      eee <- (id=1:length(reg_cros[,1])) #vector needed to produce dataframe with matching number of rows
      eee <- paste("reg_poi", eee, sep="_") #renaming vector as preparation to do merging with
      df2 <- cbind(df2, eee)
      reg_poi <- c(ddd,eee) # using vectors to produce dataframe with matching number of rows
      reg_poi <- data.frame(reg_poi)
      reg_poi <- merge(reg_poi, df1,  by.x="reg_poi", by.y="ddd", all.x=TRUE) #loading in regional crossings
      reg_poi <- merge(reg_poi, df2,  by.x="reg_poi", by.y="eee", all.x=TRUE) #loading in regional beginnings
      reg_poi$s1.x <- apply(reg_poi[,c("s1.x","s1.y")],1,paste,collapse="-") #combining data from merging into one colum
      reg_poi$s2.x <- apply(reg_poi[,c("s2.x","s2.y")],1,paste,collapse="-") #combining data from merging into one colum
      reg_poi$V3.x <- apply(reg_poi[,c("V3.x","V3.y")],1,paste,collapse="-") #combining data from merging into one colum
      reg_poi <- reg_poi[,c(1:4)]#deleting unnecessary colums
      reg_poi$s1.x <- gsub('NA-', '', reg_poi$s1.x) #claring colums from merging depris
      reg_poi$s2.x <- gsub('NA-', '', reg_poi$s2.x) #claring colums from merging depris
      reg_poi$s1.x <- gsub('-NA', '', reg_poi$s1.x) #claring colums from merging depris
      reg_poi$s2.x <- gsub('-NA', '', reg_poi$s2.x) #claring colums from merging depris
      reg_poi$V3.x <- gsub('-NA', '', reg_poi$V3.x) #claring colums from merging depris
      reg_poi$V3.x <- gsub('NA-', '', reg_poi$V3.x) #claring colums from merging depris
      reg_poi$s1.x <- as.numeric(reg_poi$s1.x) #need to be set to as numeric, because from merging it was as character
      reg_poi$s2.x <- as.numeric(reg_poi$s2.x) #need to be set to as numeric, because from merging it was as character

   # Calculating a nearest Neighbourhood
      coordinates(reg_poi)=~s1.x+s2.x #reg_poi as spatial Object to check weather it worked or not
      plot(reg_poi) # as controll...
      reg_poi <- remove.duplicates(reg_poi, zero=15000, remove.second=TRUE)   #function to remove points nearer than ...m  
      plot(reg_poi)
      nb_reg <- rbind(coordinates(reg_poi)) # only plain coordinates can be used 
      nb_reg <- graph2nb(relativeneigh(nb_reg)) #second neigbourhood is needed to know which of these points are likely to be connected...
      wts <- seq(1:length(nb_reg)); wts[] <- 1
      nbl_reg <- nb2lines(nb_reg, wts, reg_poi@coords, proj4string=CRS(as.character(crs1)))
      plot(nbl_reg)#controll weather is looks good

      
    # Preparation to use neigbourhood data and coordinates in loop
      reg <- data.frame(nbl_reg) #extracing data from neighbourhood lines object
      eee <- (id=1:length(reg[,1])) #vector to count/name neighbourhood partners
      reg <- data.frame(reg, eee)
      coord<- fortify(nbl_reg) # Spatial Lines Object to convert into a dataframe.
    
      coord_from <- coord[which(coord$order==1),]
      coord_from <- coord_from[,c(1,2,6)]
      coord_to <- coord[which(coord$order==2),]
      coord_to <- coord_to[,c(1,2,6)]
      reg <- merge(reg, coord_from, by.x="eee", by.y="id")
      reg <- merge(reg, coord_to, by.x="eee", by.y="id")
      
      names(reg)[names(reg)=="long.x"] <- "xFrom"
      names(reg)[names(reg)=="lat.x"] <- "yFrom"
      names(reg)[names(reg)=="long.y"] <- "xTo"
      names(reg)[names(reg)=="lat.y"] <- "yTo"
      names(reg)[names(reg)=="i"] <- "From" #renaming....
names(reg)[names(reg)=="j"] <- "To"

    # Preparing data to use in loops
      reg$xmin <- ifelse(reg$xFrom < reg$xTo, reg$xFrom, reg$xTo) #xmin and xmax is needed to trim the raster object inside the loop properly
      reg$xmax <- ifelse(reg$xFrom > reg$xTo, reg$xFrom, reg$xTo) #xmin and xmax is needed to trim the raster object inside the loop properly
      reg$ymin <- ifelse(reg$yFrom < reg$yTo, reg$yFrom, reg$yTo) #ymin and ymax is needed to trim the raster object inside the loop properly
      reg$ymax <- ifelse(reg$yFrom > reg$yTo, reg$yFrom, reg$yTo) #ymin and ymax is needed to trim the raster object inside the loop properly




# Interregional crossings
#-----------------------------------------------------------------------------------------------
# reducing crossings further....
# Calculating a nearest Neighbourhood
int <- reg_poi
coordinates(int)=~s1+s2 #reg_poi as spatial Object to check weather it worked or not
plot(int) # as controll...
nb_cros <- rbind(coordinates(int)) # only plain coordinates can be used 
nb_cros <- graph2nb(relativeneigh(nb_cros)) #second neigbourhood is needed to know which of these points are likely to be connected...
wts <- seq(1:length(nb_cros)); wts[] <- 1
nbl_cros <- nb2lines(nb_cros, wts, int@coords, proj4string=CRS(as.character(crs1)))
plot(nbl_cros)#controll weather is looks good 

cros <- data.frame(nbl_cros) # restore neighbourhood informtaion into a data frame
i <- cros$i # seperating points from
j <- cros$j # separating points to
ddd <- c(i, j) # into one vector
num_nei <- count(ddd) # counting how many connections started from one point
cros <- num_nei[which(num_nei$freq >=3),1] # detecting all crossings
end <- num_nei[which(num_nei$freq == 1),1] # detecting all single ends
coords <- data.frame(int)
eee <- c(id=1:length(int[,1]))
coords <- cbind(coords, eee) #extracting coordinates
coords_int <- coords[coords$eee%in% cros,]#extracting coordinates from crossings
coords_end <- coords[coords$eee %in% end,]#extracting coordinates from beginnings
 
#Combining single points and crossings
df1 <- data.frame(coords_end)
df1 <- df1[,c(2:4)]
ddd <- (id=1:length(coords_end[,1])) #vector needed to produce dataframe with matching number of rows
ddd <- paste("coords_end", ddd, sep="_") #renaming vector as preparation to do merging with
df1 <- cbind(df1, ddd) 
df2 <- data.frame(coords_int)
df2 <- df2[,c(2:4)]
eee <- (id=1:length(coords_int[,1])) #vector needed to produce dataframe with matching number of rows
eee <- paste("coords_int", eee, sep="_") #renaming vector as preparation to do merging with
df2 <- cbind(df2, eee)
int_poi <- c(ddd,eee) # using vectors to produce dataframe with matching number of rows
int_poi <- data.frame(int_poi)
int_poi <- merge(int_poi, df1,  by.x="int_poi", by.y="ddd", all.x=TRUE) #loading in regional crossings
int_poi <- merge(int_poi, df2,  by.x="int_poi", by.y="eee", all.x=TRUE) #loading in regional beginnings
int_poi$s1.x <- apply(int_poi[,c("s1.x.x","s1.x.y")],1,paste,collapse="-") #combining data from merging into one colum
int_poi$s2.x <- apply(int_poi[,c("s2.x.x","s2.x.y")],1,paste,collapse="-") #combining data from merging into one colum
int_poi$V3.x <- apply(int_poi[,c("V3.x.x","V3.x.y")],1,paste,collapse="-") #combining data from merging into one colum
int_poi <- int_poi[,c(1,8:10)]#deleting unnecessary colums
int_poi$s1.x <- gsub('NA-', '', int_poi$s1.x) #claring colums from merging depris
int_poi$s2.x <- gsub('NA-', '', int_poi$s2.x) #claring colums from merging depris
int_poi$s1.x <- gsub('-NA', '', int_poi$s1.x) #claring colums from merging depris
int_poi$s2.x <- gsub('-NA', '', int_poi$s2.x) #claring colums from merging depris
int_poi$V3.x <- gsub('-NA', '', int_poi$V3.x) #claring colums from merging depris
int_poi$V3.x <- gsub('NA-', '', int_poi$V3.x) #claring colums from merging depris
int_poi$s1.x <- as.numeric(int_poi$s1.x) #need to be set to as numeric, because from merging it was as character
int_poi$s2.x <- as.numeric(int_poi$s2.x) #need to be set to as numeric, because from merging it was as character
coordinates(int_poi)=~s1.x+s2.x #reg_poi as spatial Object to check weather it worked or not
plot(int_poi) # as controll...
int_poi <- remove.duplicates(int_poi, zero=15000, remove.second=TRUE)  #function to remove points nearer than ...m                                 
plot(int_poi)
# Calculating a nearest Neighbourhood
nb_int <- rbind(coordinates(int_poi)) # only plain coordinates can be used 
nb_int <- graph2nb(relativeneigh(nb_int)) #second neigbourhood is needed to know which of these points are likely to be connected...
wts <- seq(1:length(nb_int)); wts[] <- 1
nbl_int <- nb2lines(nb_int, wts, int_poi@coords, proj4string=CRS(as.character(crs1)))
plot(nbl_int)#controll weather is looks good

# Preparation to use neigbourhood data and coordinates in loop
int <- data.frame(nbl_int) #extracing data from neighbourhood lines object
eee <- (id=1:length(int[,1])) #vector to count/name neighbourhood partners
int <- data.frame(int, eee)
coord<- fortify(nbl_int) # Spatial Lines Object to convert into a dataframe.

coord_from <- coord[which(coord$order==1),]
coord_from <- coord_from[,c(1,2,6)]
coord_to <- coord[which(coord$order==2),]
coord_to <- coord_to[,c(1,2,6)]
int <- merge(int, coord_from, by.x="eee", by.y="id")
int <- merge(int, coord_to, by.x="eee", by.y="id")

names(int)[names(int)=="long.x"] <- "xFrom"
names(int)[names(int)=="lat.x"] <- "yFrom"
names(int)[names(int)=="long.y"] <- "xTo"
names(int)[names(int)=="lat.y"] <- "yTo"
names(int)[names(int)=="i"] <- "From" #renaming....
names(int)[names(int)=="j"] <- "To" 

# Preparing data to use in loops
int$xmin <- ifelse(int$xFrom < int$xTo, int$xFrom, int$xTo) #xmin and xmax is needed to trim the raster object inside the loop properly
int$xmax <- ifelse(int$xFrom > int$xTo, int$xFrom, int$xTo) #xmin and xmax is needed to trim the raster object inside the loop properly
int$ymin <- ifelse(int$yFrom < int$yTo, int$yFrom, int$yTo) #ymin and ymax is needed to trim the raster object inside the loop properly
int$ymax <- ifelse(int$yFrom > int$yTo, int$yFrom, int$yTo) #ymin and ymax is needed to trim the raster object inside the loop properly

      
#===================================================================================================================================================
save.image("7ws/ws_theo_2.rws")                                           #load("7ws/ws_theo_2.rws")
#===================================================================================================================================================
      
      # Remove all unnessesairy variables from last caluclation
      rm(df1); rm(df2); rm(coords_cros); rm(coords_sing); rm(jun); rm(num_nei); rm(cros); rm(ddd); rm(df); rm(eee); rm(i); rm(id); rm(j)
      rm(sing_poi); rm(wts); rm(coord); rm(coord_from); rm(coord_to)
      rm(coords); rm(sing);rm(reg_cros)
      rm(coords_end); rm(coords_cros); rm(nb_cros); rm(nbl_cros); 


#===================================================================================================================================================
save.image("7ws/ws_theo_3.rws") 
#===================================================================================================================================================

      ################################################################################
      ##################################################
      ################################
      ####################
      #########
      ######
      ##
      #          2. Raster!!! Using R3.1.1
      ##
      ######
      #########
      ####################
      ################################
      ##################################################
      ################################################################################
      library(stats)
      library(sp)
      library(proj4)
      library(rgdal)
      library(spatstat)
      library(RSQLite)
      library(gdata)
      library(KernSmooth)
      library(maptools)
      library(raster)
      library(deldir)
      library(graphics)
      library(gdistance) 
      library(spdep)
      library(plyr)

load("7ws/ws_theo_3.rws")

      # Preparation before Random Walk: Provide empty rasterobjekt to store random walk values
      ras_ai <- raster(sgdf_ai)
      grid_ai <- as(ras_ai, 'SpatialGridDataFrame')
      grid_ai@data$X3geodata.bw_srtm.asc <- NA
      emp_ai <- raster(grid_ai)

#===================================================================================================================================================
save.image("7ws/ws_theo_4.rws") 
#===================================================================================================================================================
################################################################################
##################################################
################################
####################
#########
######
##
#          2. HERZOG 2012 WALKING --> t1
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_theo_4.rws")

rWalk_wher1_sum <- emp_ai
rWalk_wher1_max <- emp_ai

for(i in 1:length(int[,1])){
  xmin <- int[i,11] - p 
  xmax <- int[i,12] + p
  ymin <- int[i,13] - p
  ymax <- int[i,14] + p
  ras <- crop(ras_ai, extent(xmin, xmax, ymin,ymax), filename= "corssings" , snap='near', overwrite=TRUE)
  tran_hdiff<- transition(ras, transitionFunction=hdiff, directions=8, symm=TRUE)
  slope <- geoCorrection(tran_hdiff,scl=FALSE)  
  adj <- adjacent(ras, cells=1:ncell(ras), direction=8)                     # Setting up a adjacent matrix using 8 directions to construct a transition Layer to store cost values calculated in the next step
  cost_her2012wi <- slope                                                   # storing Geo-Information in later TransitionLayer Object
  cost_her2012wi[adj] <- herzog2012_wi(slope[adj])                          # Calculating cost surface using different functions/equation
  cost_her2012wi <- geoCorrection(cost_her2012wi, scl=FALSE)                       # conductivity=cost/dist; time=1/conductivity
  from <- data.frame(int[i, c(7,8)])
  to <- data.frame(int[i, c(9,10)])
  coordinates(from)=~xFrom+yFrom; coordinates(to)=~xTo+yTo; 
  
  random <- passage(cost_her2012wi,from, to, theta = 0.0005, totalNet="total") # rWalk: expactaion for a passage of a cell
  
  rWalk_wher1_sum <- mosaic(rWalk_wher1_sum, random, fun=sum) # summing up expactation of single rWalk connections
  rWalk_wher1_max <- mosaic(rWalk_wher1_max, random, fun=max)  # using maximum value of expactation of single rWalk connections
}

plot(rWalk_wher1_sum)
plot(rWalk_wher1_max)

#===================================================================================================================================================
save.image("7ws/ws_theo_5.rws")                                           #load("7ws/ws_theo_5.rws")
#===================================================================================================================================================
################################################################################
##################################################
################################
####################
#########
######
##
#          3. HERZOG 2012 WALKING - Time is no problem --> t2
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_theo_5.rws")

rWalk_wher2_sum <- emp_ai
rWalk_wher2_max <- emp_ai

for(i in 1:length(int[,1])){
  xmin <- int[i,11] - p 
  xmax <- int[i,12] + p
  ymin <- int[i,13] - p
  ymax <- int[i,14] + p
  ras <- crop(ras_ai, extent(xmin, xmax, ymin,ymax), filename= "corssings" , snap='near', overwrite=TRUE)
  tran_hdiff<- transition(ras, transitionFunction=hdiff, directions=8, symm=TRUE)
  slope <- geoCorrection(tran_hdiff,scl=FALSE)  
  adj <- adjacent(ras, cells=1:ncell(ras), direction=8)                     # Setting up a adjacent matrix using 8 directions to construct a transition Layer to store cost values calculated in the next step
  cost_her2012wi <- slope                                                   # storing Geo-Information in later TransitionLayer Object
  cost_her2012wi[adj] <- herzog2012_wi(slope[adj])                          # Calculating cost surface using different functions/equation
  from <- data.frame(int[i,c(7,8)])
  to <- data.frame(int[i,c(9,10)])
  coordinates(from)=~xFrom+yFrom; coordinates(to)=~xTo+yTo; 
  
  random <- passage(cost_her2012wi,from, to, theta = 0.0005, totalNet="total") # rWalk: expactaion for a passage of a cell
  
  rWalk_wher2_sum <- mosaic(rWalk_wher2_sum, random, fun=sum) # summing up expactation of single rWalk connections
  rWalk_wher2_max <- mosaic(rWalk_wher2_max, random, fun=max)  # using maximum value of expactation of single rWalk connections
}

plot(rWalk_wher2_sum)
plot(rWalk_wher2_max)
#===================================================================================================================================================
save.image("7ws/ws_theo_6.rws") 
#===================================================================================================================================================
################################################################################
##################################################
################################
####################
#########
######
##
#          4. HERZOG 2012 DRIVING --> t3
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_theo_6.rws")

rWalk_dher_sum <- emp_ai
rWalk_dher_max <- emp_ai

for(i in 1:length(int[,1])){
  xmin <- int[i,11] - p 
  xmax <- int[i,12] + p
  ymin <- int[i,13] - p
  ymax <- int[i,14] + p
  ras <- crop(ras_ai, extent(xmin, xmax, ymin,ymax), filename= "corssings" , snap='near', overwrite=TRUE)
  tran_hdiff<- transition(ras, transitionFunction=hdiff, directions=8, symm=TRUE)
  slope <- geoCorrection(tran_hdiff,scl=FALSE)  
  adj <- adjacent(ras, cells=1:ncell(ras), direction=8)                     # Setting up a adjacent matrix using 8 directions to construct a transition Layer to store cost values calculated in the next step
  cost_her2012d8 <- slope                                                   # storing Geo-Information in later TransitionLayer Object
  cost_her2012d8[adj] <- herzog2012_d8(slope[adj])                          # Calculating cost surface using different functions/equation
  cost_her2012d8 <- geoCorrection(cost_her2012d8, scl=FALSE)                       # conductivity=cost/dist; time=1/conductivity
  from <- data.frame(int[i,c(7,8)])
  to <- data.frame(int[i,c(9,10)])
  coordinates(from)=~xFrom+yFrom; coordinates(to)=~xTo+yTo; 
  
  random <- passage(cost_her2012d8,from, to, theta = 0.0005, totalNet="total") # rWalk: expactaion for a passage of a cell
  
  rWalk_dher_sum <- mosaic(rWalk_dher_sum, random, fun=sum) # summing up expactation of single rWalk connections
  rWalk_dher_max <- mosaic(rWalk_dher_max, random, fun=max)  # using maximum value of expactation of single rWalk connections
  
}

plot(rWalk_dher_sum)
plot(rWalk_dher_max)
#===================================================================================================================================================
save.image("7ws/ws_theo_7.rws") 
#===================================================================================================================================================
################################################################################
##################################################
################################
####################
#########
######
##
#          5. HERZOG 2012 WALK: PRefering lower levels under 400m
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_theo_7.rws")

r <- emp_ai
r[] <-  ifelse(r@data@values <= 400, 0,r@data@values )

rWalk_wher3_sum <- emp_ai
rWalk_wher3_max <- emp_ai

################################################################################
##################################################
################################
####################
#########
######
##
#          4. Tober WALKING --> t4
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_theo_8.rws")

rWalk_wtob_sum <- emp_ai
rWalk_wtob_max <- emp_ai

for(i in 1:length(int[,1])){
  xmin <- int[i,11] - p 
  xmax <- int[i,12] + p
  ymin <- int[i,13] - p
  ymax <- int[i,14] + p
  ras <- crop(ras_ai, extent(xmin, xmax, ymin,ymax), filename= "corssings" , snap='near', overwrite=TRUE)
  tran_hdiff<- transition(ras, transitionFunction=hdiff, directions=8, symm=TRUE)
  slope <- geoCorrection(tran_hdiff,scl=FALSE)  
  adj <- adjacent(ras, cells=1:ncell(ras), direction=8)                     # Setting up a adjacent matrix using 8 directions to construct a transition Layer to store cost values calculated in the next step
  cost_tob1993b <- slope                                                   # storing Geo-Information in later TransitionLayer Object
  cost_tob1993b[adj] <- tobler1993b(slope[adj])                          # Calculating cost surface using different functions/equation
  cost_tob1993b<- geoCorrection(cost_tob1993b, scl=FALSE)                       # conductivity=cost/dist; time=1/conductivity
  from <- data.frame(int[i,c(7,8)])
  to <- data.frame(int[i,c(9,10)])
  coordinates(from)=~xFrom+yFrom; coordinates(to)=~xTo+yTo; 
  
  random <- passage(cost_tob1993b,from, to, theta = 0.0005, totalNet="total") # rWalk: expactaion for a passage of a cell
  
  rWalk_wtob_sum <- mosaic(rWalk_wtob_sum, random, fun=sum) # summing up expactation of single rWalk connections
  rWalk_wtob_max <- mosaic(rWalk_wtob_max, random, fun=max)  # using maximum value of expactation of single rWalk connections
  
}

plot(rWalk_wtob_sum)
plot(rWalk_wtob_max)
#===================================================================================================================================================
save.image("7ws/ws_theo_9.rws")                                           #load("7ws/ws_theo_9.rws")


#Hier dann reg!!!!!!


################################################################################
##################################################
################################
####################
#########
######
##
#          2. HERZOG 2012 WALKING --> t1
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_theo_4.rws")
      
rWalk_wher1_sum <- emp_ai
rWalk_wher1_max <- emp_ai

for(i in 1:length(reg[,1])){
  xmin <- reg[i,11] - p 
  xmax <- reg[i,12] + p
  ymin <- reg[i,13] - p
  ymax <- reg[i,14] + p
  ras <- crop(ras_ai, extent(xmin, xmax, ymin,ymax), filename= "corssings" , snap='near', overwrite=TRUE)
  tran_hdiff<- transition(ras, transitionFunction=hdiff, directions=8, symm=TRUE)
  slope <- geoCorrection(tran_hdiff,scl=FALSE)  
  adj <- adjacent(ras, cells=1:ncell(ras), direction=8)                     # Setting up a adjacent matrix using 8 directions to construct a transition Layer to store cost values calculated in the next step
  cost_her2012wi <- slope                                                   # storing Geo-Information in later TransitionLayer Object
  cost_her2012wi[adj] <- herzog2012_wi(slope[adj])                          # Calculating cost surface using different functions/equation
  cost_her2012wi <- geoCorrection(cost_her2012wi, scl=FALSE)                       # conductivity=cost/dist; time=1/conductivity
  from <- data.frame(reg[i, c(7,8)])
  to <- data.frame(reg[i, c(9,10)])
  coordinates(from)=~xFrom+yFrom; coordinates(to)=~xTo+yTo; 
  
  random <- passage(cost_her2012wi,from, to, theta = 0.0005, totalNet="total") # rWalk: expactaion for a passage of a cell
  
  rWalk_wher1_sum <- mosaic(rWalk_wher1_sum, random, fun=sum) # summing up expactation of single rWalk connections
  rWalk_wher1_max <- mosaic(rWalk_wher1_max, random, fun=max)  # using maximum value of expactation of single rWalk connections
}

plot(rWalk_wher1_sum)
plot(rWalk_wher1_max)

#===================================================================================================================================================
save.image("7ws/ws_theo_5.rws")                                           #load("7ws/ws_theo_5.rws")
#===================================================================================================================================================
################################################################################
##################################################
################################
####################
#########
######
##
#          3. HERZOG 2012 WALKING - Time is no problem --> t2
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_theo_5.rws")

rWalk_wher2_sum <- emp_ai
rWalk_wher2_max <- emp_ai

for(i in 1:length(reg[,1])){
  xmin <- reg[i,11] - p 
  xmax <- reg[i,12] + p
  ymin <- reg[i,13] - p
  ymax <- reg[i,14] + p
  ras <- crop(ras_ai, extent(xmin, xmax, ymin,ymax), filename= "corssings" , snap='near', overwrite=TRUE)
  tran_hdiff<- transition(ras, transitionFunction=hdiff, directions=8, symm=TRUE)
  slope <- geoCorrection(tran_hdiff,scl=FALSE)  
  adj <- adjacent(ras, cells=1:ncell(ras), direction=8)                     # Setting up a adjacent matrix using 8 directions to construct a transition Layer to store cost values calculated in the next step
  cost_her2012wi <- slope                                                   # storing Geo-Information in later TransitionLayer Object
  cost_her2012wi[adj] <- herzog2012_wi(slope[adj])                          # Calculating cost surface using different functions/equation
  from <- data.frame(reg[i,c(7,8)])
  to <- data.frame(reg[i,c(9,10)])
  coordinates(from)=~xFrom+yFrom; coordinates(to)=~xTo+yTo; 
  
  random <- passage(cost_her2012wi,from, to, theta = 0.0005, totalNet="total") # rWalk: expactaion for a passage of a cell
  
  rWalk_wher2_sum <- mosaic(rWalk_wher2_sum, random, fun=sum) # summing up expactation of single rWalk connections
  rWalk_wher2_max <- mosaic(rWalk_wher2_max, random, fun=max)  # using maximum value of expactation of single rWalk connections
}

plot(rWalk_wher2_sum)
plot(rWalk_wher2_max)
#===================================================================================================================================================
save.image("7ws/ws_theo_6.rws") 
#===================================================================================================================================================
################################################################################
##################################################
################################
####################
#########
######
##
#          4. HERZOG 2012 DRIVING --> t3
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_theo_6.rws")

rWalk_dher_sum <- emp_ai
rWalk_dher_max <- emp_ai

for(i in 1:length(reg[,1])){
  xmin <- reg[i,7] - p 
  xmax <- reg[i,8] + p
  ymin <- reg[i,9] - p
  ymax <- reg[i,10] + p
  ras <- crop(ras_ai, extent(xmin, xmax, ymin,ymax), filename= "corssings" , snap='near', overwrite=TRUE)
  tran_hdiff<- transition(ras, transitionFunction=hdiff, directions=8, symm=TRUE)
  slope <- geoCorrection(tran_hdiff,scl=FALSE)  
  adj <- adjacent(ras, cells=1:ncell(ras), direction=8)                     # Setting up a adjacent matrix using 8 directions to construct a transition Layer to store cost values calculated in the next step
  cost_her2012d8 <- slope                                                   # storing Geo-Information in later TransitionLayer Object
  cost_her2012d8[adj] <- herzog2012_d8(slope[adj])                          # Calculating cost surface using different functions/equation
  cost_her2012d8 <- geoCorrection(cost_her2012_d8, scl=FALSE)                       # conductivity=cost/dist; time=1/conductivity
  from <- data.frame(reg[i,c(3,4)])
  to <- data.frame(reg[i,c(5,6)])
  coordinates(from)=~xFrom+yFrom; coordinates(to)=~xTo+yTo; 
  
  random <- passage(cost_her2012d8,from, to, theta = 0.0005, totalNet="total") # rWalk: expactaion for a passage of a cell
  
  rWalk_dher_sum <- mosaic(rWalk_dher_sum, random, fun=sum) # summing up expactation of single rWalk connections
  rWalk_dher_max <- mosaic(rWalk_dher_max, random, fun=max)  # using maximum value of expactation of single rWalk connections
  
}

plot(rWalk_dher_sum)
plot(rWalk_dher_max)
#===================================================================================================================================================
save.image("7ws/ws_theo_7.rws") 
#===================================================================================================================================================
################################################################################
##################################################
################################
####################
#########
######
##
#          5. HERZOG 2012 WALK: PRefering lower levels under 400m
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_theo_7.rws")

r <- emp_ai
r[] <-  ifelse(r@data@values <= 400, 0,r@data@values )

rWalk_wher3_sum <- emp_ai
rWalk_wher3_max <- emp_ai

################################################################################
##################################################
################################
####################
#########
######
##
#          4. Tober WALKING --> t4
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_theo_8.rws")

rWalk_wtob_sum <- emp_ai
rWalk_wtob_max <- emp_ai

for(i in 1:length(reg[,1])){
  xmin <- reg[i,7] - p 
  xmax <- reg[i,8] + p
  ymin <- reg[i,9] - p
  ymax <- reg[i,10] + p
  ras <- crop(ras_ai, extent(xmin, xmax, ymin,ymax), filename= "corssings" , snap='near', overwrite=TRUE)
  tran_hdiff<- transition(ras, transitionFunction=hdiff, directions=8, symm=TRUE)
  slope <- geoCorrection(tran_hdiff,scl=FALSE)  
  adj <- adjacent(ras, cells=1:ncell(ras), direction=8)                     # Setting up a adjacent matrix using 8 directions to construct a transition Layer to store cost values calculated in the next step
  cost_tob1993b <- slope                                                   # storing Geo-Information in later TransitionLayer Object
  cost_tob1993b[adj] <- tobler1993b(slope[adj])                          # Calculating cost surface using different functions/equation
  cost_tob1993b<- geoCorrection(cost_tob1993b, scl=FALSE)                       # conductivity=cost/dist; time=1/conductivity
  from <- data.frame(reg[i,c(3,4)])
  to <- data.frame(reg[i,c(5,6)])
  coordinates(from)=~xFrom+yFrom; coordinates(to)=~xTo+yTo; 
  
  random <- passage(cost_tob1993b,from, to, theta = 0.0005, totalNet="total") # rWalk: expactaion for a passage of a cell
  
  rWalk_wtob_sum <- mosaic(rWalk_wtob_sum, random, fun=sum) # summing up expactation of single rWalk connections
  rWalk_wtob_max <- mosaic(rWalk_wtob_max, random, fun=max)  # using maximum value of expactation of single rWalk connections
  
}

plot(rWalk_wtob_sum)
plot(rWalk_wtob_max)
#===================================================================================================================================================
save.image("7ws/ws_theo_9.rws")                                           #load("7ws/ws_theo_9.rws")
#===================================================================================================================================================
################################################################################
##################################################
################################
####################
#########
######
##
#          5. EXPORT OF RANDOM WALK ANd RASTER DATA
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_theo_9.rws")

# 5.1 Points
#--------------------------------------------------------------------------------------------------------------------------


# 5.2 Lines
##--------------------------------------------------------------------------------------------------------------------------

# 5.3 Rasterdata
#--------------------------------------------------------------------------------------------------------------------------

# 5.4 Workspace
#--------------------------------------------------------------------------------------------------------------------------
    # Remove all unnessesairy variables from last caluclation
rm(herzog2012_d8)
rm(herzog2012_wi)
rm(kernel.par)
rm(tobler1993b)
rm(kernel1)
rm(kernel1d)
rm(edist)
rm(factor)
rm(gau1)
rm(hdiff)
rm(cost_her2012wi)
rm(cost_tob1993b)
rm(crs1)
rm(crs2)
rm(crs3)
rm(date)
rm(time1)
rm(time2)
rm(rw)
rm(p)
rm(sd)
rm(slpoe)
rm(tran_hdiff)
rm(xmax)
rm(xmin)
rm(ymax)
rm(ymin)
rm(ras)
rm(random)
rm(i)
rm(from)
rm(to)
rm(iter)
rm(slope)
rm(factor_i)
rm(adj)

rm(ras_ridges_corr)
rm(emp_ai)
rm(coord_rid)
rm(rid_cle)
rm(lsgdf_srtm)
rm(nb_cle)
rm(nbl_cle)
rm(nb_reg)
rm(sgdf_ai)


#===================================================================================================================================================
save.image("7ws/ws_theo_10.rws")                                           #load("7ws/ws_theo_10.rws")
#===================================================================================================================================================


savehistory(file="6report/random_pathe02.Rhistory")


################################################################################
##################################################
################################
####################
#########
######
##
#          OLD SKRIPT PARTS
##
######
#########
####################
################################
##################################################
################################################################################


#save.image("7ws/ws_theo_99.rws")
load("7ws/ws_theo_99.rws")



## Export Shortest Path of cost_tob1993b



## Export Shortest Path of cost_her2012d8



## Export Shortest Path of sam_poi7, SINGLE DIRECTIONS

ID <- 1
ID <- as.character(ID)
ID <- as.data.frame(ID)

filename <- paste(arbverz, "/4result/sp_single/sp_sin_21", sep = "")
sp_sin_21 <- SpatialLinesDataFrame(sp_sin_21, ID)
writeLinesShape(sp_sin_21, filename)

## Export Shortest Path of cost_her2012wi
ID <- c(1:2)
ID <- as.character(ID)
ID <- as.data.frame(ID)

filename <- paste(arbverz, "/4result/shortestPath/csp_poi7_b ", sep = "")
csp_poi7_b  <- SpatialLinesDataFrame(csp_poi7_b, ID)
writeLinesShape(csp_poi7_b, filename)

