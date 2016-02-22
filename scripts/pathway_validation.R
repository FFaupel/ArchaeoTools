################################################################################
## R-skript for Validation of Theoteritc Pathwaymodells in the Upper Rhine Valley
## =============================================================================
## Project: IMODEZ
##
####
######   ALL SITES, BUT NOTHING LATER THAN Lt B 
######    VALIDATION OF THEORETICAL MODELL
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
##                
## Licence data: -
## Licence Script: GPL
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
date <- Sys.date()
date

time1 <- Sys.time()  # start time
Sys.time()

# Prepare packages
#--------------------------------------------------------------------------------------------------------------
#install.packages('stats') 
#install.packages('sp')
#install.packages('proj4')
#install.packages('rgdal')
#install.packages('spatstat')
#install.packages('RSQLite')
#install.packages('gdata')
#install.packages('KernSmooth')
#install.packages('maptools')
#install.packages('deldir')
#install.packages('graphics')
#install.packages('gdistance')
#install.packages('raster')
#install.packages('spdep')


library(sp)
library(proj4)
library(rgdal)
library(spatstat)
#library(RSQLite)
library(gdata)
#library(KernSmooth)
library(maptools)

library(deldir)
library(graphics)
library(gdistance) 
library(spdep)
library(plyr)

#Sicher bnötigte packages
library(raster)
library(stats)
library(FNN)

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

# define values
#---------------------------------------------------------------------------------------------------------------
pd <- position_dodge(0.1)

# Import Data
#--------------------------------------------------------------------------------------------------------------
avs4 <- paste(arbverz,"/4result/",sep="")  
avs3 <- paste(arbverz,"/3geodata/",sep="")  

# theoretical modell, reconstructed in seperate skript
load("7ws/ws_theo_10.rws")

coordinates(g_ai)=~x+y

#============================================================================================================
save.image("7ws/ws_val_1.rws")
#============================================================================================================
################################################################################
##################################################
################################
####################
#########
######
##
#          1. Calculating shortest distance to grave mounds: Theoretical Modell 1
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_val_1.rws")

# Calculating the Quantile
#------------------------------------------------------------------------------------------------------------------------------------------------
t1_m <-  rWalk_wher1_max # Abbreviate strings
t1_s <-  rWalk_wher1_sum # Abbreviate strings
quan_t1_m <- quantile(t1_m, na.rm = TRUE)
quan_t1_s <- quantile(t1_s, na.rm = TRUE)
grid_t1_m <- as(t1_m, 'SpatialGridDataFrame')
grid_t1_s <- as(t1_s, 'SpatialGridDataFrame')

#Using quantile values to subdevide rWalk values
#------------------------------------------------------------------------------------------------------------------------------------------------
grid_t1_m@data$layer <- ifelse(grid_t1_m@data$layer == quan_t1_m[1], NA, 
                                    ifelse(grid_t1_m@data$layer <= quan_t1_m[2], 2000,
                                           ifelse(grid_t1_m@data$layer <= quan_t1_m[3], 3000, 
                                                  ifelse(grid_t1_m@data$layer <= quan_t1_m[4],4000, 
                                                         ifelse(grid_t1_m@data$layer > quan_t1_m[4],5000, NA)))))

grid_t1_s@data$layer <- ifelse(grid_t1_s@data$layer == quan_t1_s[1], NA, 
                               ifelse(grid_t1_s@data$layer <= quan_t1_s[2], 2000,
                                      ifelse(grid_t1_s@data$layer <= quan_t1_s[3], 3000, 
                                             ifelse(grid_t1_s@data$layer <= quan_t1_s[4],4000, 
                                                    ifelse(grid_t1_s@data$layer > quan_t1_s[4],5000, NA)))))

grid_t1_m <- raster(grid_t1_m)
grid_t1_s <- raster(grid_t1_s)
t1_m_q25 <- rbind(coordinates(grid_t1_m)[which(grid_t1_m@data@values == 2000 ),]) # extracting coordinates of first quantile
t1_m_q50 <- rbind(coordinates(grid_t1_m)[which(grid_t1_m@data@values == 3000 ),]) # extracting coordinates of second quantile
t1_m_q75 <- rbind(coordinates(grid_t1_m)[which(grid_t1_m@data@values == 4000 ),]) # extracting coordinates of third quantile
t1_m_q100 <- rbind(coordinates(grid_t1_m)[which(grid_t1_m@data@values == 5000 ),]) # extracting coordinates of forth quantile
t1_s_q25 <- rbind(coordinates(grid_t1_s)[which(grid_t1_s@data@values == 2000 ),]) # extracting coordinates of first quantile
t1_s_q50 <- rbind(coordinates(grid_t1_s)[which(grid_t1_s@data@values == 3000 ),]) # extracting coordinates of second quantile
t1_s_q75 <- rbind(coordinates(grid_t1_s)[which(grid_t1_s@data@values == 4000 ),]) # extracting coordinates of third quantile
t1_s_q100 <- rbind(coordinates(grid_t1_s)[which(grid_t1_s@data@values == 5000 ),]) # extracting coordinates of forth quantile
poi_t1_m_q25 <- data.frame(t1_m_q25)
poi_t1_s_q25 <- data.frame(t1_s_q25)
coordinates(poi_t1_m_q25)=~x+y #creating a SPatial Points object
coordinates(poi_t1_s_q25)=~x+y #creating a SPatial Points object
poi_t1_m_q50 <- data.frame(t1_m_q50)
poi_t1_s_q50 <- data.frame(t1_s_q50)
coordinates(poi_t1_m_q50)=~x+y #creating a SPatial Points object
coordinates(poi_t1_s_q50)=~x+y #creating a SPatial Points object
poi_t1_m_q75 <- data.frame(t1_m_q75)
poi_t1_s_q75 <- data.frame(t1_s_q75)
coordinates(poi_t1_m_q75)=~x+y #creating a SPatial Points object
coordinates(poi_t1_s_q75)=~x+y #creating a SPatial Points object
poi_t1_m_q100 <- data.frame(t1_m_q100)
poi_t1_s_q100 <- data.frame(t1_s_q100)
coordinates(poi_t1_m_q100)=~x+y #creating a SPatial Points object
coordinates(poi_t1_s_q100)=~x+y #creating a SPatial Points object

#Calculating (and storing) shortest distance
#-----------------------------------------------------------------------------------------------------------------------------------------------
dis_t1_m_q25 <- get.knnx(coordinates(poi_t1_m_q25), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t1_m_q50 <- get.knnx(coordinates(poi_t1_m_q50), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t1_m_q75 <- get.knnx(coordinates(poi_t1_m_q75), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t1_m_q100 <- get.knnx(coordinates(poi_t1_m_q100), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t1_s_q25 <- get.knnx(coordinates(poi_t1_s_q25), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t1_s_q50 <- get.knnx(coordinates(poi_t1_s_q50), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t1_s_q75 <- get.knnx(coordinates(poi_t1_s_q75), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t1_s_q100 <- get.knnx(coordinates(poi_t1_s_q100), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
# nn.index uses index numbers from poi_t... not from grave mounds, BUT distances are ordered, so first row, is first grave mound!!!!

# Removing unnecessariy data frames etc.
#------------------------------------------------------------------------------------------------------------------------------------------------
rm(rWalk_wher1_max); rm(rWalk_wher1_sum)
rm(grid_t1_m); rm(grid_t1_s)
rm(t1_m_q25); rm(t1_s_q25)
rm(t1_m_q50); rm(t1_s_q50)
rm(t1_m_q75); rm(t1_s_q75)
rm(t1_m_q100); rm(t1_s_q100)
rm(poi_t1_m_q25); rm(poi_t1_s_q25)
rm(poi_t1_m_q50); rm(poi_t1_s_q50)
rm(poi_t1_m_q75); rm(poi_t1_s_q75)
rm(poi_t1_m_q100); rm(poi_t1_s_q100)
#===================================================================================================================================================
save.image("7ws/ws_val_2.rws")                                           #load("7ws/ws_val_2.rws")
#===================================================================================================================================================

################################################################################
##################################################
################################
####################
#########
######
##
#          2. Calculating shortest distance to grave mounds: Theoretical Modell 2
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_val_2.rws")

# Calculating the Quantile
#------------------------------------------------------------------------------------------------------------------------------------------------
t2_m <-  rWalk_wher2_max # Abbreviate strings
t2_s <-  rWalk_wher2_sum # Abbreviate strings
quan_t2_m <- quantile(t2_m, na.rm = TRUE)
quan_t2_s <- quantile(t2_s, na.rm = TRUE)
grid_t2_m <- as(t2_m, 'SpatialGridDataFrame')
grid_t2_s <- as(t2_s, 'SpatialGridDataFrame')

#Using quantile values to subdevide rWalk values
#------------------------------------------------------------------------------------------------------------------------------------------------
grid_t2_m@data$layer <- ifelse(grid_t2_m@data$layer == quan_t2_m[1], NA, 
                               ifelse(grid_t2_m@data$layer <= quan_t2_m[2], 2000,
                                      ifelse(grid_t2_m@data$layer <= quan_t2_m[3], 3000, 
                                             ifelse(grid_t2_m@data$layer <= quan_t2_m[4],4000, 
                                                    ifelse(grid_t2_m@data$layer > quan_t2_m[4],5000, NA)))))
grid_t2_s@data$layer <- ifelse(grid_t2_s@data$layer == quan_t2_s[1], NA, 
                               ifelse(grid_t2_s@data$layer <= quan_t2_s[2], 2000,
                                      ifelse(grid_t2_s@data$layer <= quan_t2_s[3], 3000, 
                                             ifelse(grid_t2_s@data$layer <= quan_t2_s[4],4000, 
                                                    ifelse(grid_t2_s@data$layer > quan_t2_s[4],5000, NA)))))
grid_t2_m <- raster(grid_t2_m)
grid_t2_s <- raster(grid_t2_s)
t2_m_q25 <- rbind(coordinates(grid_t2_m)[which(grid_t2_m@data@values == 2000 ),]) # extracting coordinates of first quantile
t2_m_q50 <- rbind(coordinates(grid_t2_m)[which(grid_t2_m@data@values == 3000 ),]) # extracting coordinates of second quantile
t2_m_q75 <- rbind(coordinates(grid_t2_m)[which(grid_t2_m@data@values == 4000 ),]) # extracting coordinates of third quantile
t2_m_q100 <- rbind(coordinates(grid_t2_m)[which(grid_t2_m@data@values == 5000 ),]) # extracting coordinates of forth quantile
t2_s_q25 <- rbind(coordinates(grid_t2_s)[which(grid_t2_s@data@values == 2000 ),]) # extracting coordinates of first quantile
t2_s_q50 <- rbind(coordinates(grid_t2_s)[which(grid_t2_s@data@values == 3000 ),]) # extracting coordinates of second quantile
t2_s_q75 <- rbind(coordinates(grid_t2_s)[which(grid_t2_s@data@values == 4000 ),]) # extracting coordinates of third quantile
t2_s_q100 <- rbind(coordinates(grid_t2_s)[which(grid_t2_s@data@values == 5000 ),]) # extracting coordinates of forth quantile
poi_t2_m_q25 <- data.frame(t2_m_q25)
poi_t2_s_q25 <- data.frame(t2_s_q25)
coordinates(poi_t2_m_q25)=~x+y #creating a SPatial Points object
coordinates(poi_t2_s_q25)=~x+y #creating a SPatial Points object
poi_t2_m_q50 <- data.frame(t2_m_q50)
poi_t2_s_q50 <- data.frame(t2_s_q50)
coordinates(poi_t2_m_q50)=~x+y #creating a SPatial Points object
coordinates(poi_t2_s_q50)=~x+y #creating a SPatial Points object
poi_t2_m_q75 <- data.frame(t2_m_q75)
poi_t2_s_q75 <- data.frame(t2_s_q75)
coordinates(poi_t2_m_q75)=~x+y #creating a SPatial Points object
coordinates(poi_t2_s_q75)=~x+y #creating a SPatial Points object
poi_t2_m_q100 <- data.frame(t2_m_q100)
poi_t2_s_q100 <- data.frame(t2_s_q100)
coordinates(poi_t2_m_q100)=~x+y #creating a SPatial Points object
coordinates(poi_t2_s_q100)=~x+y #creating a SPatial Points object

#Calculating (and storing) shortest distance
#-----------------------------------------------------------------------------------------------------------------------------------------------
dis_t2_m_q25 <- get.knnx(coordinates(poi_t2_m_q25), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t2_m_q50 <- get.knnx(coordinates(poi_t2_m_q50), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t2_m_q75 <- get.knnx(coordinates(poi_t2_m_q75), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t2_m_q100 <- get.knnx(coordinates(poi_t2_m_q100), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t2_s_q25 <- get.knnx(coordinates(poi_t2_s_q25), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t2_s_q50 <- get.knnx(coordinates(poi_t2_s_q50), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t2_s_q75 <- get.knnx(coordinates(poi_t2_s_q75), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t2_s_q100 <- get.knnx(coordinates(poi_t2_s_q100), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
# nn.index uses index numbers from poi_t... not from grave mounds, BUT distances are ordered, so first row, is first grave mound!!!!

# Removing unnecessariy data frames etc.
#------------------------------------------------------------------------------------------------------------------------------------------------
rm(rWalk_wher2_max); rm(rWalk_wher2_sum)
rm(grid_t2_m); rm(grid_t2_s)
rm(t2_m_q25); rm(t2_s_q25)
rm(t2_m_q50); rm(t2_s_q50)
rm(t2_m_q75); rm(t2_s_q75)
rm(t2_m_q100); rm(t2_s_q100)
rm(poi_t2_m_q25); rm(poi_t2_s_q25)
rm(poi_t2_m_q50); rm(poi_t2_s_q50)
rm(poi_t2_m_q75); rm(poi_t2_s_q75)
rm(poi_t2_m_q100); rm(poi_t2_s_q100)
#===================================================================================================================================================
save.image("7ws/ws_val_3.rws")                                           #load("7ws/ws_val_3.rws")
#===================================================================================================================================================
################################################################################
##################################################
################################
####################
#########
######
##
#          3. Calculating shortest distance to grave mounds: Theoretical Modell 3
##
######
#########
####################
################################
##################################################
################################################################################

# Calculating the Quantile
#------------------------------------------------------------------------------------------------------------------------------------------------
t3_m <-  rWalk_dher_max # Abbreviate strings
t3_s <-  rWalk_dher_sum # Abbreviate strings
quan_t3_m <- quantile(t3_m, na.rm = TRUE)
quan_t3_s <- quantile(t3_s, na.rm = TRUE)
grid_t3_m <- as(t3_m, 'SpatialGridDataFrame')
grid_t3_s <- as(t3_s, 'SpatialGridDataFrame')

#Using quantile values to subdevide rWalk values
#------------------------------------------------------------------------------------------------------------------------------------------------
grid_t3_m@data$layer <- ifelse(grid_t3_m@data$layer == quan_t3_m[1], NA, 
                               ifelse(grid_t3_m@data$layer <= quan_t3_m[2], 2000,
                                      ifelse(grid_t3_m@data$layer <= quan_t3_m[3], 3000, 
                                             ifelse(grid_t3_m@data$layer <= quan_t3_m[4],4000, 
                                                    ifelse(grid_t3_m@data$layer > quan_t3_m[4],5000, NA)))))
grid_t3_s@data$layer <- ifelse(grid_t3_s@data$layer == quan_t3_s[1], NA, 
                               ifelse(grid_t3_s@data$layer <= quan_t3_s[2], 2000,
                                      ifelse(grid_t3_s@data$layer <= quan_t3_s[3], 3000, 
                                             ifelse(grid_t3_s@data$layer <= quan_t3_s[4],4000, 
                                                    ifelse(grid_t3_s@data$layer > quan_t3_s[4],5000, NA)))))
grid_t3_m <- raster(grid_t3_m)
grid_t3_s <- raster(grid_t3_s)
t3_m_q25 <- rbind(coordinates(grid_t3_m)[which(grid_t3_m@data@values == 2000 ),]) # extracting coordinates of first quantile
t3_m_q50 <- rbind(coordinates(grid_t3_m)[which(grid_t3_m@data@values == 3000 ),]) # extracting coordinates of second quantile
t3_m_q75 <- rbind(coordinates(grid_t3_m)[which(grid_t3_m@data@values == 4000 ),]) # extracting coordinates of third quantile
t3_m_q100 <- rbind(coordinates(grid_t3_m)[which(grid_t3_m@data@values == 5000 ),]) # extracting coordinates of forth quantile
t3_s_q25 <- rbind(coordinates(grid_t3_s)[which(grid_t3_s@data@values == 2000 ),]) # extracting coordinates of first quantile
t3_s_q50 <- rbind(coordinates(grid_t3_s)[which(grid_t3_s@data@values == 3000 ),]) # extracting coordinates of second quantile
t3_s_q75 <- rbind(coordinates(grid_t3_s)[which(grid_t3_s@data@values == 4000 ),]) # extracting coordinates of third quantile
t3_s_q100 <- rbind(coordinates(grid_t3_s)[which(grid_t3_s@data@values == 5000 ),]) # extracting coordinates of forth quantile
poi_t3_m_q25 <- data.frame(t3_m_q25)
poi_t3_s_q25 <- data.frame(t3_s_q25)
coordinates(poi_t3_m_q25)=~x+y #creating a SPatial Points object
coordinates(poi_t3_s_q25)=~x+y #creating a SPatial Points object
poi_t3_m_q50 <- data.frame(t3_m_q50)
poi_t3_s_q50 <- data.frame(t3_s_q50)
coordinates(poi_t3_m_q50)=~x+y #creating a SPatial Points object
coordinates(poi_t3_s_q50)=~x+y #creating a SPatial Points object
poi_t3_m_q75 <- data.frame(t3_m_q75)
poi_t3_s_q75 <- data.frame(t3_s_q75)
coordinates(poi_t3_m_q75)=~x+y #creating a SPatial Points object
coordinates(poi_t3_s_q75)=~x+y #creating a SPatial Points object
poi_t3_m_q100 <- data.frame(t3_m_q100)
poi_t3_s_q100 <- data.frame(t3_s_q100)
coordinates(poi_t3_m_q100)=~x+y #creating a SPatial Points object
coordinates(poi_t3_s_q100)=~x+y #creating a SPatial Points object

#Calculating (and storing) shortest distance
#-----------------------------------------------------------------------------------------------------------------------------------------------
dis_t3_m_q25 <- get.knnx(coordinates(poi_t3_m_q25), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t3_m_q50 <- get.knnx(coordinates(poi_t3_m_q50), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t3_m_q75 <- get.knnx(coordinates(poi_t3_m_q75), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t3_m_q100 <- get.knnx(coordinates(poi_t3_m_q100), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t3_s_q25 <- get.knnx(coordinates(poi_t3_s_q25), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t3_s_q50 <- get.knnx(coordinates(poi_t3_s_q50), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t3_s_q75 <- get.knnx(coordinates(poi_t3_s_q75), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t3_s_q100 <- get.knnx(coordinates(poi_t3_s_q100), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
# nn.index uses index numbers from poi_t... not from grave mounds, BUT distances are ordered, so first row, is first grave mound!!!!

# Removing unnecessariy data frames etc.
#------------------------------------------------------------------------------------------------------------------------------------------------
rm(rWalk_dher_max); rm(rWalk_dher_sum)
rm(grid_t3_m); rm(grid_t3_s)
rm(t3_m_q25); rm(t3_s_q25)
rm(t3_m_q50); rm(t3_s_q50)
rm(t3_m_q75); rm(t3_s_q75)
rm(t3_m_q100); rm(t3_s_q100)
rm(poi_t3_m_q25); rm(poi_t3_s_q25)
rm(poi_t3_m_q50); rm(poi_t3_s_q50)
rm(poi_t3_m_q75); rm(poi_t3_s_q75)
rm(poi_t3_m_q100); rm(poi_t3_s_q100)
#===================================================================================================================================================
save.image("7ws/ws_val_4.rws")                                           #load("7ws/ws_val_4.rws")
#===================================================================================================================================================
################################################################################
##################################################
################################
####################
#########
######
##
#          4. Calculating shortest distance to grave mounds: Theoretical Modell 4
##
######
#########
####################
################################
##################################################
################################################################################

# Calculating the Quantile
#------------------------------------------------------------------------------------------------------------------------------------------------
t4_m <-  rWalk_wtob_max # Abbreviate strings
t4_s <-  rWalk_wtob_sum # Abbreviate strings
quan_t4_m <- quantile(t4_m, na.rm = TRUE)
quan_t4_s <- quantile(t4_s, na.rm = TRUE)
grid_t4_m <- as(t4_m, 'SpatialGridDataFrame')
grid_t4_s <- as(t4_s, 'SpatialGridDataFrame')

#Using quantile values to subdevide rWalk values
#------------------------------------------------------------------------------------------------------------------------------------------------
grid_t4_m@data$layer <- ifelse(grid_t4_m@data$layer == quan_t4_m[1], NA, 
                               ifelse(grid_t4_m@data$layer <= quan_t4_m[2], 2000,
                                      ifelse(grid_t4_m@data$layer <= quan_t4_m[3], 3000, 
                                             ifelse(grid_t4_m@data$layer <= quan_t4_m[4],4000, 
                                                    ifelse(grid_t4_m@data$layer > quan_t4_m[4],5000, NA)))))
grid_t4_s@data$layer <- ifelse(grid_t4_s@data$layer == quan_t4_s[1], NA, 
                               ifelse(grid_t4_s@data$layer <= quan_t4_s[2], 2000,
                                      ifelse(grid_t4_s@data$layer <= quan_t4_s[3], 3000, 
                                             ifelse(grid_t4_s@data$layer <= quan_t4_s[4],4000, 
                                                    ifelse(grid_t4_s@data$layer > quan_t4_s[4],5000, NA)))))
grid_t4_m <- raster(grid_t4_m)
grid_t4_s <- raster(grid_t4_s)
t4_m_q25 <- rbind(coordinates(grid_t4_m)[which(grid_t4_m@data@values == 2000 ),]) # extracting coordinates of first quantile
t4_m_q50 <- rbind(coordinates(grid_t4_m)[which(grid_t4_m@data@values == 3000 ),]) # extracting coordinates of second quantile
t4_m_q75 <- rbind(coordinates(grid_t4_m)[which(grid_t4_m@data@values == 4000 ),]) # extracting coordinates of third quantile
t4_m_q100 <- rbind(coordinates(grid_t4_m)[which(grid_t4_m@data@values == 5000 ),]) # extracting coordinates of forth quantile
t4_s_q25 <- rbind(coordinates(grid_t4_s)[which(grid_t4_s@data@values == 2000 ),]) # extracting coordinates of first quantile
t4_s_q50 <- rbind(coordinates(grid_t4_s)[which(grid_t4_s@data@values == 3000 ),]) # extracting coordinates of second quantile
t4_s_q75 <- rbind(coordinates(grid_t4_s)[which(grid_t4_s@data@values == 4000 ),]) # extracting coordinates of third quantile
t4_s_q100 <- rbind(coordinates(grid_t4_s)[which(grid_t4_s@data@values == 5000 ),]) # extracting coordinates of forth quantile
poi_t4_m_q25 <- data.frame(t4_m_q25)
poi_t4_s_q25 <- data.frame(t4_s_q25)
coordinates(poi_t4_m_q25)=~x+y #creating a SPatial Points object
coordinates(poi_t4_s_q25)=~x+y #creating a SPatial Points object
poi_t4_m_q50 <- data.frame(t4_m_q50)
poi_t4_s_q50 <- data.frame(t4_s_q50)
coordinates(poi_t4_m_q50)=~x+y #creating a SPatial Points object
coordinates(poi_t4_s_q50)=~x+y #creating a SPatial Points object
poi_t4_m_q75 <- data.frame(t4_m_q75)
poi_t4_s_q75 <- data.frame(t4_s_q75)
coordinates(poi_t4_m_q75)=~x+y #creating a SPatial Points object
coordinates(poi_t4_s_q75)=~x+y #creating a SPatial Points object
poi_t4_m_q100 <- data.frame(t4_m_q100)
poi_t4_s_q100 <- data.frame(t4_s_q100)
coordinates(poi_t4_m_q100)=~x+y #creating a SPatial Points object
coordinates(poi_t4_s_q100)=~x+y #creating a SPatial Points object

#Calculating (and storing) shortest distance
#-----------------------------------------------------------------------------------------------------------------------------------------------
dis_t4_m_q25 <- get.knnx(coordinates(poi_t4_m_q25), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t4_m_q50 <- get.knnx(coordinates(poi_t4_m_q50), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t4_m_q75 <- get.knnx(coordinates(poi_t4_m_q75), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t4_m_q100 <- get.knnx(coordinates(poi_t4_m_q100), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t4_s_q25 <- get.knnx(coordinates(poi_t4_s_q25), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t4_s_q50 <- get.knnx(coordinates(poi_t4_s_q50), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t4_s_q75 <- get.knnx(coordinates(poi_t4_s_q75), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
dis_t4_s_q100 <- get.knnx(coordinates(poi_t4_s_q100), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
# nn.index uses index numbers from poi_t... not from grave mounds, BUT distances are ordered, so first row, is first grave mound!!!!

# Removing unnecessariy data frames etc.
#------------------------------------------------------------------------------------------------------------------------------------------------
rm(rWalk_wtob_max); rm(rWalk_wtob_sum)
rm(grid_t4_m); rm(grid_t4_s)
rm(t4_m_q25); rm(t4_s_q25)
rm(t4_m_q50); rm(t4_s_q50)
rm(t4_m_q75); rm(t4_s_q75)
rm(t4_m_q100); rm(t4_s_q100)
rm(poi_t4_m_q25); rm(poi_t4_s_q25)
rm(poi_t4_m_q50); rm(poi_t4_s_q50)
rm(poi_t4_m_q75); rm(poi_t4_s_q75)
rm(poi_t4_m_q100); rm(poi_t4_s_q100)
#===================================================================================================================================================
save.image("7ws/ws_val_5.rws")                                           #load("7ws/ws_val_5.rws")
#===================================================================================================================================================

################################################################################
##################################################
################################
####################
#########
######
##
#          10. 
##
######
#########
####################
################################
##################################################
################################################################################
#load("7ws/ws_val_.rws")
# Calculating distance to empiric modell
poi_emp <- data.frame(coord_rid)
coordinates(poi_emp)=~s1+s2
emp <- get.knnx(coordinates(poi_emp), coordinates(g_ai),k=1) #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)

# Combining distances of all theoretical modells

id <- (id=1:length(g_ai@coords[,1]))
# rWalk with maximum values
val_m <- data.frame(g_ai, id, dis_t1_m_q100)
names(val_m)[names(val_m)== "nn.index"] <- "t1_m_cell_q100"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_t1_m_q100"
val_m <- data.frame(val_m, dis_t1_m_q50)
names(val_m)[names(val_m)== "nn.index"] <- "t1_m_cell_q50"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_t1_m_q50"
val_m <- data.frame(val_m, dis_t1_m_q75)
names(val_m)[names(val_m)== "nn.index"] <- "t1_m_cell_q75"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_t1_m_q75"
val_m <- data.frame(val_m, dis_t1_m_q25)
names(val_m)[names(val_m)== "nn.index"] <- "t1_m_cell_q25"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_t1_m_q25"
val_m <- data.frame(val_m, dis_t2_m_q100)
names(val_m)[names(val_m)== "nn.index"] <- "t2_m_cell_q100"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_t2_m_q100"
val_m <- data.frame(val_m, dis_t2_m_q50)
names(val_m)[names(val_m)== "nn.index"] <- "t2_m_cell_q50"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_t2_m_q50"
val_m <- data.frame(val_m, dis_t2_m_q75)
names(val_m)[names(val_m)== "nn.index"] <- "t2_m_cell_q75"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_t2_m_q75"
val_m <- data.frame(val_m, dis_t2_m_q25)
names(val_m)[names(val_m)== "nn.index"] <- "t2_m_cell_q25"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_t2_m_q25"
val_m <- data.frame(val_m, dis_t3_m_q100)
names(val_m)[names(val_m)== "nn.index"] <- "t3_m_cell_q100"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_t3_m_q100"
val_m <- data.frame(val_m, dis_t3_m_q50)
names(val_m)[names(val_m)== "nn.index"] <- "t3_m_cell_q50"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_t3_m_q50"
val_m <- data.frame(val_m, dis_t3_m_q75)
names(val_m)[names(val_m)== "nn.index"] <- "t3_m_cell_q75"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_t3_m_q75"
val_m <- data.frame(val_m, dis_t3_m_q25)
names(val_m)[names(val_m)== "nn.index"] <- "t3_m_cell_q25"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_t3_m_q25"
val_m <- data.frame(val_m, dis_t4_m_q100)
names(val_m)[names(val_m)== "nn.index"] <- "t4_m_cell_q100"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_t4_m_q100"
val_m <- data.frame(val_m, dis_t4_m_q50)
names(val_m)[names(val_m)== "nn.index"] <- "t4_m_cell_q50"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_t4_m_q50"
val_m <- data.frame(val_m, dis_t4_m_q75)
names(val_m)[names(val_m)== "nn.index"] <- "t4_m_cell_q75"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_t4_m_q75"
val_m <- data.frame(val_m, dis_t4_m_q25)
names(val_m)[names(val_m)== "nn.index"] <- "t4_m_cell_q25"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_t4_m_q25"
val_m <- data.frame(val_m, emp)
names(val_m)[names(val_m)== "nn.index"] <- "cell_emp"
names(val_m)[names(val_m)== "nn.dist"] <- "dis_emp"
val_m <- rbind(val_m, c(1))
val_m$dis_t1_m_q100[2790] <- median(val_m$dis_t1_m_q100[c(1:2789)], na.rm = FALSE)
val_m$dis_t1_m_q50[2790] <- median(val_m$dis_t1_m_q50[c(1:2789)], na.rm = FALSE)
val_m$dis_t1_m_q75[2790] <- median(val_m$dis_t1_m_q75[c(1:2789)], na.rm = FALSE)
val_m$dis_t1_m_q25[2790] <- median(val_m$dis_t1_m_q25[c(1:2789)], na.rm = FALSE)
val_m$dis_t2_m_q100[2790] <- median(val_m$dis_t2_m_q100[c(1:2789)], na.rm = FALSE)
val_m$dis_t2_m_q50[2790] <- median(val_m$dis_t2_m_q50[c(1:2789)], na.rm = FALSE)
val_m$dis_t2_m_q75[2790] <- median(val_m$dis_t2_m_q75[c(1:2789)], na.rm = FALSE)
val_m$dis_t2_m_q25[2790] <- median(val_m$dis_t2_m_q25[c(1:2789)], na.rm = FALSE)
val_m$dis_t3_m_q100[2790] <- median(val_m$dis_t3_m_q100[c(1:2789)], na.rm = FALSE)
val_m$dis_t3_m_q50[2790] <- median(val_m$dis_t3_m_q50[c(1:2789)], na.rm = FALSE)
val_m$dis_t3_m_q75[2790] <- median(val_m$dis_t3_m_q75[c(1:2789)], na.rm = FALSE)
val_m$dis_t3_m_q25[2790] <- median(val_m$dis_t3_m_q25[c(1:2789)], na.rm = FALSE)
val_m$dis_t4_m_q100[2790] <- median(val_m$dis_t4_m_q100[c(1:2789)], na.rm = FALSE)
val_m$dis_t4_m_q50[2790] <- median(val_m$dis_t4_m_q50[c(1:2789)], na.rm = FALSE)
val_m$dis_t4_m_q75[2790] <- median(val_m$dis_t4_m_q75[c(1:2789)], na.rm = FALSE)
val_m$dis_t4_m_q25[2790] <- median(val_m$dis_t4_m_q25[c(1:2789)], na.rm = FALSE)
val_m$dis_emp[2790] <- median(val_m$dis_emp[c(1:2789)], na.rm = FALSE)
# preparing df to plot
val_max <- data.frame(c(1:4))
val_max$Quantile[1] <- "q25"
val_max$Quantile[2] <- "q50"
val_max$Quantile[3] <- "q75"
val_max$Quantile[4] <- "q100"
val_max$t1[4]<- val_m$dis_t1_m_q100[2790]
val_max$t1[3]<- val_m$dis_t1_m_q75[2790]
val_max$t1[2]<- val_m$dis_t1_m_q50[2790]
val_max$t1[1]<- val_m$dis_t1_m_q25[2790]
val_max$t2[4]<- val_m$dis_t2_m_q100[2790]
val_max$t2[3]<- val_m$dis_t2_m_q75[2790]
val_max$t2[2]<- val_m$dis_t2_m_q50[2790]
val_max$t2[1]<- val_m$dis_t2_m_q25[2790]
val_max$t3[4]<- val_m$dis_t3_m_q100[2790]
val_max$t3[3]<- val_m$dis_t3_m_q75[2790]
val_max$t3[2]<- val_m$dis_t3_m_q50[2790]
val_max$t3[1]<- val_m$dis_t3_m_q25[2790]
val_max$t4[4]<- val_m$dis_t4_m_q100[2790]
val_max$t4[3]<- val_m$dis_t4_m_q75[2790]
val_max$t4[2]<- val_m$dis_t4_m_q50[2790]
val_max$t4[1]<- val_m$dis_t4_m_q25[2790]

quan <- c(val_max$Quantile, val_max$Quantile, val_max$Quantile, val_max$Quantile)
quan <- sort(quan)
df_m <- data.frame(quan)
ddd <- c("t1", "t2", "t3", "t4")
df_m$mod <- paste(ddd)
df_m$dis[1] <- val_m$dis_t1_m_q100[2790]
df_m$dis[2] <- val_m$dis_t2_m_q100[2790]
df_m$dis[3] <- val_m$dis_t3_m_q100[2790]
df_m$dis[4] <- val_m$dis_t4_m_q100[2790]
df_m$dis[5] <- val_m$dis_t1_m_q25[2790]
df_m$dis[6] <- val_m$dis_t2_m_q25[2790]
df_m$dis[7] <- val_m$dis_t3_m_q25[2790]
df_m$dis[8] <- val_m$dis_t4_m_q25[2790]
df_m$dis[9] <- val_m$dis_t1_m_q50[2790]
df_m$dis[10] <- val_m$dis_t2_m_q50[2790]
df_m$dis[11] <- val_m$dis_t3_m_q50[2790]
df_m$dis[12] <- val_m$dis_t4_m_q50[2790]
df_m$dis[13] <- val_m$dis_t1_m_q75[2790]
df_m$dis[14] <- val_m$dis_t2_m_q75[2790]
df_m$dis[15] <- val_m$dis_t3_m_q75[2790]
df_m$dis[16] <- val_m$dis_t4_m_q75[2790]

filename <- paste(arbverz, "/5pictures/Vergleich_max", "_", "theoretische Modelle", ".jpeg", sep = "") 
jpeg(file=filename) 
ggplot(df_m, aes(x=mod, y=dis, group = quan, colour=quan),width=.1)+
  geom_line(position=pd)+
  geom_point(position=pd, size=3, shape=21, fill="white") +
  #geom_line(aes(colour =quan))+
  xlab("Theoretische Modelle")+
  ylab("Distanz (Median in m)")+
  ggtitle("Vergleich der theoretischen Modelle")+
  expand_limits(y=0)+
  theme_bw()+
  theme(legend.justification =c(1,0),
        legend.position=c(1,0))
dev.off()

#===================================================================================================================================================
save.image("7ws/ws_val_6.rws")                                           #load("7ws/ws_val_6.rws")
#===================================================================================================================================================



# rWalk with summed up values
val_s <- data.frame(g_ai, id, dis_t1_s_q100)
names(val_s)[names(val_s)== "nn.index"] <- "t1_s_cell_q100"
names(val_s)[names(val_s)== "nn.dist"] <- "dis_t1_s_q100"
val_s <- data.frame(val_s, dis_t1_s_q50)
names(val_s)[names(val_s)== "nn.index"] <- "t1_s_cell_q50"
names(val_s)[names(val_s)== "nn.dist"] <- "dis_t1_s_q50"
val_s <- data.frame(val_s, dis_t1_s_q75)
names(val_s)[names(val_s)== "nn.index"] <- "t1_s_cell_q75"
names(val_s)[names(val_s)== "nn.dist"] <- "dis_t1_s_q75"
val_s <- data.frame(val_s, dis_t1_s_q25)
names(val_s)[names(val_s)== "nn.index"] <- "t1_s_cell_q25"
names(val_s)[names(val_s)== "nn.dist"] <- "dis_t1_s_q25"
val_s <- data.frame(val_s, dis_t2_s_q100)
names(val_s)[names(val_s)== "nn.index"] <- "t2_s_cell_q100"
names(val_s)[names(val_s)== "nn.dist"] <- "dis_t2_s_q100"
val_s <- data.frame(val_s, dis_t2_s_q50)
names(val_s)[names(val_s)== "nn.index"] <- "t2_s_cell_q50"
names(val_s)[names(val_s)== "nn.dist"] <- "dis_t2_s_q50"
val_s <- data.frame(val_s, dis_t2_s_q75)
names(val_s)[names(val_s)== "nn.index"] <- "t2_s_cell_q75"
names(val_s)[names(val_s)== "nn.dist"] <- "dis_t2_s_q75"
val_s <- data.frame(val_s, dis_t2_s_q25)
names(val_s)[names(val_s)== "nn.index"] <- "t2_s_cell_q25"
names(val_s)[names(val_s)== "nn.dist"] <- "dis_t2_s_q25"
val_s <- data.frame(val_s, dis_t3_s_q100)
names(val_s)[names(val_s)== "nn.index"] <- "t3_s_cell_q100"
names(val_s)[names(val_s)== "nn.dist"] <- "dis_t3_s_q100"
val_s <- data.frame(val_s, dis_t3_s_q50)
names(val_s)[names(val_s)== "nn.index"] <- "t3_s_cell_q50"
names(val_s)[names(val_s)== "nn.dist"] <- "dis_t3_s_q50"
val_s <- data.frame(val_s, dis_t3_s_q75)
names(val_s)[names(val_s)== "nn.index"] <- "t3_s_cell_q75"
names(val_s)[names(val_s)== "nn.dist"] <- "dis_t3_s_q75"
val_s <- data.frame(val_s, dis_t3_s_q25)
names(val_s)[names(val_s)== "nn.index"] <- "t3_s_cell_q25"
names(val_s)[names(val_s)== "nn.dist"] <- "dis_t3_s_q25"
val_s <- data.frame(val_s, dis_t4_s_q100)
names(val_s)[names(val_s)== "nn.index"] <- "t4_s_cell_q100"
names(val_s)[names(val_s)== "nn.dist"] <- "dis_t4_s_q100"
val_s <- data.frame(val_s, dis_t4_s_q50)
names(val_s)[names(val_s)== "nn.index"] <- "t4_s_cell_q50"
names(val_s)[names(val_s)== "nn.dist"] <- "dis_t4_s_q50"
val_s <- data.frame(val_s, dis_t4_s_q75)
names(val_s)[names(val_s)== "nn.index"] <- "t4_s_cell_q75"
names(val_s)[names(val_s)== "nn.dist"] <- "dis_t4_s_q75"
val_s <- data.frame(val_s, dis_t4_s_q25)
names(val_s)[names(val_s)== "nn.index"] <- "t4_s_cell_q25"
names(val_s)[names(val_s)== "nn.dist"] <- "dis_t4_s_q25"
val_s <- rbind(val_s, c(1))
val_s$dis_t1_s_q100[2790] <- median(val_s$dis_t1_s_q100[c(1:2789)], na.rm = FALSE)
val_s$dis_t1_s_q50[2790] <- median(val_s$dis_t1_s_q50[c(1:2789)], na.rm = FALSE)
val_s$dis_t1_s_q75[2790] <- median(val_s$dis_t1_s_q75[c(1:2789)], na.rm = FALSE)
val_s$dis_t1_s_q25[2790] <- median(val_s$dis_t1_s_q25[c(1:2789)], na.rm = FALSE)
val_s$dis_t2_s_q100[2790] <- median(val_s$dis_t2_s_q100[c(1:2789)], na.rm = FALSE)
val_s$dis_t2_s_q50[2790] <- median(val_s$dis_t2_s_q50[c(1:2789)], na.rm = FALSE)
val_s$dis_t2_s_q75[2790] <- median(val_s$dis_t2_s_q75[c(1:2789)], na.rm = FALSE)
val_s$dis_t2_s_q25[2790] <- median(val_s$dis_t2_s_q25[c(1:2789)], na.rm = FALSE)
val_s$dis_t3_s_q100[2790] <- median(val_s$dis_t3_s_q100[c(1:2789)], na.rm = FALSE)
val_s$dis_t3_s_q50[2790] <- median(val_s$dis_t3_s_q50[c(1:2789)], na.rm = FALSE)
val_s$dis_t3_s_q75[2790] <- median(val_s$dis_t3_s_q75[c(1:2789)], na.rm = FALSE)
val_s$dis_t3_s_q25[2790] <- median(val_s$dis_t3_s_q25[c(1:2789)], na.rm = FALSE)
val_s$dis_t4_s_q100[2790] <- median(val_s$dis_t4_s_q100[c(1:2789)], na.rm = FALSE)
val_s$dis_t4_s_q50[2790] <- median(val_s$dis_t4_s_q50[c(1:2789)], na.rm = FALSE)
val_s$dis_t4_s_q75[2790] <- median(val_s$dis_t4_s_q75[c(1:2789)], na.rm = FALSE)
val_s$dis_t4_s_q25[2790] <- median(val_s$dis_t4_s_q25[c(1:2789)], na.rm = FALSE)
# preparing df to plot
val_sum <- data.frame(c(1:4))
val_sum$Quantile[1] <- "q25"
val_sum$Quantile[2] <- "q50"
val_sum$Quantile[3] <- "q75"
val_sum$Quantile[4] <- "q100"
val_sum$t1[4]<- val_s$dis_t1_s_q100[2790]
val_sum$t1[3]<- val_s$dis_t1_s_q75[2790]
val_sum$t1[2]<- val_s$dis_t1_s_q50[2790]
val_sum$t1[1]<- val_s$dis_t1_s_q25[2790]
val_sum$t2[4]<- val_s$dis_t2_s_q100[2790]
val_sum$t2[3]<- val_s$dis_t2_s_q75[2790]
val_sum$t2[2]<- val_s$dis_t2_s_q50[2790]
val_sum$t2[1]<- val_s$dis_t2_s_q25[2790]
val_sum$t3[4]<- val_s$dis_t3_s_q100[2790]
val_sum$t3[3]<- val_s$dis_t3_s_q75[2790]
val_sum$t3[2]<- val_s$dis_t3_s_q50[2790]
val_sum$t3[1]<- val_s$dis_t3_s_q25[2790]
val_sum$t4[4]<- val_s$dis_t4_s_q100[2790]
val_sum$t4[3]<- val_s$dis_t4_s_q75[2790]
val_sum$t4[2]<- val_s$dis_t4_s_q50[2790]
val_sum$t4[1]<- val_s$dis_t4_s_q25[2790]

quan <- c(val_sum$Quantile, val_sum$Quantile, val_sum$Quantile, val_sum$Quantile)
quan <- sort(quan)
df_s <- data.frame(quan)
ddd <- c("t1", "t2", "t3", "t4")
df_s$mod <- paste(ddd)
df_s$dis[1] <- val_s$dis_t1_s_q100[2790]
df_s$dis[2] <- val_s$dis_t2_s_q100[2790]
df_s$dis[3] <- val_s$dis_t3_s_q100[2790]
df_s$dis[4] <- val_s$dis_t4_s_q100[2790]
df_s$dis[5] <- NA
df_s$dis[6] <- NA
df_s$dis[7] <- NA
df_s$dis[8] <- NA

df_s$dis[9] <- val_s$dis_t1_s_q50[2790]
df_s$dis[10] <- val_s$dis_t2_s_q50[2790]
df_s$dis[11] <- val_s$dis_t3_s_q50[2790]
df_s$dis[12] <- val_s$dis_t4_s_q50[2790]

df_s$dis[13] <- val_s$dis_t1_s_q75[2790]
df_s$dis[14] <- val_s$dis_t2_s_q75[2790]
df_s$dis[15] <- val_s$dis_t3_s_q75[2790]
df_s$dis[16] <- val_s$dis_t4_s_q75[2790]

