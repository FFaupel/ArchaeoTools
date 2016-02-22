################################################################################
## amod2014
## ========
## R-Skript Agentenmodell Siedlungsagglomeration
## =============================================================================
## Projekt: Vorsingen Kiel
## Skriptbearbeiter: Oliver Nakoinz
## Version: 04
## letzte Bearbeitung: 15.01.2014
## Daten: Wegerekonstruktion aus SHKR-Daten
## Datenbearbeiter: O. Nakoinz 2013
## Zweck: Agentenmodell Siedlungsagglomeration für Vorsingen Kiel 2014
## Inhalt:  0. Vorbereitung, 1. Daten laden etc., 2. Umgebung vorbereiten,
##          3. Schleife, 4. Ausgabe
## Lizenz Daten: -
## Lizenz Skript: GPL (http://www.gnu.de/documents/gpl.en.html)
################################################################################
               ##################################################
                        ################################
                              ####################
                                   #########
                                     ######
                                       ##
# @@ RGEDIT LANDMARK @@: 0. Vorbereitung =======================================
# 0.1 Variablen setzen ---------------------------------------------------------
arbverz <- "/home/fon/daten/analyse/shkr2014/amodell/amod1"
setwd(arbverz)
crs1 <- "+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=WGS84 +units=m +no_defs" # gk3
#crs2 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"  # crs-Text fuer geogr. mit WGS84

rndp     <- 30        # Anzahl der Zufallspunkte, die jedem Agenten offeriert werden.
steps    <- 50        # Schleifendurchläufe je run
tf       <- 0.5       # Handelsfaktor, der Anteil, des Reichtums, den ein Händler verhandelt
sd_adens <- 600       # sd für Bewertung der Agentendichte (für Flächenbedarf)
sd_dis_indig  <- 20000   # sd für Distanzgewichtung
sd_dis_merch  <- 30000   # sd für Distanzgewichtung
sd1      <- 6000      # sd für Punktdichte KDE
rb       <- 5         # Randbreite für "gehe nach Süden"
sf       <- 2250000   # Skalierungsfaktor für kde (punkte je Zelle statt Punkte je qm)

nag_merch   <-  20  # Anzahl Händer
nag_indig   <- 100  # Anzahl Einheimische

dist_exp_indig <- 1
road_exp_indig <- 0.3
weal_exp_indig <- 0.5
area_exp_indig <- 1
nort_exp_indig <- 0

dist_exp_merch <- 0.3
road_exp_merch <- 1
weal_exp_merch <- 0.7
area_exp_merch <- 0
nort_exp_merch <- 1

# Später einbauen:
#   Krisenanfälligkeit  
#   Lernfähigkeit 
#   Sättigungseffekte
#   Distanzen reduziern Güter/Reichtum
#   Produktion und Verbrauch von Agrarpodukten und Handelsgütern

# 0.2 Pakete -------------------------------------------------------------------
library(sqldf)
library(lattice)
library(rgdal) 
library(sp)
library(spatstat)
library(raster)
library(maptools)

# 0.3 Funktionen ---------------------------------------------------------------
edist  <- function(x){sqrt((x[1] - x[3])^2 + ((x[2] - x[4])^2))}  # euklidische Distanz (x1,y1,x2,y2
gau1   <- function(x, sd){dnorm(edist(x), mean=0, sd=sd)}       # euklidische Distanzen gewichtet mit normalverteilung des statischen kernels
kde <- function(pp,grid,sd1,f=1){    # pp: zwei spalten cbind(x,y); grid: spgdf; sd1: stdev für kde; f=skalierungsfaktor
    grid2 <- grid
    for (i in seq(along=grid2@data[,1])){       
        pdist <- cbind(coordinates(grid2)[i,1],coordinates(grid2)[i,2],pp[,1],pp[,2])### x,y,x1,y1
        grid2@data[i,1]<- sum(apply(pdist,1,gau1,sd=sd1)) * f
        }
    return(grid2)
    }

# kde für ausgewählte Punkte im Raster, deren index im vektor smp steht
kdep <- function(pp,grid,sd1,smp,f=1){    # pp: zwei spalten cbind(x,y); grid: spgdf; sd1: stdev für kde; smp: index der Probenpunkte
    ret <- seq(1:length(smp))
    j <- 1
    for (i in smp){       
        pdist <- cbind(coordinates(grid)[i,1],coordinates(grid)[i,2],pp[,1],pp[,2]) ### x,y,x1,y1
        ret[j]<- sum(apply(pdist,1,gau1,sd=sd1)) * f
        j <- j+1
        }
    return(ret)
    }
                                       ##
                                     ######
                                   #########
                              ####################
                        ################################
               ##################################################
################################################################################
# @@ RGEDIT LANDMARK @@: 1. Daten laden etc. ===================================
#wegedistc <- readGDAL("wegedist.tif")  
wegedistc <- readGDAL("wegedist3000.tif") 
zz        <- length(wegedistc@data[,1])   # Zellenanzahl des Rasters
spalten   <- wegedistc@grid@cells.dim[1]
zeilen    <- wegedistc@grid@cells.dim[2]

                                       ##
                                     ######
                                   #########
                              ####################
                        ################################
               ##################################################
################################################################################
# @@ RGEDIT LANDMARK @@: 2. Umgebung einrichten ================================
# Spielbrett wegedistc dient als Vorlage
#===========
no <- wegedistc    # Nordraster
    w <- rep(zeilen,spalten)
    for (i in zeilen:2){
        w <- c(w, rep(i,spalten))
        }
    no@data[,1] <- w/zeilen
    
# Agenten (id,atyp,rzelle, reichtum) 
#========
# Einheimischer, typ 1
id          <- 0
rzelle      <- 0
atyp        <- 1  
reichtum    <- 0.05
agenttemplate1 <- data.frame(cbind(id=rep(id,nag_indig),atyp=rep(atyp,nag_indig),rzelle=rep(rzelle,nag_indig),reichtum=rep(reichtum,nag_indig)))
agenttemplate1[,3] <- floor(runif(nag_indig,min=0,max=zz))
# Händler, typ 2
atyp        <- 2 
reichtum    <- 1
agenttemplate2 <- data.frame(cbind(id=rep(id,nag_merch),atyp=rep(atyp,nag_merch),rzelle=rep(rzelle,nag_merch),reichtum=rep(reichtum,nag_merch)))
agenttemplate2[,3] <- floor(runif(nag_merch,min=spalten*zeilen-spalten,max=spalten*zeilen))
a <- rbind(agenttemplate1,agenttemplate2)

# Transformtationsfunktionen wandelt Parameter in Bewertungskriterien um
#===========================
# Distanzgewichtung 
eval.dist.indig  <- function(dist){(sd_dis_indig*dnorm(dist, mean=0, sd=sd_dis_indig))}
eval.dist.merch  <- function(dist){(sd_dis_merch*dnorm(dist, mean=0, sd=sd_dis_merch))}

# Wegerdistanzraster dist aus wegedistc  0 km = 1; 63 km = 0; linear
x  <- c(63000,0)
y  <- c(0,1)
fm <- lm(y ~ x)
a.wdist <- fm$coefficients[1]
b.wdist <- fm$coefficients[2]
eval.wdist  <- function(dist){a.wdist + b.wdist*dist}  

# Reichtumsgewichtung, linear
x <- c(0,30)
y <- c(0.1,1)
fm <- lm(y ~ x)
a.wlt <- fm$coefficients[1]
b.wlt <- fm$coefficients[2]
eval.wlt  <- function(dist){a.wlt + b.wlt*dist}  

# Agentendichtegewichtung: Flächenbedarfgrenze (normalverteilt)
eval.agdens.gau  <- function(dist){sd_adens*dnorm(dist, mean=0, sd=sd_adens)}

# Karten der Ausgangslage
# =======================
ergebnis_0  <- kde(coordinates(wegedistc)[a[1:nag_indig,3],], wegedistc, sd=3000,f=sf)
ergebnis_t0 <- kde(coordinates(wegedistc)[a[(nag_indig+1):(nag_indig+nag_merch),3],], wegedistc, sd=3000,f=sf)
                                       ##
                                     ######
                                   #########
                              ####################
                        ################################
               ##################################################
################################################################################
# @@ RGEDIT LANDMARK @@: 3. Schleife ===========================================
for (i in 1:steps){ 
    for (ai in seq(along=a[,3])){   # durch alle Agenten
        air <- a[ai,3] # Agentenrasterindex
            # Zielpunkte und deren Eigenschaften
            zpi     <- floor(runif(rndp,1,zz))  # RZellnummern von Zufallspunkten
            zpi.x   <- coordinates(wegedistc)[zpi,1]
            zpi.y   <- coordinates(wegedistc)[zpi,2]
            zpi.wd  <- wegedistc@data[zpi,1]      # Distanz zum Weg
            zpi.d   <-   apply(cbind(rep(coordinates(wegedistc)[air,1],rndp),rep(coordinates(wegedistc)[air,2],rndp),zpi.x,zpi.y),1,edist)       # Distanz zu den neuen Punkten 
            zpi.kde <- kdep(coordinates(wegedistc)[a[1:nag_indig,3],], wegedistc, sd=sd1,smp=zpi,f=sf)  # Dichte für Zufallspunkte  
            zpi.n   <- no@data[zpi,1]     # Nordwert für Zufallspunkte
            zpi.wlt <- rep(0,length(zpi))  # Reichtum
                iii <- 1
            index <- c(0)
            for (zzi in zpi){
                index <- intersect(which(a[,2] == 1), which(a[,3] == zzi))
                zpi.wlt[iii] <- sum(a[index,4])
                iii <- iii + 1
                }
            #Standorteigenschaften
            stp.wd <- wegedistc@data[air,1]
            stp.d <- 0
            stp.kde <- kdep(coordinates(wegedistc)[a[1:nag_indig,3],], wegedistc, sd=sd1,smp=c(air,air),f=sf)[1]
            stp.n   <- no@data[air,1]
            stp.wlt <- sum(a[intersect(which(a[,2] == 1), which(a[,3] == air)),4])

            if (a[ai,2]==1){  # für Einheimische
                # Bewertung
                val <- ((eval.dist.indig(zpi.d)  )^dist_exp_indig * 
                        (eval.wdist(zpi.wd)      )^road_exp_indig *
                        (eval.wlt(zpi.wlt)       )^weal_exp_indig *
                        (eval.agdens.gau(zpi.kde))^area_exp_indig)

                # Neues Ziel ermitteln
                ziel <- zpi[which(val==max(val))]
                ziel <- ziel[floor(runif(1,min=1,max=length(ziel)+0.9999999))]
              
                # Bleibe wenn Standpunkt besser

                valdatstp <- ((1                  )^dist_exp_indig * 
                        (eval.wdist(stp.wd)      )^road_exp_indig *
                        (eval.wlt(stp.wlt)       )^weal_exp_indig *
                        (eval.agdens.gau(stp.kde))^area_exp_indig)
                if (valdatstp > max(val)) ziel <- air

                # agiere
                #     schreite:
                a[ai,3] <- ziel
                }

            if (a[ai,2]==2){  # für Händler
                val <- ((eval.dist.indig(zpi.d)  )^dist_exp_merch * 
                        (eval.wdist(zpi.wd)      )^road_exp_merch *
                        (eval.wlt(zpi.wlt)       )^weal_exp_merch *
                        (eval.agdens.gau(zpi.kde))^area_exp_merch *
                                             zpi.n^nort_exp_merch)

                # Neues Ziel ermitteln
                ziel <- zpi[which(val==max(val))]
                ziel <- ziel[floor(runif(1,min=1,max=length(ziel)+0.9999999))]

                # agiere
                #     schreite:
                a[ai,3] <- ziel
                index <- c(0)
                #     handele
                index <- intersect(which(a[,2] == 1), which(a[,3] == a[ai,3]))
                if (length(index) > 0){
                        r <- a[index,4] + 1
                        verteilung <- tf*a[ai,4]
                        a[ai,4] <- a[ai,4]-verteilung
                        anteil <- r*verteilung/sum(r)
                        a[index,4] <- a[index,4] + anteil
                        }

                #      gehe zurück nach Süden wenn du am Nordrand bist 
                if (ziel < spalten*rb){
                    a[ai,3] <- floor(runif(1,min=spalten*zeilen-(rb*spalten),max=spalten*zeilen))
                    a[ai,4] <- 1
                    }
               }
       }
}
                                       ##
                                     ######
                                   #########
                              ####################
                        ################################
               ##################################################
################################################################################
# @@ RGEDIT LANDMARK @@: 4. Ausgabe ============================================
ergebnis_1  <- kde(coordinates(wegedistc)[a[1:nag_indig,3],], wegedistc, sd=3000,f=sf)
ergebnis_t1 <- kde(coordinates(wegedistc)[a[(nag_indig+1):(nag_indig+nag_merch),3],], wegedistc, sd=3000,f=sf)

par(mfrow=c(2,2))
image(ergebnis_0)
image(ergebnis_1)
image(ergebnis_t0)
image(ergebnis_t1)
par(mfrow=c(1,1))

################################################################################
               ##################################################
                        ################################
                              ####################
                                   #########
                                     #ENDE#
                                       ##


