################################################################################
## amod2014
## ========
## R-Skript zur PBF-Kulturraum-Analyse
## =============================================================================
## Projekt: PBF-Tagung Mainz 2014
## Skriptbearbeiter: Oliver Nakoinz
## Version: 01
## letzte Bearbeitung: 11.09.2014
## Daten: PBF: Laux (Beile, Nadeln, Lanzenspitzen)
## Datenbearbeiter: Karten digitalisiert durch O. Nakoinz
## Zweck: Vortrag
## Inhalt:  0. Vorbereitung, 1. Daten laden etc., 2. Analyseumgebung
##      vorbereiten, 3. Dichte, 4. Typenspektren, 5. Clusteranalyse, 
##      6. Distanzgraphik, 7. xxx, 8. Darstellung
## Lizenz Daten: -
## Lizenz Skript: GPL (http://www.gnu.de/documents/gpl.en.html)
################################################################################
			   ##################################################
						################################
							  ####################
								   #########
									 ######
									   ##

									   ##
									 ######
								   #########
							  ####################
						################################
			   ##################################################
################################################################################
# @@ RGEDIT LANDMARK @@: 0. preparation =======================================
# 0.1 working directory -------------------------------------------------------
arbverz <- "/home/fon/daten/analyse/pbf/pbf2014"
setwd(arbverz)

# 0.2 set variables -----------------------------------------------------------
crs1 <- "+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
crs2 <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
findesfile <- "/home/fon/daten/analyse/pbf/pbf2014/2data/pbfpkt.shp"
typesfile <- "2data/typenliste.csv"
agfile <- "/home/fon/daten/analyse/pbf/pbf2014/3geodata/ag.shp"
ufile <- "/home/fon/daten/analyse/pbf/pbf2014/3geodata/umgebung.shp"
pfile <- "/home/fon/daten/analyse/pbf/pbf2014/3geodata/profp.shp"

sdev <- 15000
rw <- 1000


# 0.3 load packages -----------------------------------------------------------
library(sp)      
library(rgdal)    
library(proj4)      
library(spatstat)
library(maptools)
library(cluster)

# 0.4 define functions --------------------------------------------------------
edist <- function(a,b){sqrt(sum((a-b) ^ 2))}
									   ##
									 ######
								   #########
							  ####################
						################################
			   ##################################################
################################################################################
# @@ RGEDIT LANDMARK @@: 1. load data =======================================
# 1.1 types -----------------------------------------------------------
types <- read.table(typesfile, header=TRUE, sep=";", dec=".") 
typvec <- as.character(types[,2])

# 1.2 finds ----------------------------------------------------------
ogrListLayers(findesfile) 
find <- readOGR(findesfile, layer="pbfpkt")

# 1.3 area of interest ------------------------------------------------
ogrListLayers(agfile) 
ag <- readOGR(agfile, layer="ag")
bbag <- bbox(ag)
win <- owin(xrange=c(bbox(ag)[1,1],bbox(ag)[1,2]), yrange=c(bbox(ag)[2,1],bbox(ag)[2,2]), unitname="m") 

ogrListLayers(ufile) 
umgebung <- readOGR(ufile, layer="umgebung")

ogrListLayers(pfile) 
profil <- readOGR(pfile, layer="profp")
proj4string(profil) = CRS(as.character(crs1))

									   ##
									 ######
								   #########
							  ####################
						################################
			   ##################################################
################################################################################
# @@ RGEDIT LANDMARK @@: 2. set up environment =================================

# 2.1 sample points ------------------------------------------------
sampp <- spsample(ag, n=100, type='regular')
proj4string(sampp) = CRS(as.character(crs1))
#sampp = SpatialPointsDataFrame(sampp, data.frame(1:length(sampp@coords[,1])))
sampp = SpatialPointsDataFrame(sampp, data.frame(cbind(1:length(sampp@coords[,1]),0)))
colnames(sampp@data) <- c("id","w")

prof <- spsample(profil, n=10, type='regular')
proj4string(prof) = CRS(as.character(crs1))
#sampp = SpatialPointsDataFrame(sampp, data.frame(1:length(sampp@coords[,1])))
prof = SpatialPointsDataFrame(prof, data.frame(cbind(1:length(prof@coords[,1]),0)))
colnames(prof@data) <- c("id","w")

# 2.2 ts ------------------------------------------------
ts1 <- 1:length(typvec)
ts1[] <- 0
ts <- ts1
for (i in 1:(length(sampp@coords[,1])-1)) { ts <-  rbind(ts,ts1) }
#ts$attributes(ts)$dimnames[[1]] <- list(1:length(sampp@coords[,1]))[[1]]

ts2 <- 1:length(typvec)
ts2[] <- 0
ts3 <- ts2
for (i in 1:(length(profil@coords[,1])-1)) { ts3 <-  rbind(ts3,ts2) }
#ts$attributes(ts)$dimnames[[1]] <- list(1:length(sampp@coords[,1]))[[1]]




							    	   ##
									 ######
								   #########
							  ####################
						################################
			   ##################################################
################################################################################
# @@ RGEDIT LANDMARK @@: 3. density =======================================
denslist <- typvec
for (i in 1:length(typvec) ) {
    t0 <- typvec[i]
    t1 <- paste("^", t0, ".*", sep = "")
    denslist[i] <- paste("dens", t0, sep = "")
    f1 <- subset(find, grepl(t1, find@data$tnr))
    pp1 <- ppp(f1@coords[,1], f1@coords[,2], window=win)  
    d1 <- density(pp1, sdev, eps=rw, edge=TRUE)
    d2 <- as.SpatialGridDataFrame.im(d1)
    proj4string(d2) = CRS(as.character(crs1))
    assign(denslist[i], d2)
}




									   ##
									 ######
								   #########
							  ####################
						################################
			   ##################################################
################################################################################
# @@ RGEDIT LANDMARK @@: 4. Typenspektren =======================================
# 4.1 sampling ------------------------------------------------
for (i in 1:length(typvec) ) {
    ts[,i] <- over(sampp, get(denslist[i]))$v
}

for (i in 1:length(typvec) ) {
    ts3[,i] <- over(profil, get(denslist[i]))$v
}


# 4.2 normalization ------------------------------------------------
# summe in probenpunkt=1
s <- apply(ts,1,sum)
for (i in 1:length(s) ) {
    ts[i,]	<- ts[i,] / s[i]
    }

s <- apply(ts3,1,sum)
for (i in 1:length(s) ) {
    ts3[i,]	<- ts3[i,] / s[i]
    }



					         		   ##
									 ######
								   #########
							  ####################
						################################
			   ##################################################
################################################################################
# @@ RGEDIT LANDMARK @@: 5. cluster analysis ====================================
# 5.1 hclust ------------------------------------------------
dis <- dist(ts, method = "euclidean")
ch1 <- hclust(dis, method="average")
plot(ch1)

# 5.2 pclust ------------------------------------------------
cp2 <- pam(ts, 2)
cp2sil <- silhouette(cp2)  
plot(cp2sil)
cp2_msil <- mean(silhouette(cp2))

cp3 <- pam(ts, 3)
cp3sil <- silhouette(cp3)  
plot(cp3sil)
cp3_msil <- mean(silhouette(cp3))

cp4 <- pam(ts, 4)
cp4sil <- silhouette(cp4)  
plot(cp4sil)
cp4_msil <- mean(silhouette(cp4))

cp5 <- pam(ts, 5)
cp5sil <- silhouette(cp5)  
plot(cp5sil)
cp5_msil <- mean(silhouette(cp5))

cp6 <- pam(ts, 6)
cp6sil <- silhouette(cp6)  
plot(cp6sil)
cp6_msil <- mean(silhouette(cp6))

cp7 <- pam(ts, 7)
cp7sil <- silhouette(cp7)  
plot(cp7sil)
cp7_msil <- mean(silhouette(cp7))

cp8 <- pam(ts, 8)
cp8sil <- silhouette(cp8)  
plot(cp8sil)
cp8_msil <- mean(silhouette(cp8))

cp9 <- pam(ts, 9)
cp9sil <- silhouette(cp9)  
plot(cp9sil)
cp9_msil <- mean(silhouette(cp9))

plot(c(cp2_msil,cp3_msil,cp4_msil,cp5_msil,cp6_msil,cp7_msil,cp8_msil,cp9_msil))

# 5.3 cluster ident ------------------------------------------------
cp2.1 <- (cp2$clustering[] == 1)
cp2.2 <- (cp2$clustering[] == 2)

cp3.1 <- (cp3$clustering[] == 1)
cp3.2 <- (cp3$clustering[] == 2)
cp3.3 <- (cp3$clustering[] == 3)

cp4.1 <- (cp4$clustering[] == 1)
cp4.2 <- (cp4$clustering[] == 2)
cp4.3 <- (cp4$clustering[] == 3)
cp4.4 <- (cp4$clustering[] == 4)

cp5.1 <- (cp5$clustering[] == 1)
cp5.2 <- (cp5$clustering[] == 2)
cp5.3 <- (cp5$clustering[] == 3)
cp5.4 <- (cp5$clustering[] == 4)
cp5.5 <- (cp5$clustering[] == 5)

cp6.1 <- (cp6$clustering[] == 1)
cp6.2 <- (cp6$clustering[] == 2)
cp6.3 <- (cp6$clustering[] == 3)
cp6.4 <- (cp6$clustering[] == 4)
cp6.5 <- (cp6$clustering[] == 5)
cp6.6 <- (cp6$clustering[] == 6)

cp7.1 <- (cp7$clustering[] == 1)
cp7.2 <- (cp7$clustering[] == 2)
cp7.3 <- (cp7$clustering[] == 3)
cp7.4 <- (cp7$clustering[] == 4)
cp7.5 <- (cp7$clustering[] == 5)
cp7.6 <- (cp7$clustering[] == 6)
cp7.7 <- (cp7$clustering[] == 7)

cp8.1 <- (cp8$clustering[] == 1)
cp8.2 <- (cp8$clustering[] == 2)
cp8.3 <- (cp8$clustering[] == 3)
cp8.4 <- (cp8$clustering[] == 4)
cp8.5 <- (cp8$clustering[] == 5)
cp8.6 <- (cp8$clustering[] == 6)
cp8.7 <- (cp8$clustering[] == 7)
cp8.8 <- (cp8$clustering[] == 8)

cp9.1 <- (cp9$clustering[] == 1)
cp9.2 <- (cp9$clustering[] == 2)
cp9.3 <- (cp9$clustering[] == 3)
cp9.4 <- (cp9$clustering[] == 4)
cp9.5 <- (cp9$clustering[] == 5)
cp9.6 <- (cp9$clustering[] == 6)
cp9.7 <- (cp9$clustering[] == 7)
cp9.8 <- (cp9$clustering[] == 8)
cp9.9 <- (cp9$clustering[] == 9)

									   ##
									 ######
								   #########
							  ####################
						################################
			   ##################################################
################################################################################
# @@ RGEDIT LANDMARK @@: 6. distanzgraphik ===================================


##
##
##
kdist <- c(
edist(ts3[1,],ts3[2,]),
edist(ts3[1,],ts3[3,]),
edist(ts3[1,],ts3[4,]),
edist(ts3[1,],ts3[5,]),
edist(ts3[1,],ts3[6,]),
edist(ts3[1,],ts3[7,]),
edist(ts3[1,],ts3[8,]),
edist(ts3[1,],ts3[9,]),
edist(ts3[1,],ts3[10,]),
edist(ts3[1,],ts3[11,]),
edist(ts3[1,],ts3[12,]),
edist(ts3[1,],ts3[13,]),
edist(ts3[1,],ts3[14,]),
edist(ts3[1,],ts3[15,]),
edist(ts3[1,],ts3[16,]),
edist(ts3[1,],ts3[17,]),
edist(ts3[1,],ts3[18,]),
edist(ts3[1,],ts3[19,]),				
edist(ts3[1,],ts3[20,]),
edist(ts3[1,],ts3[21,]),
edist(ts3[1,],ts3[22,]),
edist(ts3[1,],ts3[23,]),
edist(ts3[1,],ts3[24,])
)
						
									 ######
								   #########
							  ####################
						################################
			   ##################################################
################################################################################
# @@ RGEDIT LANDMARK @@: 7. meta ===================================
clustera1 <- cbind(
as.numeric(cp2.1),
as.numeric(cp2.2),
as.numeric(cp3.1),
as.numeric(cp3.2),
as.numeric(cp3.3),
as.numeric(cp4.1),
as.numeric(cp4.2),
as.numeric(cp4.3),
as.numeric(cp4.4),
as.numeric(cp5.1),
as.numeric(cp5.2),
as.numeric(cp5.3),
as.numeric(cp5.4),
as.numeric(cp5.5),
as.numeric(cp6.1),
as.numeric(cp6.2),
as.numeric(cp6.3),
as.numeric(cp6.4),
as.numeric(cp6.5),
as.numeric(cp6.6),
as.numeric(cp7.1),
as.numeric(cp7.2),
as.numeric(cp7.3),
as.numeric(cp7.4),
as.numeric(cp7.5),
as.numeric(cp7.6),
as.numeric(cp7.7),
as.numeric(cp8.1),
as.numeric(cp8.2),
as.numeric(cp8.3),
as.numeric(cp8.4),
as.numeric(cp8.5),
as.numeric(cp8.6),
as.numeric(cp8.7),
as.numeric(cp8.8),
as.numeric(cp9.1),
as.numeric(cp9.2),
as.numeric(cp9.3),
as.numeric(cp9.4),
as.numeric(cp9.5),
as.numeric(cp9.6),
as.numeric(cp9.7),
as.numeric(cp9.8),
as.numeric(cp9.9)
)
write.table(clustera1, file = "4result/clustera1", append = FALSE, sep = ";", dec = ".")

mclus <- kmeans(t(clustera1),3)
mclus1 <- which(mclus$cluster == 1)
mclus2 <- which(mclus$cluster == 2)
mclus3 <- which(mclus$cluster == 3)

mclus1s <- apply(clustera1[,mclus1], FUN="sum",MARGIN=1)
mclus2s <- apply(clustera1[,mclus2], FUN="sum",MARGIN=1)
mclus3s <- apply(clustera1[,mclus3], FUN="sum",MARGIN=1)

samppmc1 <- sampp
samppmc1@data <- data.frame(mclus1s) 
sampppmc1 <- as.ppp(samppmc1)
mclus1idw <- idw(sampppmc1, power=0.05, at="pixels")

samppmc2 <- sampp
samppmc2@data <- data.frame(mclus2s) 
sampppmc2 <- as.ppp(samppmc2)
mclus2idw <- idw(sampppmc2, power=0.05, at="pixels")

samppmc3 <- sampp
samppmc3@data <- data.frame(mclus3s) 
sampppmc3 <- as.ppp(samppmc3)
mclus3idw <- idw(sampppmc3, power=0.05, at="pixels")


writeAsciiGrid(as(mclus1idw, "SpatialGridDataFrame"), "4result/meta1idw.asc", attr = 1, na.value = -9999) 
writeAsciiGrid(as(mclus2idw, "SpatialGridDataFrame"), "4result/meta2idw.asc", attr = 1, na.value = -9999) 
writeAsciiGrid(as(mclus3idw, "SpatialGridDataFrame"), "4result/meta3idw.asc", attr = 1, na.value = -9999) 

									   ##
									 ######
								   #########
							  ####################
						################################
			   ##################################################
################################################################################
# @@ RGEDIT LANDMARK @@: 8. export & graphics ===================================

# 8.1 Punktkarte ------------------------------------------------
farbe <- colorRampPalette(c("white", "red"))
png("5pictures/fundkarte.png", height=12, width=15, units="cm", res=300, bg = "white")
    plot(ag, col=colors()[340], mar = c(0, 0, 0, 0),oma = c(0, 0, 0, 0))
    plot(umgebung, add=TRUE); plot(ag, col=colors()[340], add=TRUE)
    points(find, cex=0.5, pch=19, col=colors()[24])
    dev.off()

# 8.2 Dichte ------------------------------------------------
farbe <- colorRampPalette(c("yellow", "red"))
farbe <- heat.colors(n, alpha = 1)
png("5pictures/dichtekarte221b.png", height=12, width=19, units="cm", res=300, bg = "white")
    image(dens221, mar = c(0, 0, 0, 0),oma = c(0, 0, 0, 0))
    plot(ag, add=TRUE)
    points(subset(find, grepl(paste("^", "221", ".*", sep = ""), find@data$tnr)), pch=19, col=colors()[24], cex=0.5)
    dev.off()

# 8.3 hclust ------------------------------------------------
png("5pictures/dendrogram.png", height=24, width=40, units="cm", res=300, bg = "white")
    plot(ch1)
    dev.off()

# 8.4 pclust ------------------------------------------------
#641 562  657   39    31       33      142    
#rot blau grün  ocker violett  orange  gelb

png("5pictures/clusterkarten.png", height=12, width=19, units="cm", res=300, bg = "white") 
    par(mfrow=c(2,3),mar = c(0, 0, 0, 0),oma = c(0, 0, 0, 0))
        plot(ag)
            points(sampp[cp2.1,], cex=3, pch=19, col=colors()[641])
            points(sampp[cp2.2,], cex=3, pch=19, col= colors()[657])
        plot(ag)
            points(sampp[cp3.1,], cex=3 ,pch=19, col=colors()[641])
            points(sampp[cp3.2,], cex=3 ,pch=19, col= colors()[562])
            points(sampp[cp3.3,], cex=3 ,pch=19, col= colors()[657])
        plot(ag)
            points(sampp[cp4.1,], cex=3 ,pch=19, col=colors()[641])
            points(sampp[cp4.2,], cex=3 ,pch=19, col= colors()[562])
            points(sampp[cp4.3,], cex=3 ,pch=19, col= colors()[657])
            points(sampp[cp4.4,], cex=3 ,pch=19, col= colors()[39])
        plot(ag)
            points(sampp[cp5.1,], cex=3 ,pch=19, col=colors()[641])
            points(sampp[cp5.2,], cex=3 ,pch=19, col= colors()[562])
            points(sampp[cp5.3,], cex=3 ,pch=19, col= colors()[657])
            points(sampp[cp5.4,], cex=3 ,pch=19, col= colors()[31])
            points(sampp[cp5.5,], cex=3 ,pch=19, col= colors()[39])
        plot(ag)
            points(sampp[cp6.1,], cex=3 ,pch=19, col=colors()[641])
            points(sampp[cp6.2,], cex=3 ,pch=19, col= colors()[562])
            points(sampp[cp6.3,], cex=3 ,pch=19, col= colors()[33])
            points(sampp[cp6.4,], cex=3 ,pch=19, col= colors()[31])
            points(sampp[cp6.5,], cex=3 ,pch=19, col= colors()[657]) 
            points(sampp[cp6.6,], cex=3 ,pch=19, col= colors()[39])
        plot(ag)
            points(sampp[cp7.1,], cex=3 ,pch=19, col=colors()[142]) 
            points(sampp[cp7.2,], cex=3 ,pch=19, col= colors()[641])
            points(sampp[cp7.3,], cex=3 ,pch=19, col= colors()[562])
            points(sampp[cp7.4,], cex=3 ,pch=19, col= colors()[33])
            points(sampp[cp7.5,], cex=3 ,pch=19, col= colors()[31]) 
            points(sampp[cp7.6,], cex=3 ,pch=19, col= colors()[657])
            points(sampp[cp7.7,], cex=3 ,pch=19, col= colors()[39]) 
    par(mfrow=c(1,1))
    dev.off()



# 8.5 metaanalyse ------------------------------------------------
farbe <- colorRampPalette(c("white", "red"))
png("5pictures/metaclus1.png", height=12, width=19, units="cm", res=300, bg = "white")
    image(mclus1idw,col=farbe, xlab="", ylab="", xaxt="n", yaxt="n")
    plot(ag, add=TRUE)
    points(find, cex=0.1, pch=19, col=colors()[24])
    dev.off()

png("5pictures/metaclus2.png", height=12, width=19, units="cm", res=300, bg = "white")
    image(mclus2idw,col=farbe, xlab="", ylab="", xaxt="n", yaxt="n")
    plot(ag, add=TRUE)
    points(find, cex=0.1, pch=19, col=colors()[24])
    dev.off()

png("5pictures/metaclus3.png", height=12, width=19, units="cm", res=300, bg = "white")
    image(mclus3idw,col=farbe, xlab="", ylab="", xaxt="n", yaxt="n")
    plot(ag, add=TRUE)
    points(find, cex=0.1, pch=19, col=colors()[24])
    dev.off()


# 8.6 kulturdistanz ------------------------------------------------

png("5pictures/kdist.png", height=10, width=15, units="cm", res=300, bg = "white") 
    plot(1:length(kdist),kdist, col="black", pch=16, )
    lines(1:length(kdist),kdist)
    dev.off() 

png("5pictures/kdistkarte.png", height=12, width=19, units="cm", res=300, bg = "white")
    plot(ag)
    points(find, cex=0.1, pch=19, col=colors()[24])
    points(profil, cex=0.7, pch=19, col=colors()[24])
    dev.off()

