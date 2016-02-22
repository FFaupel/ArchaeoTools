
install.packages('spdep')
library(spdep)
install.packages('maptools')
library(maptools)
crs1 <- "+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=WGS84 +units=m +no_defs" # gk3
fs <- read.table('/home/fon/fs_ag.csv', sep=';', dec = ".", header=TRUE, stringsAsFactors=FALSE)
fs_nb_relativ <- graph2nb(relativeneigh(cbind(fs[,3],fs[,4]))) 
plot.nb(fs_nb_relativ, cbind(fs[,3],fs[,4]),  pch= 16, col="black", points=TRUE, add=FALSE, arrows=FALSE)
wts <- seq(1:length(fs_nb_relativ)); wts[] <- 1
nbsldf_relativ <- nb2lines(fs_nb_relativ, wts, cbind(fs[,3],fs[,4]), proj4string=CRS(as.character(crs1)))
writeLinesShape(nbsldf_relativ, "nbexport")




library(igraph)
library(spdep)

ids <- row.names(as(fs, "data.frame"))
fs_nb_del <- tri2nb(fs, row.names=ids)       # Delaunay-Nachbarschaftsobjekt erzeugne
m <- nb2mat(fs_nb_del)                       # Adjetanz-Matrix
#index <- which(m[,] > 0);   m[index] <- 1    # Kanten in Matrix auf 1 setzen
g <- graph.adjacency(m, mode="lower", weighted=T) # Graph erzeugen
g <- set.vertex.attribute(g, "x", index=V(g), fs[,3])   
g <- set.vertex.attribute(g, "y", index=V(g), fs[,4]) 

x1 <- V(g)$x[unlist(g[[3]][])+1] 
x2 <- V(g)$x[unlist(g[[4]][])+1]
y1 <- V(g)$y[unlist(g[[3]][])+1]
y2 <- V(g)$y[unlist(g[[4]][])+1]

g <- set.edge.attribute(g, "x1", index=E(g), x1) 
g <- set.edge.attribute(g, "x2", index=E(g), x2) 
g <- set.edge.attribute(g, "y1", index=E(g), y1) 
g <- set.edge.attribute(g, "y2", index=E(g), y2) 

plot(g)

        # Zentralitaet, Kanten
c.betweenness.kante <- edge.betweenness(g, e=E(g,directed=F))
g <- set.edge.attribute(g, "ebetween", index=E(g), c.betweenness.kante) 

# delaunaykanten mit Gewicht in .shp exportieren
LinesList <- list()         # leere LinesList Liste erstellen
for(i in seq(along= E(g))) {        # durch alle Punkte i
    m <- matrix(data = c(E(g)$x1[i],E(g)$x2[i],E(g)$y1[i],E(g)$y2[i]), nrow=2, ncol=2)    # Matrix aus zwei Koordinatenpunkten der Linienenden
    #m <- matrix(data = c(E(g)$x[i],E(g)$x[i],E(g)$y[i],E(g)$y2[i]), nrow=2, ncol=2) 
    L <- Line(m)
    LL <- list(L)
    name  <- paste("Kante", "_", i, sep="")
    LLL <- Lines(LL, ID = name)                 # Lines Objekt mit eindeutigen ID
    LinesList[length(LinesList)+1] <- LLL       # An LinesList anfuegen
    }
deldf <- data.frame(id=seq(along= E(g)),x1=E(g)$x1,y1=E(g)$y1,x2=E(g)$x2,y2=E(g)$y2, ebetween=E(g)$ebetween, weight=E(g)$weight)
sl <- SpatialLines(LinesList, proj4string = CRS(as.character(crs1)))
sdf <- SpatialLinesDataFrame(sl, deldf, match.ID = FALSE)
writeOGR(sdf, avs, "13_4_1_del_ebetween", driver="ESRI Shapefile")



