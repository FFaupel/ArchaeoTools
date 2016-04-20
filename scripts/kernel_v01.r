# Kernel-Funktionen
# Zusammengesetzt aus Kosinus und 1/x

#Soll-Parameter
xp <- 5     # x-wert von Punkt 1
s <- -0.1   # Steigung von Punkt 1

# Hilfsfunktionen

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

edist  <- function(x){sqrt((x[1] - x[3])^2 + ((x[2] - x[4])^2))}  # euklidische Distanz


#Kernelfunktionen
gau1   <- function(x, sd){dnorm(edist(x), mean=0, sd=sd)}       # euklidische Distanzen gewichtet mit normalverteilung des statischen kernels

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
    if (x <= xp) {y <- cos(x) + a}
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
    if (d <= xp) {y <- cos(x) + a}
    if (d > xp)  {y <- 1/(x + b)}
    y <- y*int
    return(y)
    }


x=pdist[22,]
x=edist(pdist[22,])
kernel1d(x=x, kp=kerpar)
kernel1d(x=pdist[22,], kp=kerpar)

#Plot
kerpar <- kernel.par(xp,s)
x <- seq(from=0, to=30000,by=1)
y <- x
for (i in 1:length(x)) {
    y[i] <- kernel1(x[i],kerpar)
    }

plot(x,y)












## gnu plot
#plot [-5:5] sin(x),cos(x),asin(x)


