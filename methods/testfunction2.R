edist  <- function(x){sqrt((x[1] - x[3])^2 + ((x[2] - x[4])^2))}  # euklidische Distanz
gau1   <- function(x, sd){dnorm(edist(x), mean=0, sd=sd)}         # euklidische Distanzen gewichtet mit normalverteilung des statischen kernels

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


