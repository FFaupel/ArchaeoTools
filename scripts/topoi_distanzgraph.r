dist <- c(7,9.5,15.5,23,25.5,32.5,35)
interaction <- c(2,6,10,2,2,1,0)
cum_interaction <- c(24,21,15,5,3,1,0)

png("siz_interaktionsdist.png", height=10, width=15, units="cm", res=300, bg = "white")       
    plot(dist,cum_interaction, col="black", pch=16, xlab="Distanz (km)", ylab="Zahl ausgetauschter Scherben")
    lines(dist,interaction,lty=2,col="red")    # gestichelte Linie in Rot
    lines(dist,cum_interaction,lty=2,col="black")    # gestichelte Linie in Rot
    points(dist,cum_interaction, col="black", pch=16)
    points(dist,interaction, col="red", pch=16)          # gefuellte Punkte in Rot
    title("Interaktionsdistanzen")
    dev.off()                                                                            # png
