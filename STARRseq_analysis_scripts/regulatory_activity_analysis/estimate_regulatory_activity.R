setwd("~/Desktop/sapelo2/atac-starr-seq_maize/v5/regulatory_fragment_activity")

# load libraries
library(MASS)
library(viridis)
library(scales)
library(densityClust)
library(RColorBrewer)

### functions ###

# scale 0-1
normalize <- function(x){
    (x - min(x))/(max(x)-min(x))
}

# plot density scatter
plot_colorByDensity = function(x1,x2,
                               ylim=c(min(x2),max(x2)),
                               xlim=c(min(x1),max(x1)),
                               xlab="",ylab="",main="") {
    
    df <- data.frame(x1,x2)
    x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")),
                  nbin=500, bandwidth=0.25)
    df$dens <- col2rgb(x)[1,] + 1L
    cols <-  colorRampPalette(c("grey75","darkorchid4","firebrick3","darkorange",
                                "gold1","yellow"), bias=1)(256)
    df$col <- cols[df$dens]
    df$opa <- normalize(df$dens)
    plot(x2~x1, data=df[order(df$dens),], 
         ylim=ylim,xlim=xlim,pch=16,col=alpha(col,1),
         cex=0.75,xlab=xlab,ylab=ylab,
         main=main)
    
    model1 <- lm(x2~x1)
    term1 <- summary(model1)
    text(5,15,labels=paste("Regression coefficient=",signif(term1$coefficients[2,1], digits=3), sep=""))
    grid()
}

# load data
a <- read.table("acr.dist.reg.frag.bed")
b <- read.table("matched_regulatory_fragments.bed")

# process data

# reads per million
con <- 1
rf <- 8.810088
rr <- 8.658067
df <- 2.950587
dr <- 2.909018

a$V12 <- a$V12/rf
a$V13 <- a$V13/rr
a$V14 <- a$V14/df
a$V15 <- a$V15/dr

b$V8 <- b$V8/rf
b$V9 <- b$V9/rr
b$V10 <- b$V10/df
b$V11 <- b$V11/dr

a$logR <- log2(a$V12+a$V13+con)
a$logD <- log2(a$V14+a$V15+con)
a$ratio <- log2((a$V12+a$V13+con)/(a$V14+a$V15+con))

b$logR <- log2(b$V8+b$V9+con)
b$logD <- log2(b$V10+b$V11+con)
b$ratio <- log2((b$V8+b$V9+con)/(b$V10+b$V10+con))

# split by acr type
a$ids <- paste(a$V1,a$V2,a$V3, sep="_")
b$ids <- paste(b$V1,b$V2,b$V3, sep="_")
acrs <- aggregate(list(R1=a$V12, R2=a$V13, D1=a$V14, D2=a$V15), 
                   by=list(sites=a$ids), FUN=sum)
cntr <- aggregate(list(R1=b$V8, R2=b$V9, D1=b$V10, D2=b$V11), 
                  by=list(sites=b$ids), FUN=sum)
dists <- aggregate(list(dist=a$V4), by=list(sites=a$ids), FUN=mean)
acrs$dists <- dists$dist
acrs$logR <- log2(acrs$R1+acrs$R2+con)
acrs$logD <- log2(acrs$D1+acrs$D2+con)
cntr$logR <- log2(cntr$R1+cntr$R2+con)
cntr$logD <- log2(cntr$D1+cntr$D2+con)
dacrs <- subset(acrs, abs(acrs$dist)>2000)
pacrs <- subset(acrs, abs(acrs$dist)<2000 & abs(acrs$dist)>0)
p500acrs <- subset(acrs, acrs$dist > 0 & acrs$dist <= 200)
p1kacrs <- subset(acrs, acrs$dist > 200 & acrs$dist < 800)
p2kacrs <- subset(acrs, acrs$dist > 800 & acrs$dist < 2000)
gacrs <- subset(acrs, acrs$dist == 0)

## make plot ##
layout(matrix(c(1:3), nrow=1, byrow=T))
maxx <- c(0,max(acrs$logD,acrs$logR,cntr$logD,cntr$logR))
max2 <- max(acrs$logD,acrs$logR,cntr$logD,cntr$logR)+1

# distal
plot_colorByDensity(dacrs$logD, dacrs$logR, 
                    xlim=maxx, 
                    ylim=maxx,
                    xlab="log2(Input[FPM]+1)",
                    ylab="log2(RNA[FPM]+1)", main="dACRs")
segments(-2,-2,max2,max2, col="black", lwd=1)


# proximal
plot_colorByDensity(pacrs$logD, pacrs$logR, 
                    xlim=maxx, 
                    ylim=maxx,
                    xlab="log2(Input[FPM]+1)",
                    ylab="log2(RNA[FPM]+1)", main="pACRs")
segments(-2,-2,max2,max2, col="black", lwd=1)

# genic
plot_colorByDensity(cntr$logD, cntr$logR, 
                    xlim=maxx, 
                    ylim=maxx,
                    xlab="log2(Input[FPM]+1)",
                    ylab="log2(RNA[FPM]+1)", main="integenic controls")
segments(-2,-2,max2,max2, col="black", lwd=1)

