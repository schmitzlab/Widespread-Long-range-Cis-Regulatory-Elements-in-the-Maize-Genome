setwd("~/Desktop/sapelo2/atac-starr-seq_maize/v5/directional_enhancer/")

# load libraries
library(MASS)
library(viridis)
library(scales)
library(densityClust)
library(RColorBrewer)
library(mgcv)

### functions ###
normalize <- function(x){
    (x - min(x))/(max(x)-min(x))
}
plot_colorByDensity = function(x1,x2, scores=NULL,
                               ylim=c(min(x2, na.rm=T),max(x2, na.rm=T)),
                               xlim=c(min(x1, na.rm=T),max(x1, na.rm=T)),
                               xlab="",ylab="",main="") {
    
    df <- data.frame(x1,x2,scores)
    x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")),
                  nbin=500, bandwidth=c(0.25,0.025))
    df$dens <- col2rgb(x)[1,] + 1L
    cols <-  colorRampPalette(c("grey75","darkorchid4","firebrick3","darkorange",
                                "gold1","yellow"), bias=1)(256)
    #cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
    #df$col <- cols[df$dens]
    if(!is.null(scores)){
        #df$col <- ifelse(scores < 0.01, "darkorange", "grey75")
        df$col <- rep("grey75", nrow(df))
    }else{
        df$opa <- normalize(df$dens)
        df$col <- cols[df$dens]
    }
    plot(x2~x1, data=df[order(df$scores, decreasing=T),], 
         ylim=ylim,xlim=xlim,pch=16,col=alpha(col,1),
         cex=1,xlab=xlab,ylab=ylab,
         main=main)
    
    # fit GAM
    df <- df[order(df$x1, decreasing=F),]
    fit <- gam(x2~s(x1), data=df, method="REML")
    print(summary(fit))
    
    # predict values
    newx <- data.frame(x1=seq(min(x1), max(x1), length.out=length(x1)))
    pred <- predict(fit, newx, interval="confidence", level=0.95, se.fit=T)
    
    # plot
    lines(newx$x1, pred$fit, col="darkorange", lwd=2)
    lines(newx$x1, pred$se.fit+pred$fit, col="darkorange", lty=2)
    lines(newx$x1, pred$fit-pred$se.fit, col="darkorange", lty=2)
    
}
fisher.combine <- function(x){
    x2 <- -2*sum(log(x))
    df <- length(x)
    pval <- pchisq(x2, 2*df, lower.tail=F)
    pval
}

# load beta binomial data
a <- read.table("distal.betabinomial.bed", header=F)
b <- read.table("proximal.betabinomial.bed", header=F)
d <- read.table("genic.betabinomial.bed", header=F)
c <- read.table("matched_regions.betabinomial.bed", header=F)
header <- c("dACRchr", "dACRstart", "dACRend", "dACRactivity", 
            "FRAGchr", "FRAGstart", "FRAGend",
            "RNA_F", "RNA_R", "INPUT_F", "INPUT_R", 
            "Empirical_bias", "Sim_bias", 
            "BinomialTest", "Prob", "Sim_Prob")
colnames(a) <- header
colnames(b) <- header
colnames(c) <- header
colnames(d) <- header

# filter for input counts
a$input <- a$INPUT_F+a$INPUT_R
b$input <- b$INPUT_F+b$INPUT_R
c$input <- c$INPUT_F+c$INPUT_R
d$input <- d$INPUT_F+d$INPUT_R

#a <- subset(a, a$input>1)
#b <- subset(b, b$input>1)
#c <- subset(c, c$input>1)

a$act <- log2(a$dACRactivity+1)
b$act <- log2(b$dACRactivity+1)
c$act <- log2(c$dACRactivity+1)

# plot layout
layout(matrix(c(1:3), nrow=1))

# plot dACR
plot_colorByDensity(a$act,a$Sim_bias, main="dACRs",
                    ylab="RNA F/R - input F/R",
                    xlab="Regulatory activity",
                    ylim=c(0, max(a$Sim_bias,b$Sim_bias,c$Sim_bias)),
                    xlim=c(0, max(a$act, b$act, c$act)))

# plot pACR
plot_colorByDensity(b$act,b$Sim_bias, main="pACRs",
                    ylim=c(0, max(a$Sim_bias,b$Sim_bias,c$Sim_bias)),
                    xlim=c(0, max(a$act, b$act, c$act)))

# plot gACR
plot_colorByDensity(c$act,c$Sim_bias, main="gACRs",
                    ylim=c(0, max(a$Sim_bias,b$Sim_bias, c$Sim_bias)),
                    xlim=c(0, max(a$act, b$act, c$act)))

# average
a$id <- paste(a$dACRchr,a$dACRstart,a$dACRend,a$dACRactivity, sep="_")
b$id <- paste(b$dACRchr,b$dACRstart,b$dACRend,b$dACRactivity, sep="_")
c$id <- paste(c$dACRchr,c$dACRstart,c$dACRend,c$dACRactivity, sep="_")
aa <- aggregate(list(rnaF=a$RNA_F, rnaR=a$RNA_R,inputF=a$INPUT_F,inputR=a$INPUT_R), 
                by=list(id=a$id), FUN=sum)
aa1 <- aggregate(list(Sim_bias=a$Sim_bias), by=list(id=a$id), FUN=mean)
aa2 <- aggregate(list(prob=a$Sim_Prob), by=list(id=a$id), FUN=fisher.combine)
bb <- aggregate(list(rnaF=b$RNA_F, rnaR=b$RNA_R,inputF=b$INPUT_F,inputR=b$INPUT_R), 
                by=list(id=b$id), FUN=sum)
bb1 <- aggregate(list(Sim_bias=b$Sim_bias), by=list(id=b$id), FUN=mean)
bb2 <- aggregate(list(prob=b$Sim_Prob), by=list(id=b$id), FUN=fisher.combine)
cc <- aggregate(list(rnaF=c$RNA_F, rnaR=c$RNA_R,inputF=c$INPUT_F,inputR=c$INPUT_R), 
                by=list(id=c$id), FUN=sum)
cc1 <- aggregate(list(Sim_bias=c$Sim_bias), by=list(id=c$id), FUN=mean)
cc2 <- aggregate(list(prob=c$Sim_Prob), by=list(id=c$id), FUN=fisher.combine)
testbb <- data.frame(do.call('rbind', strsplit(as.character(bb$id),'_',fixed=TRUE)))
testaa <- data.frame(do.call('rbind', strsplit(as.character(aa$id),'_',fixed=TRUE)))
testcc <- data.frame(do.call('rbind', strsplit(as.character(cc$id),'_',fixed=TRUE)))
aa$act <- as.numeric(as.character(testaa$X4))
bb$act <- as.numeric(as.character(testbb$X4))
cc$act <- as.numeric(as.character(testcc$X4))
aa$bias <- aa1$Sim_bias
bb$bias <- bb1$Sim_bias
cc$bias <- cc1$Sim_bias
aa$actl2 <- log2(aa$act+1)
bb$actl2 <- log2(bb$act+1)
cc$actl2 <- log2(cc$act+1)
aa$prob <- p.adjust(aa2$prob, method="fdr")
bb$prob <- p.adjust(bb2$prob, method="fdr")
cc$prob <- p.adjust(cc2$prob, method="fdr")


layout(matrix(c(1:3), nrow=1))
plot_colorByDensity(aa$actl2,aa$bias, main="dACRs", scores=aa$prob,
                    ylab="RNA F/R - input F/R",
                    xlab="Regulatory activity",
                    ylim=c(0, max(aa$bias,bb$bias,cc$bias)),
                    xlim=c(0, max(aa$actl2, bb$actl2, cc$actl2)))
plot_colorByDensity(bb$actl2,bb$bias, main="pACRs", scores=bb$prob,
                    ylim=c(0, max(aa$bias,bb$bias,cc$bias)),
                    xlim=c(0, max(aa$actl2, bb$actl2, cc$actl2)))
plot_colorByDensity(cc$actl2,cc$bias, main="gACRs", scores=cc$prob,
                    ylim=c(0, max(aa$bias,bb$bias,cc$bias)),
                    xlim=c(0, max(aa$actl2, bb$actl2, cc$actl2)))












# quick plot
layout(matrix(c(1:1), nrow=1))

breaks <- c(0, 0.05, 1)
labels <- c("<5e-2", ">5e-2")
a$abins <- cut(a$BinomialTest, breaks, include.lowest = T, right=FALSE, labels=labels)
b$bbins <- cut(b$BinomialTest, breaks, include.lowest = T, right=FALSE, labels=labels)
aprop <- prop.table(table(a$abins))
bprop <- prop.table(table(b$bbins))
both <- data.frame(dACRs=aprop, input_peaks=bprop)
both$input_peaks.Var1 <- NULL
colnames(both) <- c("q-value range", "dACRs", "Input_peaks")
rownames(both) <- both$`q-value range`
both$`q-value range` <- NULL
barplot(t(as.matrix(both)), beside=T, ylim=c(0,1), xlab="Fragment bias (p-value bin)",
        ylab="Proportion")





#################################################################################
### Bias heat map 
#################################################################################

# load functions
revrow <- function(x){
    len <- length(x)
    half <- as.integer(len/2)
    left <- sum(x[1:half], na.rm=T)
    right <- sum(x[(1+half):len], na.rm=T)
    if(left > right){
        x <- rev(x)    
    }else{
        x <- x
    }
    return(x)
}

revran <- function(x){
    ran <- runif(1, 0, 1)
    if(ran > 0.5){
        x <- rev(x)    
    }else{
        x <- x
    }
    return(x)
}

modscale <- function(x){
    t1 <- ((4)*(x-0.5))+(-1) 
}

modscale0 <- function(x){
    t1 <- (((1-0)/(max(x, na.rm=T)-min(x,na.rm=T)))*(x-max(x, na.rm=T))) + 1
}

normalizeM <- function(x){
    quant1 <- quantile(x, 0.01)
    quant2 <- quantile(x, 0.99)
    x[x<quant1] <- quant1
    x[x>quant2] <- quant2
    return(x)
}

# load data
adj <- 0.0227
drnaf <- as.matrix(read.table("peakDistalACRs.bed.RNA_F.mat"))
drnar <- as.matrix(read.table("peakDistalACRs.bed.RNA_R.mat"))
dinpf <- as.matrix(read.table("peakDistalACRs.bed.Input_F.mat"))
dinpr <- as.matrix(read.table("peakDistalACRs.bed.Input_R.mat"))

crnaf <- as.matrix(read.table("matched_regions.sort.bed.RNA_F.mat"))
crnar <- as.matrix(read.table("matched_regions.sort.bed.RNA_R.mat"))
crnaf <- t(apply(crnaf, 1, revran))
crnar <- t(apply(crnar, 1, revran))

cinpf <- as.matrix(read.table("matched_regions.sort.bed.Input_F.mat"))
cinpr <- as.matrix(read.table("matched_regions.sort.bed.Input_R.mat"))

# plot column means
avedrf <- colMeans(drnaf, na.rm=T)
avedrr <- colMeans(drnar, na.rm=T)
avedif <- colMeans(dinpf, na.rm=T)
avedir <- colMeans(dinpr, na.rm=T)
avecrf <- colMeans(crnaf, na.rm=T)
avecrr <- colMeans(crnar, na.rm=T)
avecif <- colMeans(cinpf, na.rm=T)
avecir <- colMeans(cinpr, na.rm=T)

# plot
layout(matrix(c(1:2), nrow=1))
plot(avecrf, type="l", col="forestgreen", main="STARR-RNA", xaxt="n", xlab="", 
     ylim=c(min(-avedrr,-avecrr,-avedir,-avecir), 
            max(avedrf, avecrf, avedif, avecif)))
lines(avedrf, col="darkorchid4")
lines(-avecrr, col="forestgreen")
lines(-avedrr, col="darkorchid4")
abline(h=0, col="grey35", lty=2)
axis(1, at=c(0, 200, 400), labels=c("-1kb",0, "1kb"))

plot(avecif, type="l", col="forestgreen", main="STARR-input", xaxt="n", xlab="", 
     ylim=c(min(-avedrr,-avecrr,-avedir,-avecir), 
            max(avedrf, avecrf, avedif, avecif)))
lines(avedif, col="darkorchid4")
lines(-avecir, col="forestgreen")
lines(-avedir, col="darkorchid4")
abline(h=0, col="grey35", lty=2)
axis(1, at=c(0, 200, 400), labels=c("-1kb",0, "1kb"))







## for heatmap matrix
rnaf[is.na(rnaf)] <- 0
rnar[is.na(rnar)] <- 0
inpf[is.na(inpf)] <- 0
inpr[is.na(inpr)] <- 0

rnaf<-rnaf+adj
rnar<-rnar+adj
inpf<-inpf+adj
inpr<-inpr+adj

rnas <- rnaf + rnar
inps <- inpf + inpr

# process
rna <- (1-abs((rnaf/(rnar+rnaf))-0.5))
input <- (1-abs((inpf/(inpf+inpr))-0.5))
rna <- modscale(rna)
input <- modscale(input)
rnarev <- rna*(modscale0(rnas))
inputrev <- input*(modscale0(inps))

rownames(rnarev) <- seq(1:nrow(rnarev))
rownames(inputrev) <- seq(1:nrow(inputrev))
rnarev <- normalizeM(rnarev[order(rowSums(rnarev)),])
inputrev <- normalizeM(inputrev[rownames(rnarev),])

# plot parameters
nums1 <- seq(from=min(rnarev,inputrev), to=0, length.out=50)
nums2 <- seq(from=0, to=max(rnarev,inputrev), length.out=50)
nums <- c(nums1, nums2)
cols <- colorRampPalette(c("darkorange","white","darkorchid4"))(99)
layout(matrix(c(1:2), nrow=1))

# plots
image(t(rnarev), useRaster=T, breaks=nums, col=cols, xaxt="n", yaxt="n", main="dACR_rna")
image(t(inputrev), useRaster=T, breaks=nums, col=cols, xaxt="n", yaxt="n", main="dACR_input")












