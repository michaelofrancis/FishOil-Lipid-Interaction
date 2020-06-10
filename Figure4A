##This script generates figure 4A. 
##Written by Kaixiong Ye

Df  <- paste("GJB2.eQTL", ".", "all_tissues.hg19", ".pdf", sep="", collapse="");
pdf(file=Df, w=6,h=7.5)
par( mar=c(2.1,4.1,2.1,1.1), mgp=c(2, 0.8, 0))
layout(m=matrix(1:3, 3, 1))
## eQTL

tissue <- read.table("GJB2_eQTLs.all_tissues.for_plot.hg19.txt", header=TRUE)
plot(tissue[,2]/10^6, -log10(tissue[,3])+5, pch=20, xlab="Chr13", ylab="-log10(p)", bty="n", yaxt="n", ylim=c(0,35), main="", col="black")
axis(side=2, at=c(5,15,25,35), labels=c("0","10","20","30"))

arrows(x0=20761609/10^6, y0=(1.5), x1=20767037/10^6, y1=(1.5), lwd=2, code=1, length=0.02)
text(x=20767037/10^6+0.1, y=(1.5), labels="GJB2", cex=1.2);


plot(tissue[,2]/10^6, -log10(tissue[,3])+5, pch=20, xlab="Chr13", ylab="-log10(p)", bty="n", yaxt="n", ylim=c(0,35), xlim=c(20.72,20.82), main="", col="black")
axis(side=2, at=c(5,15,25,35), labels=c("0","10","20","30"))

arrows(x0=20761609/10^6, y0=(1.5), x1=20767037/10^6, y1=(1.5), lwd=2, code=1, length=0.05)
text(x=20767037/10^6+0.005, y=(1.5), labels="GJB2", cex=1.2);


data <- read.table("GJB2_eQTLs.all_tissues.100Kb.hg19.txt", header=TRUE, sep="\t")

plot(data[,1]/10^6, -log10(data[,4])+5, pch=20, xlab="Chr13", ylab="-log10(p)", bty="n", yaxt="n", ylim=c(0,35), xlim=c(20.72,20.82), main="", col="black")
axis(side=2, at=c(5,15,25,35), labels=c("0","10","20","30"))

#arrows(x0=20761609/10^6, y0=(1.5), x1=20767037/10^6, y1=(1.5), lwd=2, code=1, length=0.05)
#text(x=20767037/10^6+0.005, y=(1.5), labels="GJB2", cex=1.2);

## Draw  CDS     20763043        20763720

rect(xleft=c(20763043/10^6, 20763043/10^6), ybottom=c(0.5, 0.5), xright=c(20763720/10^6, 20763720/10^6), ytop=c(2.5, 2.5), col="black", border="black")

## Draw UTR
#UTR     20766922        20767037
#UTR     20763721        20763742
#UTR     20761609        20763039

rect(xleft=c(20766922/10^6, 20766922/10^6), ybottom=c(0.9, 0.9), xright=c(20767037/10^6, 20767037/10^6), ytop=c(2.1, 2.1), col="black", border="black")
rect(xleft=c(20763721/10^6, 20763721/10^6), ybottom=c(0.9, 0.9), xright=c(20763742/10^6, 20763742/10^6), ytop=c(2.1, 2.1), col="black", border="black")
rect(xleft=c(20761609/10^6, 20761609/10^6), ybottom=c(0.9, 0.9), xright=c(20763039/10^6, 20763039/10^6), ytop=c(2.1, 2.1), col="black", border="black")
## Draw Introns
rect(xleft=c(20761609/10^6, 20761609/10^6), ybottom=c(1.4, 1.4), xright=c(20767037/10^6, 20767037/10^6), ytop=c(1.6, 1.6), col="black", border="black")
arrows(x0=20767037/10^6-0.005, y0=(4), x1=20767037/10^6, y1=(4), lwd=1, code=1, length=0.05)
text(x=20767037/10^6+0.0045, y=(1.5), labels="GJB2", cex=1.2);



adipose_idx = data[,5] == "Adipose - Subcutaneous"
brain_idx   = data[,5] == "Brain - Cerebellum"
skin_idx    = data[,5] == "Skin - Sun Exposed (Lower leg)"
muscle_idx  = data[,5] == "Muscle - Skeletal"
Hypothalamus= data[,5] == "Brain - Hypothalamus"
fibroblast_idx = data[,5] == "Cells - Cultured fibroblasts"
thyroid_idx = data[,5] == "Thyroid"

points(data[adipose_idx,1]/10^6, -log10(data[adipose_idx,4])+5, pch=20, col="red")
points(data[brain_idx,1]/10^6, -log10(data[brain_idx,4])+5, pch=20, col="blue")
points(data[skin_idx,1]/10^6, -log10(data[skin_idx,4])+5, pch=20, col="dodgerblue")
points(data[muscle_idx,1]/10^6, -log10(data[muscle_idx,4])+5, pch=20, col="orange")
points(data[Hypothalamus,1]/10^6, -log10(data[Hypothalamus,4])+5, pch=20, col="orchid2")
points(data[fibroblast_idx,1]/10^6, -log10(data[fibroblast_idx,4])+5, pch=20, col="purple")
points(data[thyroid_idx,1]/10^6, -log10(data[thyroid_idx,4])+5, pch=20, col="coral2")

allcolors = c("red","blue", "dodgerblue","orange");
legend("topright", legend=c("Adipose - Subcutaneous (39)", "Brain - Cerebellum (37)","Skin - Sun Exposed (Lower leg) (32)", "Muscle - Skeletal (9)"), col=allcolors, bty="n", pch=16, text.col =allcolors)

allcolors = c("orchid2","purple","coral2","black");
legend("topleft", legend=c("Brain - Hypothalamus (6)", "Cells - Cultured fibroblasts (4)", "Thyroid (4)", "Other Tissues (3)"), col=allcolors, bty="n", pch=16, text.col =allcolors)

GWAS_hit <- read.table("significant_GWAS_hit.plot.hg19.txt", header=TRUE)
hit_idx <- data[,1] %in% GWAS_hit[,3]
points(data[hit_idx,1]/10^6, -log10(data[hit_idx,4])+5, pch=0, cex=1.2, col="black")
legend("top", legend=c("Interaction hits"), col="black", bty="n", pch=0, text.col="black")

text(x=20768144/10^6+0.0052, y=(-log10(2.6e-25)+5), labels="rs7987144", cex=0.8, col="black");
text(x=20790451/10^6+0.006, y=(-log10(7.7e-14)+5), labels="rs112803755", cex=0.8, col="black")


dev.off()
