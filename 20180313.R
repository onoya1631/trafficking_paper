library(exactRankTests)
library(seqinr)
library(ape)
library(gtools)
library(ggplot2)
library(e1071) 

MyGray <- rgb(t(col2rgb("black")), alpha=60, maxColorValue=255)
MyTransparent <- rgb(t(col2rgb("black")), alpha=0, maxColorValue=255)
MyRed <- rgb(t(col2rgb("red")), alpha=150, maxColorValue=255)
MyBlue <- rgb(t(col2rgb("blue")), alpha=150, maxColorValue=255)

vivo <- read.csv("20170303.vivo.csv")
oOR <- subset(vivo, vivo=="oOR")
uOR <- subset(vivo, vivo=="uOR")

facs <- head(merge(read.csv("20180313.FACS.csv"),read.csv("20170303.vivo.csv"),by = "Olfr", all=T),242)
facs <- facs[rev(order(facs$Normalized)),]

OR.aln <- read.alignment(file="20170318.facs.yue.replace_20170402.txt",format="clustal")
OR.aln$nam <- as.character(lapply(strsplit(OR.aln$nam, split="/"), "[", 1))
OR.dist <- dist.alignment(OR.aln,matrix="similarity")
OR.tree <- nj(OR.dist)

all.OR <- NULL
for (i in 1:length(OR.aln$nam)) {
  all.OR <- rbind(all.OR, strsplit(OR.aln$seq[[i]],NULL)[[1]])
};rownames(all.OR) <- OR.aln$nam

sub.facs <- subset(facs, vivo=="oOR"|vivo=="uOR")

library(mixtools); facs.p <- 0.001
hist(sub.facs$Normalized, breaks=150, xlab="Fold Change (normarized by Olfr539 and Olfr541)",ylab="density" ,freq=F)
set.seed(1);mixmodl.facs <- normalmixEM(sub.facs$Normalized);plot(mixmodl.facs,which=2,breaks=100)
facs.cut <- qnorm(p=(1-facs.p), mixmodl.facs$mu[1],mixmodl.facs$sigma[1])
hist(sub.facs$Normalized,breaks=100,freq=F,xlab="Fold Change",ylab="density",main="FACS normalized by Olfr539")
abline(v=facs.cut,col="blue"); text(1,3,paste("facs cut-off:",round(facs.cut,3)))

hist(subset(sub.facs, vivo=="uOR")$Normalized, col = MyBlue, border = MyBlue, ylim=c(0,60), breaks=c(seq(-0.1,2.6,0.05)),xlab="Fold Change (normarized by Olfr539 and Olfr541)",ylab="freq" ,freq=T, xlim=c(-0.1,2.6));abline(v=0.155, lty=2, col="blue")
hist(subset(sub.facs, vivo=="oOR")$Normalized, col = MyRed, add = T, border = MyRed, ylim=c(0,60),breaks=c(seq(-0.1,2.6,0.05)), xlab="Fold Change (normarized by Olfr539 and Olfr541)",ylab="freq" ,freq=T);abline(v=0.155, lty=2, col="blue", main="oOR")
legend("right",c("oOR","uOR"),pch=15,col=c(MyRed,MyBlue), bty="n")
wilcox.exact(subset(sub.facs, vivo=="uOR")$Normalized,subset(sub.facs, vivo=="oOR")$Normalized, paired = F)#, alternative = "g")

ind <- subset(sub.facs, vivo=="oOR"&Normalized>0.153)
dep <- subset(sub.facs, vivo=="uOR"&Normalized<0.153)

positive <- subset(sub.facs, Normalized>0.153)
negative <- subset(sub.facs, Normalized<0.153) 

##### phylogenetic tree ###################################################
plot.OR.tree <- function(OR.tree,edgecol=MyGray2,edgewid=0.5) { #label2 overrides label1
  #  length <- 1111
  #  plot(OR.tree,cex=1,type="fan",edge.color="azure4",pch=21,show.tip.label=F,tip.color="black",x.lim=c(-0.2,0.4),y.lim=c(-0.41,0.29))
  #  tiplabels(pch=21,cex=0.9,col="azure4",bg="azure3",lwd=0.7)
  plot(OR.tree,cex=1,type="unrooted",edge.color=edgecol,pch=21,show.tip.label=F,tip.color="black",edge.width=edgewid,no.margin=T)#,y.lim=c(-0.41,0.29)) x.lim=c(-0.2,0.4),y.lim=c(0.4,-0.3
  tiplabels(pch=21,cex=0.9,col=MyTransparent,bg=MyTransparent,lwd=0.7)
}
add.OR.tree.label <- function(OR.tree,label,pch=21,pt.col="azure2",bg.col="azure1",lwd=0.7,cex=1,frame="none") {
  length <- length(OR.tree$tip.label)
  pt.color <- rep(MyTransparent,length)
  bg.color <- rep(MyTransparent,length)
  tip.color <- rep(MyTransparent,length)
  for (i in 1:length) {
    if (OR.tree$tip.label[i] %in% label)
    { pt.color[i]<-pt.col;bg.color[i]<-bg.col;tip.color[i]<-bg.col }
  }   
  tiplabels(pch=pch,cex=cex,col=pt.color,bg=bg.color,lwd=lwd)
  #  tiplabels(OR.tree$tip.label,cex=cex,col=tip.color,bg=MyTransparent,frame="none")
}
##### independent vs dependent ###################################################
edge.ind <- which.edge(OR.tree,as.character(ind$Olfr))
edge.dep <- which.edge(OR.tree,as.character(dep$Olfr))

edge.clr <- rep(MyGray,length(OR.tree$edge))
edge.clr[edge.ind] <- "red"
edge.clr[edge.dep] <- "black"

edge.wdth <- rep(0.5,length(OR.tree$edge))
edge.wdth[edge.ind] <- 3
edge.wdth[edge.dep] <- 3

plot.OR.tree(OR.tree,edge.clr,edge.wdth)
add.OR.tree.label(OR.tree,c(OR.aln$nam),21,"gray31","azure3",0.5,1)
add.OR.tree.label(OR.tree,as.character(ind$Olfr),21,"gray31","red",0.5,1)
add.OR.tree.label(OR.tree,as.character(dep$Olfr),21,"gray31","black",0.5,1)
legend("topright",c("independent","dependent"), pch=21,col="gray31",pt.bg=c ("red","black"),bty="n",pt.lwd=0.5,pt.cex=1,cex=1) 
##### independent ###################################################
edge.ind <- which.edge(OR.tree,as.character(ind$Olfr))

edge.clr <- rep(MyGray,length(OR.tree$edge))
edge.clr[edge.ind] <- "red"

edge.wdth <- rep(0.5,length(OR.tree$edge))
edge.wdth[edge.ind] <- 3

plot.OR.tree(OR.tree,edge.clr,edge.wdth)
add.OR.tree.label(OR.tree,c(OR.aln$nam),21,"gray31","azure3",0.5,1)
add.OR.tree.label(OR.tree,as.character(ind$Olfr),21,"gray31","red",0.5,1)
legend("topright",c("independent"), pch=21,col="gray31",pt.bg=c ("red"),bty="n",pt.lwd=0.5,pt.cex=1,cex=1) 
##### dependent ###################################################
edge.dep <- which.edge(OR.tree,as.character(dep$Olfr))

edge.clr <- rep(MyGray,length(OR.tree$edge))
edge.clr[edge.dep] <- "black"

edge.wdth <- rep(0.5,length(OR.tree$edge))
edge.wdth[edge.dep] <- 3

plot.OR.tree(OR.tree,edge.clr,edge.wdth)
add.OR.tree.label(OR.tree,c(OR.aln$nam),21,"gray31","azure3",0.5,1)
add.OR.tree.label(OR.tree,as.character(dep$Olfr),21,"gray31","black",0.5,1)
legend("topright",c("dependent"), pch=21,col="gray31",pt.bg=c ("black"),bty="n",pt.lwd=0.5,pt.cex=1,cex=1) 
###### positive vs negative ##################################################
edge.ind <- which.edge(OR.tree,as.character(positive$Olfr))
edge.dep <- which.edge(OR.tree,as.character(negative$Olfr))

edge.clr <- rep(MyGray,length(OR.tree$edge))
edge.clr[edge.ind] <- "red"
edge.clr[edge.dep] <- "black"

edge.wdth <- rep(0.5,length(OR.tree$edge))
edge.wdth[edge.ind] <- 3
edge.wdth[edge.dep] <- 3

plot.OR.tree(OR.tree,edge.clr,edge.wdth)
add.OR.tree.label(OR.tree,c(OR.aln$nam),21,"gray31","azure3",0.5,1)
add.OR.tree.label(OR.tree,as.character(positive$Olfr),21,"gray31","red",0.5,1)
add.OR.tree.label(OR.tree,as.character(negative$Olfr),21,"gray31","black",0.5,1)
legend("topright",c("positive","negative"), pch=21,col="gray31",pt.bg=c ("red","black"),bty="n",pt.lwd=0.5,pt.cex=1,cex=1) 
###### positive ##################################################
edge.ind <- which.edge(OR.tree,as.character(positive$Olfr))

edge.clr <- rep(MyGray,length(OR.tree$edge))
edge.clr[edge.ind] <- "red"

edge.wdth <- rep(0.5,length(OR.tree$edge))
edge.wdth[edge.ind] <- 3

plot.OR.tree(OR.tree,edge.clr,edge.wdth)
add.OR.tree.label(OR.tree,c(OR.aln$nam),21,"gray31","azure3",0.5,1)
add.OR.tree.label(OR.tree,as.character(positive$Olfr),21,"gray31","red",0.5,1)
legend("topright",c("positive"), pch=21,col="gray31",pt.bg=c ("red"),bty="n",pt.lwd=0.5,pt.cex=1,cex=1) 
###### negative ##################################################
edge.dep <- which.edge(OR.tree,as.character(negative$Olfr))

edge.clr <- rep(MyGray,length(OR.tree$edge))
edge.clr[edge.dep] <- "black"

edge.wdth <- rep(0.5,length(OR.tree$edge))
edge.wdth[edge.dep] <- 3

plot.OR.tree(OR.tree,edge.clr,edge.wdth)
add.OR.tree.label(OR.tree,c(OR.aln$nam),21,"gray31","azure3",0.5,1)
add.OR.tree.label(OR.tree,as.character(negative$Olfr),21,"gray31","black",0.5,1)
legend("topright",c("negative"), pch=21,col="gray31",pt.bg=c ("black"),bty="n",pt.lwd=0.5,pt.cex=1,cex=1) 
###### oOR vs uOR ##################################################
edge.ind <- which.edge(OR.tree,as.character(oOR$Olfr))
edge.dep <- which.edge(OR.tree,as.character(uOR$Olfr))

edge.clr <- rep(MyGray,length(OR.tree$edge))
edge.clr[edge.ind] <- "red"
edge.clr[edge.dep] <- "black"

edge.wdth <- rep(0.5,length(OR.tree$edge))
edge.wdth[edge.ind] <- 3
edge.wdth[edge.dep] <- 3

plot.OR.tree(OR.tree,edge.clr,edge.wdth)
add.OR.tree.label(OR.tree,c(OR.aln$nam),21,"gray31","azure3",0.5,1)
add.OR.tree.label(OR.tree,as.character(oOR$Olfr),21,"gray31","red",0.5,1)
add.OR.tree.label(OR.tree,as.character(uOR$Olfr),21,"gray31","black",0.5,1)
legend("topright",c("oOR","uOR"), pch=21,col="gray31",pt.bg=c ("red","black"),bty="n",pt.lwd=0.5,pt.cex=1,cex=1) 
###### oOR ##################################################
edge.ind <- which.edge(OR.tree,as.character(oOR$Olfr))

edge.clr <- rep(MyGray,length(OR.tree$edge))
edge.clr[edge.ind] <- "red"

edge.wdth <- rep(0.5,length(OR.tree$edge))
edge.wdth[edge.ind] <- 3

plot.OR.tree(OR.tree,edge.clr,edge.wdth)
add.OR.tree.label(OR.tree,c(OR.aln$nam),21,"gray31","azure3",0.5,1)
add.OR.tree.label(OR.tree,as.character(oOR$Olfr),21,"gray31","red",0.5,1)
legend("topright",c("oOR"), pch=21,col="gray31",pt.bg=c ("red"),bty="n",pt.lwd=0.5,pt.cex=1,cex=1) 
###### uOR ##################################################
edge.dep <- which.edge(OR.tree,as.character(uOR$Olfr))

edge.clr <- rep(MyGray,length(OR.tree$edge))
edge.clr[edge.dep] <- "black"

edge.wdth <- rep(0.5,length(OR.tree$edge))
edge.wdth[edge.dep] <- 3

plot.OR.tree(OR.tree,edge.clr,edge.wdth)
add.OR.tree.label(OR.tree,c(OR.aln$nam),21,"gray31","azure3",0.5,1)
add.OR.tree.label(OR.tree,as.character(uOR$Olfr),21,"gray31","black",0.5,1)
legend("topright",c("uOR"), pch=21,col="gray31",pt.bg=c ("black"),bty="n",pt.lwd=0.5,pt.cex=1,cex=1) 
########################################################
##### PCA ###################################################
AAprop <- read.csv("Yue_AAprop.txt"); AAprop <- AAprop[,2:5]

dim(all.OR)
AA2prop <- function(aa,AAprop) {
  return(as.numeric(subset(AAprop,AA==aa,select=c(c,p,v))))
}

Prop.OR.raw <- NULL 
for (i in 1:dim(all.OR)[1]) {
  row <- NULL
  for (j in (1:dim(all.OR)[2])) {
    row <- c(row,AA2prop(all.OR[i,j],AAprop))
  }
  Prop.OR.raw <- rbind(Prop.OR.raw,row)
}
#Prop.OR.raw <- Prop.OR
colNam <- NULL
for (i in 1:dim(all.OR)[2]) {colNam <- c(colNam,paste(i,"c",sep="_"),paste(i,"p",sep="_"),paste(i,"v",sep="_"))}
keep.ind <- (colMeans(all.OR=="-") < 0.2)
Prop.OR.raw <- as.data.frame(Prop.OR.raw)
colnames(Prop.OR.raw) <- colNam
row.names(Prop.OR.raw) <- OR.aln$nam

Prop.OR <- Prop.OR.raw[,as.vector(sapply(keep.ind, function (x) rep(x,3)))]
for (i in 1:dim(Prop.OR)[2]) {Prop.OR[is.na(Prop.OR[,i]),i] <- mean(Prop.OR[,i],na.rm=T)}
var.ind <- apply(Prop.OR,2,var)!=0
Prop.OR <- Prop.OR[var.ind]
#OR.pca <- prcomp(~ ., data = prop.OR, na.action = na.omit)
OR.pca <- prcomp(Prop.OR,scale.=T)
plot(OR.pca$x[,1],OR.pca$x[,2],pch=3,cex=0.7)
summary(OR.pca)$importance[2:3,]

library(rgl)
col=rep("gray",length(OR.pca$x[,1]))
col[which(OR.aln$nam%in%ind$Olfr)] <- "red"
col[which(OR.aln$nam%in%dep$Olfr)] <- "black"
plot3d(OR.pca$x[,1],OR.pca$x[,2],OR.pca$x[,3],col=col,pch=16,xlab="PC1",ylab="PC2",zlab="PC3",size=5)
rgl.postscript("pca.pdf","pdf")
pdf("legend.pdf",5,5)
#plot.new();legend("bottomright",pch=16,col=c("cyan","magenta","gray"),c("Activated","Not activated","Not tested"),bty = "n")
plot.new();legend("bottomright",pch=16,col=c("red","black","gray"),c("RTP-independent ORs","RTP-dependent ORs","Undetermined"),bty = "n")
########################################################
##### distance ###################################################
keep <- NULL; for(i in 1:dim(all.OR)[2]){keep <- c(keep, mean(all.OR[,i]=="-")<0.2)}

ind.aa <- subset(all.OR, rownames(all.OR)%in%ind$Olfr)
dep.aa <- subset(all.OR, rownames(all.OR)%in%dep$Olfr)

grantham <- read.table("grantham2.txt"); colnames(grantham) <- rownames(grantham)

func.dist <- function(aa){
  idx <- combinations(nrow(aa),2)
  t1 <- NULL
  for(j in 1:ncol(aa)){
    t2 <- NULL
    for(i in 1:nrow(idx)){
      c1 <- idx[i,1]; c2 <- idx[i,2] 
      t2 <- c(t2,grantham[which(rownames(grantham)==aa[c1,j]),which(colnames(grantham)==aa[c2,j])])
    }; t1 <- cbind(t1,t2)
  };  return(t1)  
}
ind.dist <- func.dist(ind.aa)
dep.dist <- func.dist(dep.aa)

mix.i <- data.frame(ind=rep(1:26, each=nrow(dep.aa)), dep=rep(1:nrow(dep.aa), 26)); dist.mix <- NULL
for(j in 1:ncol(all.OR)){
  dist.mix.1 <- NULL
  for(i in 1:nrow(mix.i)){
    c1 <- mix.i[i,1]; c2 <- mix.i[i,2] 
    dist.mix.1 <- c(dist.mix.1,grantham[which(rownames(grantham)==ind.aa[c1,j]),which(colnames(grantham)==dep.aa[c2,j])])
  }; dist.mix <- cbind(dist.mix,dist.mix.1)
}; mix.dist <- dist.mix

all.dist <- rbind(mix.dist, ind.dist, dep.dist)

pval.ind <- NULL; for(i in which(keep)){pval.ind <- c(pval.ind,as.numeric(t.test(all.dist[,i], ind.dist[,i], alternative = "g")[3]))}
pval.dep <- NULL; for(i in which(keep)){pval.dep <- c(pval.dep,as.numeric(t.test(all.dist[,i], dep.dist[,i], alternative = "g")[3]))}

########################################################
##### consensus ###################################################
consensus <- NULL; for (i in 1:ncol(all.OR)){k1 <- NULL; for (k in 1:20){k1 <- c(k1,sum(all.OR[,i]==rownames(grantham)[k]))}; consensus <- c(consensus,rownames(grantham)[which.max(k1)])}
freq.cons <- NULL; for(i in 1:ncol(all.OR)){freq.cons <- c(freq.cons, (sum(all.OR[,i]==consensus[i])/sum(all.OR[,i]!="-"))*100)}
freq.cons.ind <- NULL; for(i in 1:ncol(ind.aa)){freq.cons.ind <- c(freq.cons.ind, (sum(ind.aa[,i]==consensus[i])/sum(ind.aa[,i]!="-"))*100)}
freq.cons.dep <- NULL; for(i in 1:ncol(dep.aa)){freq.cons.dep <- c(freq.cons.dep, (sum(dep.aa[,i]==consensus[i])/sum(dep.aa[,i]!="-"))*100)}

par(mar=c(4,3,3,2), mgp=c(2,0.8,0))

plot(1, type="n", xlim=c(1,66), ylim=c(0,100), bty="n")
points(1:66, freq.cons[keep][which(pval.ind<(0.05/307))], type="l", col="black")
points(1:66, freq.cons.dep[keep][which(pval.ind<(0.05/307))], type="l", col="blue")
points(1:66, freq.cons.ind[keep][which(pval.ind<(0.05/307))], type="l", col="red")

pdf("consensus.pdf")
od <- rev(order(freq.cons[keep][which(pval.ind<(0.05/307))]))
plot(1, type="n", xlim=c(1,66), ylim=c(0,100), bty="n")
points(1:66, freq.cons[keep][which(pval.ind<(0.05/307))][od], type="l", col="black")
points(1:66, freq.cons.dep[keep][which(pval.ind<(0.05/307))][od], type="l", col="blue")
points(1:66, freq.cons.ind[keep][which(pval.ind<(0.05/307))][od], type="l", col="red")
dev.off()

write.csv(all.OR[,keep][,which(pval.ind<0.05/307)][,od],"./logo/allOR.66sites.csv")
########################################################
##### MD ###################################################
Olfr539.allTM.MD1 <- read.table("./MD/20170717/Olfr539/MD1/ken_rmsF-allTM-Olfr539-MD1.xm")[1:1050,2]
Olfr539.allTM.MD2 <- read.table("./MD/20170717/Olfr539/MD2/ken_rmsF-allTM-Olfr539-MD2.xm")[1:1050,2]
Olfr539.allTM.MD3 <- read.table("./MD/20170717/Olfr539/MD3/ken_rmsF-allTM-Olfr539-MD3.xm")[1:1050,2]
Olfr539.allTM.MD4 <- read.table("./MD/20170717/Olfr539/MD4/ken_rmsF-allTM-Olfr539-MD4.xm")[1:1050,2]
Olfr539.allTM.MD5 <- read.table("./MD/20170717/Olfr539/MD5/ken_rmsF-allTM-Olfr539-MD5.xm")[1:1050,2]
Olfr539.allTM.MD6 <- read.table("./MD/20170717/Olfr539/MD6/ken_rmsF-allTM-Olfr539-MD6.xm")[1:1050,2]

Olfr539G154C.allTM.MD1 <- read.table("./MD/20170717/Olfr539G154C/MD1/ken_rmsF-allTM-Olfr539C-MD1.xm")[1:1050,2]
Olfr539G154C.allTM.MD2 <- read.table("./MD/20170717/Olfr539G154C/MD2/ken_rmsF-allTM-Olfr539C-MD2.xm")[1:1050,2]
Olfr539G154C.allTM.MD3 <- read.table("./MD/20170717/Olfr539G154C/MD3/ken_rmsF-allTM-Olfr539C-MD3.xm")[1:1050,2]
Olfr539G154C.allTM.MD4 <- read.table("./MD/20170717/Olfr539G154C/MD4/ken_rmsF-allTM-Olfr539C-MD4.xm")[1:1050,2]
Olfr539G154C.allTM.MD5 <- read.table("./MD/20170717/Olfr539G154C/MD5/ken_rmsF-allTM-Olfr539C-MD5.xm")[1:1050,2]
Olfr539G154C.allTM.MD6 <- read.table("./MD/20170717/Olfr539G154C/MD6/ken_rmsF-allTM-Olfr539C-MD6.xm")[1:1050,2]

Olfr539V209G.allTM.MD1 <- read.table("./MD/20170717/Olfr539V209G/MD1/ken_rmsF-allTM-Olfr539V209G-MD1.xm")[1:1050,2]
Olfr539V209G.allTM.MD2 <- read.table("./MD/20170717/Olfr539V209G/MD2/ken_rmsF-allTM-Olfr539V209G-MD2.xm")[1:1050,2]
Olfr539V209G.allTM.MD3 <- read.table("./MD/20170717/Olfr539V209G/MD3/ken_rmsF-allTM-Olfr539V209G-MD3.xm")[1:1050,2]
Olfr539V209G.allTM.MD4 <- read.table("./MD/20170717/Olfr539V209G/MD4/ken_rmsF-allTM-Olfr539V209G-MD4.xm")[1:1050,2]
Olfr539V209G.allTM.MD5 <- read.table("./MD/20170717/Olfr539V209G/MD5/ken_rmsF-allTM-Olfr539V209G-MD5.xm")[1:1050,2]
Olfr539V209G.allTM.MD6 <- read.table("./MD/20170717/Olfr539V209G/MD6/ken_rmsF-allTM-Olfr539V209G-MD6.xm")[1:1050,2]

Olfr539L155A.allTM.MD1 <- read.table("./MD/20170717/Olfr539L155A/MD1/ken_rmsF-allTM-Olfr539L155A-MD1.xm")[1:1050,2]
Olfr539L155A.allTM.MD2 <- read.table("./MD/20170717/Olfr539L155A/MD2/ken_rmsF-allTM-Olfr539L155A-MD2.xm")[1:1050,2]
Olfr539L155A.allTM.MD3 <- read.table("./MD/20170717/Olfr539L155A/MD3/ken_rmsF-allTM-Olfr539L155A-MD3.xm")[1:1050,2]
Olfr539L155A.allTM.MD4 <- read.table("./MD/20170717/Olfr539L155A/MD4/ken_rmsF-allTM-Olfr539L155A-MD4.xm")[1:1050,2]
Olfr539L155A.allTM.MD5 <- read.table("./MD/20170717/Olfr539L155A/MD5/ken_rmsF-allTM-Olfr539L155A-MD5.xm")[1:1050,2]
Olfr539L155A.allTM.MD6 <- read.table("./MD/20170717/Olfr539L155A/MD6/ken_rmsF-allTM-Olfr539L155A-MD6.xm")[1:1050,2]

Olfr541.allTM.MD1 <- read.table("./MD/20170717/Olfr541/MD1/ken_rmsF-allTM-Olfr541-MD1.xm")[1:1050,2]
Olfr541.allTM.MD2 <- read.table("./MD/20170717/Olfr541/MD2/ken_rmsF-allTM-Olfr541-MD2.xm")[1:1050,2]
Olfr541.allTM.MD3 <- read.table("./MD/20170717/Olfr541/MD3/ken_rmsF-allTM-Olfr541-MD3.xm")[1:1050,2]
Olfr541.allTM.MD4 <- read.table("./MD/20170717/Olfr541/MD4/ken_rmsF-allTM-Olfr541-MD4.xm")[1:1050,2]
Olfr541.allTM.MD5 <- read.table("./MD/20170717/Olfr541/MD5/ken_rmsF-allTM-Olfr541-MD5.xm")[1:1050,2]
Olfr541.allTM.MD6 <- read.table("./MD/20170717/Olfr541/MD6/ken_rmsF-allTM-Olfr541-MD6.xm")[1:1050,2]

Olfr541C154G.allTM.MD1 <- read.table("./MD/20170717/Olfr541C154G/MD1/ken_rmsF-allTM-Olfr541C154G-MD1.xm")[1:1050,2]
Olfr541C154G.allTM.MD2 <- read.table("./MD/20170717/Olfr541C154G/MD2/ken_rmsF-allTM-Olfr541C154G-MD2.xm")[1:1050,2]
Olfr541C154G.allTM.MD3 <- read.table("./MD/20170717/Olfr541C154G/MD3/ken_rmsF-allTM-Olfr541C154G-MD3.xm")[1:1050,2]
Olfr541C154G.allTM.MD4 <- read.table("./MD/20170717/Olfr541C154G/MD4/ken_rmsF-allTM-Olfr541C154G-MD4.xm")[1:1050,2]
Olfr541C154G.allTM.MD5 <- read.table("./MD/20170717/Olfr541C154G/MD5/ken_rmsF-allTM-Olfr541C154G-MD5.xm")[1:1050,2]
Olfr541C154G.allTM.MD6 <- read.table("./MD/20170717/Olfr541C154G/MD6/ken_rmsF-allTM-Olfr541C154G-MD6.xm")[1:1050,2]

Olfr541GV.allTM.MD1 <- read.table("./MD/20170717/Olfr541GV/MD1/ken_rmsF-allTM-Olfr541GV-MD1.xm")[1:1050,2]
Olfr541GV.allTM.MD2 <- read.table("./MD/20170717/Olfr541GV/MD2/ken_rmsF-allTM-Olfr541GV-MD2.xm")[1:1050,2]
Olfr541GV.allTM.MD3 <- read.table("./MD/20170717/Olfr541GV/MD3/ken_rmsF-allTM-Olfr541GV-MD3.xm")[1:1050,2]
Olfr541GV.allTM.MD4 <- read.table("./MD/20170717/Olfr541GV/MD4/ken_rmsF-allTM-Olfr541GV-MD4.xm")[1:1050,2]
Olfr541GV.allTM.MD5 <- read.table("./MD/20170717/Olfr541GV/MD5/ken_rmsF-allTM-Olfr541GV-MD5.xm")[1:1050,2]
Olfr541GV.allTM.MD6 <- read.table("./MD/20170717/Olfr541GV/MD6/ken_rmsF-allTM-Olfr541GV-MD6.xm")[1:1050,2]

m <- c(Olfr539.allTM.MD1,Olfr539.allTM.MD2,Olfr539.allTM.MD3,Olfr539.allTM.MD4,Olfr539.allTM.MD5,Olfr539.allTM.MD6,
       Olfr539G154C.allTM.MD1,Olfr539G154C.allTM.MD2,Olfr539G154C.allTM.MD3,Olfr539G154C.allTM.MD4,Olfr539G154C.allTM.MD5,Olfr539G154C.allTM.MD6,
       Olfr539V209G.allTM.MD1,Olfr539V209G.allTM.MD2,Olfr539V209G.allTM.MD3,Olfr539V209G.allTM.MD4,Olfr539V209G.allTM.MD5,Olfr539V209G.allTM.MD6,
       Olfr539L155A.allTM.MD1,Olfr539L155A.allTM.MD2,Olfr539L155A.allTM.MD3,Olfr539L155A.allTM.MD4,Olfr539L155A.allTM.MD5,Olfr539L155A.allTM.MD6)
       
par(mar=c(2, 2, 1, 1), xpd = TRUE, bty="n", mgp = c(1,0.5,0))
plot(1,xlim=c(0,28), ylim=c(0,max(m)+1),axes = FALSE, type="n", xlab=NA, ylab="RMSD")
axis(1, labels = NA, xaxt="n"); axis(2, las=1)
boxplot(Olfr539.allTM.MD1[950:1050], at = 1, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539.allTM.MD2[950:1050], at = 2, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539.allTM.MD3[950:1050], at = 3, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539.allTM.MD4[950:1050], at = 4, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539.allTM.MD5[950:1050], at = 5, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539.allTM.MD6[950:1050], at = 6, add = TRUE, yaxt  = "n", width = 1.5)

boxplot(Olfr539L155A.allTM.MD1[950:1050], at = 8, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539L155A.allTM.MD2[950:1050], at = 9, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539L155A.allTM.MD3[950:1050], at = 10, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539L155A.allTM.MD4[950:1050], at = 11, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539L155A.allTM.MD5[950:1050], at = 12, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539L155A.allTM.MD6[950:1050], at = 13, add = TRUE, yaxt  = "n", width = 1.5)

boxplot(Olfr539V209G.allTM.MD1[950:1050], at = 15, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539V209G.allTM.MD2[950:1050], at = 16, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539V209G.allTM.MD3[950:1050], at = 17, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539V209G.allTM.MD4[950:1050], at = 18, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539V209G.allTM.MD5[950:1050], at = 19, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539V209G.allTM.MD6[950:1050], at = 20, add = TRUE, yaxt  = "n", width = 1.5)

boxplot(Olfr539G154C.allTM.MD1[950:1050], at = 22, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539G154C.allTM.MD2[950:1050], at = 23, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539G154C.allTM.MD3[950:1050], at = 24, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539G154C.allTM.MD4[950:1050], at = 25, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539G154C.allTM.MD5[950:1050], at = 26, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr539G154C.allTM.MD6[950:1050], at = 27, add = TRUE, yaxt  = "n", width = 1.5)
segments(1,0,6,0, lwd = 2);segments(8,0,13,0, lwd = 2);segments(15,0,20,0, lwd = 2);segments(22,0,27,0, lwd = 2);
text(3.5,-1,"Olfr539");text(10.5,-1,"L155A");text(17.5,-1,"V209G");text(24.5,-1,"G154C");

m <- c(Olfr541.allTM.MD1,Olfr541.allTM.MD2,Olfr541.allTM.MD3,Olfr541.allTM.MD4,Olfr541.allTM.MD5,Olfr541.allTM.MD6,
       Olfr541C154G.allTM.MD1,Olfr541C154G.allTM.MD2,Olfr541C154G.allTM.MD3,Olfr541C154G.allTM.MD4,Olfr541C154G.allTM.MD5,Olfr541C154G.allTM.MD6,
       Olfr541GV.allTM.MD1,Olfr541GV.allTM.MD2,Olfr541GV.allTM.MD3,Olfr541GV.allTM.MD4,Olfr541GV.allTM.MD5,Olfr541GV.allTM.MD6)

plot(1,xlim=c(0,28), ylim=c(0,max(m)+1),axes = FALSE, type="n", xlab=NA, ylab="RMSD")
axis(1, labels = NA, xaxt="n"); axis(2, las=1)
boxplot(Olfr541.allTM.MD1[950:1050], at = 1, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr541.allTM.MD2[950:1050], at = 2, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr541.allTM.MD3[950:1050], at = 3, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr541.allTM.MD4[950:1050], at = 4, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr541.allTM.MD5[950:1050], at = 5, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr541.allTM.MD6[950:1050], at = 6, add = TRUE, yaxt  = "n", width = 1.5)

boxplot(Olfr541C154G.allTM.MD1[950:1050], at = 8, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr541C154G.allTM.MD2[950:1050], at = 9, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr541C154G.allTM.MD3[950:1050], at = 10, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr541C154G.allTM.MD4[950:1050], at = 11, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr541C154G.allTM.MD5[950:1050], at = 12, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr541C154G.allTM.MD6[950:1050], at = 13, add = TRUE, yaxt  = "n", width = 1.5)

boxplot(Olfr541GV.allTM.MD1[950:1050], at = 15, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr541GV.allTM.MD2[950:1050], at = 16, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr541GV.allTM.MD3[950:1050], at = 17, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr541GV.allTM.MD4[950:1050], at = 18, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr541GV.allTM.MD5[950:1050], at = 19, add = TRUE, yaxt  = "n", width = 1.5)
boxplot(Olfr541GV.allTM.MD6[950:1050], at = 20, add = TRUE, yaxt  = "n", width = 1.5)

segments(1,0,6,0, lwd = 2);segments(8,0,13,0, lwd = 2);segments(15,0,20,0, lwd = 2);
text(3.5,-1,"Olfr541");text(10.5,-1,"C154G");text(17.5,-1,"C154G/G209V");

Olfr539mutants <- read.csv("Olfr539mutants.csv")[c(17,6,2,1),]
Olfr541mutants <- read.csv("Olfr541mutants.csv")[c(5,1,3),]

mrmsd1 <- c(var(c(mean(Olfr539.allTM.MD1[950:1050]),
                  mean(Olfr539.allTM.MD2[950:1050]),
                  mean(Olfr539.allTM.MD3[950:1050]),
                  mean(Olfr539.allTM.MD4[950:1050]),
                  mean(Olfr539.allTM.MD5[950:1050]),
                  mean(Olfr539.allTM.MD6[950:1050]))),
            var(c(mean(Olfr539L155A.allTM.MD1[950:1050]),
                  mean(Olfr539L155A.allTM.MD2[950:1050]),
                  mean(Olfr539L155A.allTM.MD3[950:1050]),
                  mean(Olfr539L155A.allTM.MD4[950:1050]),
                  mean(Olfr539L155A.allTM.MD5[950:1050]),
                  mean(Olfr539L155A.allTM.MD6[950:1050]))),
            var(c(mean(Olfr539V209G.allTM.MD1[950:1050]),
                  mean(Olfr539V209G.allTM.MD2[950:1050]),
                  mean(Olfr539V209G.allTM.MD3[950:1050]),
                  mean(Olfr539V209G.allTM.MD4[950:1050]),
                  mean(Olfr539V209G.allTM.MD5[950:1050]),
                  mean(Olfr539V209G.allTM.MD6[950:1050]))),
            var(c(mean(Olfr539G154C.allTM.MD1[950:1050]),
                  mean(Olfr539G154C.allTM.MD2[950:1050]),
                  mean(Olfr539G154C.allTM.MD3[950:1050]),
                  mean(Olfr539G154C.allTM.MD4[950:1050]),
                  mean(Olfr539G154C.allTM.MD5[950:1050]),
                  mean(Olfr539G154C.allTM.MD6[950:1050]))))

mrmsd2 <- c(var(c(mean(Olfr541.allTM.MD1[950:1050]),
            mean(Olfr541.allTM.MD2[950:1050]),
            mean(Olfr541.allTM.MD3[950:1050]),
            mean(Olfr541.allTM.MD4[950:1050]),
            mean(Olfr541.allTM.MD5[950:1050]),
            mean(Olfr541.allTM.MD6[950:1050]))),
        var(c(mean(Olfr541C154G.allTM.MD1[950:1050]),
            mean(Olfr541C154G.allTM.MD2[950:1050]),
            mean(Olfr541C154G.allTM.MD3[950:1050]),
            mean(Olfr541C154G.allTM.MD4[950:1050]),
            mean(Olfr541C154G.allTM.MD5[950:1050]),
            mean(Olfr541C154G.allTM.MD6[950:1050]))),
        var(c(mean(Olfr541GV.allTM.MD1[950:1050]),
            mean(Olfr541GV.allTM.MD2[950:1050]),
            mean(Olfr541GV.allTM.MD3[950:1050]),
            mean(Olfr541GV.allTM.MD4[950:1050]),
            mean(Olfr541GV.allTM.MD5[950:1050]),
            mean(Olfr541GV.allTM.MD6[950:1050]))))

pdf("figures/20180330.mrmsd.pdf")
par(mar=c(4, 5, 4, 5), xpd = TRUE, bty="n", mgp = c(4,0.8,0))
plot(1,xlim=c(0.8,4.2), ylim=c(0,max(mrmsd1)+0.2),axes = FALSE, type="n", xlab=NA, ylab="variance of RMSD")
axis(1, labels = NA, xaxt="n"); axis(2, las=1, lwd = 3, cex.axis=2)
points(1:4, mrmsd1, type="p", pch=20, col="red", cex=2)
points(1:4, mrmsd1, type="l", lwd=4, col="red")
#legend("top", c("top half","surface expression"), pch=21, pt.bg = c("red","blue"), col=c("gray"), bty="n" )
text(1,-0.5,"Olfr539", cex=0.8);text(2,-0.5,"L155A", cex=0.8);text(3,-0.5,"V209G", cex=0.8);text(4,-0.5,"G154C", cex=0.8);
par(new = TRUE)
plot(1,xlim=c(0.8,4.2), ylim=c(0,1),axes = FALSE, type="n", xlab=NA, ylab=NA)
axis(1, labels = NA, xaxt="n"); axis(4, las=1, col="black", lwd = 3, cex.axis=2)
points(1:4, Olfr539mutants[,9]/100, type="p", pch=20, col ="blue", cex=2)
points(1:4, Olfr539mutants[,9]/100, type="l", lwd=4, col ="blue")

plot(1,xlim=c(0.8,3.2), ylim=c(0,max(mrmsd2)+0.2),axes = FALSE, type="n", xlab=NA, ylab="variance of RMSD")
axis(1, labels = NA, xaxt="n"); axis(2, las=1, lwd = 3, cex.axis=2)
points(1:3, mrmsd2, type="p", pch=20, col="red", cex=2)
points(1:3, mrmsd2, type="l", lwd=4, col="red")
#legend("top", c("top half","surface expression"), pch=21, pt.bg = c("red","blue"), col=c("gray"), bty="n" )
text(1,-0.5,"Olfr541", cex=0.8);text(2,-0.5,"C154G", cex=0.8);text(3,-0.5,"C154G/G209V", cex=0.8)
par(new = TRUE)
plot(1,xlim=c(0.8,3.2), ylim=c(0,1),axes = FALSE, type="n", xlab=NA, ylab=NA)
axis(1, labels = NA, xaxt="n"); axis(4, las=1, col="black", lwd = 3, cex.axis=2)
points(1:3, Olfr541mutants[,9]/100, type="p", pch=20, col ="blue", cex=2)
points(1:3, Olfr541mutants[,9]/100, type="l", lwd=4, col ="blue")
dev.off()
