##############################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GSVA")
BiocManager::install("GSEABase")
BiocManager::install("RTCGA")
BiocManager::install("AnnotationHub")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("edgeR")
BiocManager::install("GO.db")

library(GSVA)
library(GSEABase)
library(RTCGA)
library(dplyr)
library(AnnotationHub)
library(org.Hs.eg.db)
library(edgeR)
library(GO.db)
library(readr)

##############################################################
checkTCGA("DataSets","STAD")
cancerType <- "STAD"

#################################illumina RNAseqV2
dir.create("STADRNAseqV2")

# Download UCEC mRNA data and store it in folder
downloadTCGA(cancerType, dataSet = "STAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3",destDir="STADRNAseqV2",
date = NULL, untarFile = TRUE, removeTar = TRUE, allDataSets = FALSE)

# shortening paths and directories
list.files("STADRNAseqV2/") %>%
  file.path("STADRNAseqV2", .) %>%
  file.rename( to = substr(.,start=1,stop=50))

# reading data
list.files("STADRNAseqV2/") %>%
file.path("STADRNAseqV2", .) -> folder
folder %>%
list.files %>%
file.path( folder, .) %>%
grep( pattern = 'illuminahiseq', x = ., value = TRUE) -> pathRNA
readTCGA( path = pathRNA, dataType = 'rnaseq' ) -> my_data

######################################################
map <- select(org.Hs.eg.db,keys(org.Hs.eg.db),"SYMBOL","ENTREZID")
write.csv(map,file="EntrezID_SYMBOL.csv",row.names=F)

row.names(my_data) <- my_data[,1]
STADori <- t(my_data[,-1])

###############################################################ID convert
number <- 1
check <- row.names(STADori)

for(i in 1:nrow(STADori)){
if(sub("\\|.*","",check[i])=="?"){

 if(length(map[map$ENTREZID==sub("^.*\\|","",check[i]),2])==0){
 check[i] <- paste("Unknown",number,sep="")
 number=number+1
 }else{
 check[i] <- map[as.character(map$ENTREZID)%in%sub("^.*\\|","",check[i]),2]
 }

}else{
check[i] <- sub("\\|.*","",check[i])
}

}
#####################################################GSVA use ssGSEA
row.names(STADori) <- check
STAD <- STADori
STAD <- STAD[!duplicated(rownames(STAD)), ]


url <- "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_UP&fileType=gmt"
temp_file <- tempfile(fileext = ".gmt")
download.file(url, temp_file, mode = "wb")
gset <- readGMT(temp_file)

# Run ssGSEA
gsvaPar <- ssgseaParam(STAD, gset)
TES <- gsva(gsvaPar)

######################################################tSNE
install.packages("Rtsne")
library(Rtsne)

set.seed(19960219)
tsne <- Rtsne(t(TES),pca=FALSE,theta=0.0)

#Create a function to generate a continuous color palette
RWBcolor <- colorRampPalette(c("blue","white","red"))

#This adds a column of color values
TEScolor <- RWBcolor(tsne$N)[as.numeric(cut(colSums(TES),breaks = tsne$N))]



tiff(file="STAD_Tcell Responses Signature.tiff",width=6, height=6, units="in", res=500)
plot(tsne$Y,main=paste("STAD (n=",nrow(tsne$Y),")",sep=""),xlab="tSNE 1",ylab="tSNE 2","cex.main"=2,"cex.lab"=1.5,pch=21,cex=2,col="black",bg=TEScolor,lwd=2)
abline(v=0,lty=3,lwd=2)
abline(a=0,b=0,lty=3,lwd=2)

dev.off()
#########################################################Epigenetic members mutation
Mut <- read.csv(file="I:\\Michael_MOSTProject\\cBioprotal\\STAD_Top50_non-silent mutation and CNV.csv",row.names=1)

tsne_sample <- as.data.frame(tsne$Y)
rownames(tsne_sample) <- sub("01A-.*","01",colnames(TES))

KMT2D <- tsne_sample[rownames(tsne_sample)%in%rownames(Mut[Mut$KMT2D==1,]),]
KDM6A <- tsne_sample[rownames(tsne_sample)%in%rownames(Mut[Mut$KDM6A==1,]),]
ARID1A <- tsne_sample[rownames(tsne_sample)%in%rownames(Mut[Mut$ARID1A==1,]),]
KMT2C <- tsne_sample[rownames(tsne_sample)%in%rownames(Mut[Mut$KMT2C==1,]),]
EP300 <- tsne_sample[rownames(tsne_sample)%in%rownames(Mut[Mut$EP300==1,]),]
CREBBP <- tsne_sample[rownames(tsne_sample)%in%rownames(Mut[Mut$CREBBP==1,]),]
KMT2A <- tsne_sample[rownames(tsne_sample)%in%rownames(Mut[Mut$KMT2A==1,]),]
MIR31HG <- tsne_sample[rownames(tsne_sample)%in%rownames(Mut[Mut$MIR31HG <= -0.5,]),]

tiff(file="STAD_KMT2D.tiff",width=6, height=6, units="in", res=500)
plot(tsne$Y,main="KMT2D non-silent mutation",xlab="tSNE 1",ylab="tSNE 2","cex.main"=2,"cex.lab"=1.5,pch=21,cex=2,col="black",bg=TEScolor,lwd=1)
points(KMT2D$V1,KMT2D$V2,pch=21,cex=2,col="black",lwd=3)
abline(v=0,lty=3,lwd=2)
abline(a=0,b=0,lty=3,lwd=2)
dev.off()

tiff(file="STAD_KDM6A.tiff",width=6, height=6, units="in", res=500)
plot(tsne$Y,main="KDM6A non-silent mutation",xlab="tSNE 1",ylab="tSNE 2","cex.main"=2,"cex.lab"=1.5,pch=21,cex=2,col="black",bg=TEScolor,lwd=1)
points(KDM6A$V1,KDM6A$V2,pch=21,cex=2,col="black",lwd=3)
abline(v=0,lty=3,lwd=2)
abline(a=0,b=0,lty=3,lwd=2)
dev.off()

tiff(file="STAD_ARID1A.tiff",width=6, height=6, units="in", res=500)
plot(tsne$Y,main="ARID1A non-silent mutation",xlab="tSNE 1",ylab="tSNE 2","cex.main"=2,"cex.lab"=1.5,pch=21,cex=2,col="black",bg=TEScolor,lwd=1)
points(ARID1A$V1,ARID1A$V2,pch=21,cex=2,col="black",lwd=3)
abline(v=0,lty=3,lwd=2)
abline(a=0,b=0,lty=3,lwd=2)
dev.off()

tiff(file="STAD_KMT2C.tiff",width=6, height=6, units="in", res=500)
plot(tsne$Y,main="KMT2C non-silent mutation",xlab="tSNE 1",ylab="tSNE 2","cex.main"=2,"cex.lab"=1.5,pch=21,cex=2,col="black",bg=TEScolor,lwd=1)
points(KMT2C$V1,KMT2C$V2,pch=21,cex=2,col="black",lwd=3)
abline(v=0,lty=3,lwd=2)
abline(a=0,b=0,lty=3,lwd=2)
dev.off()

tiff(file="STAD_EP300.tiff",width=6, height=6, units="in", res=500)
plot(tsne$Y,main="EP300 non-silent mutation",xlab="tSNE 1",ylab="tSNE 2","cex.main"=2,"cex.lab"=1.5,pch=21,cex=2,col="black",bg=TEScolor,lwd=1)
points(EP300$V1,EP300$V2,pch=21,cex=2,col="black",lwd=3)
abline(v=0,lty=3,lwd=2)
abline(a=0,b=0,lty=3,lwd=2)
dev.off()

tiff(file="STAD_CREBBP.tiff",width=6, height=6, units="in", res=500)
plot(tsne$Y,main="CREBBP non-silent mutation",xlab="tSNE 1",ylab="tSNE 2","cex.main"=2,"cex.lab"=1.5,pch=21,cex=2,col="black",bg=TEScolor,lwd=1)
points(CREBBP$V1,CREBBP$V2,pch=21,cex=2,col="black",lwd=3)
abline(v=0,lty=3,lwd=2)
abline(a=0,b=0,lty=3,lwd=2)
dev.off()

tiff(file="STAD_KMT2A.tiff",width=6, height=6, units="in", res=500)
plot(tsne$Y,main="KMT2A non-silent mutation",xlab="tSNE 1",ylab="tSNE 2","cex.main"=2,"cex.lab"=1.5,pch=21,cex=2,col="black",bg=TEScolor,lwd=1)
points(KMT2A$V1,KMT2A$V2,pch=21,cex=2,col="black",lwd=3)
abline(v=0,lty=3,lwd=2)
abline(a=0,b=0,lty=3,lwd=2)
dev.off()

#########################################################Targeted Gene
TEScolor <- RWBcolor(tsne$N)[as.numeric(cut(colSums(TES),breaks = tsne$N))]

tsne_sample <- as.data.frame(tsne$Y)
rownames(tsne_sample) <- colnames(TES)
selectdata <- as.data.frame(t(STAD))

TargetGene <- tsne_sample[rownames(tsne_sample)%in%rownames(selectdata[selectdata$EHMT2>=quantile(selectdata$EHMT2,0.75),]),]



tiff(file="STAD_Tcell Responses_G9a(EHMT2).tiff",width=6, height=6, units="in", res=500)
plot(tsne$Y,main=paste("STAD (n=",nrow(tsne$Y),")",sep=""),xlab="tSNE 1",ylab="tSNE 2","cex.main"=2,"cex.lab"=1.5,pch=21,cex=2,col="black",bg=TEScolor,lwd=2)

points(TargetGene$V1,TargetGene$V2,pch=21,cex=2,col="black",lwd=3)

abline(v=0,lty=3,lwd=2)
abline(a=0,b=0,lty=3,lwd=2)

dev.off()
#########################################################DEG
AllES <- as.data.frame(colSums(TES))
colnames(AllES) <- "TES"
AllES$ID <- rownames(AllES)

LoES <- subset(AllES,AllES$TES <= quantile(AllES$TES, 0.25))
HiES <- subset(AllES,AllES$TES >= quantile(AllES$TES, 0.75))

LoES_STAD <- STADori[,colnames(STADori)%in%LoES$ID]
HiES_STAD <- STADori[,colnames(STADori)%in%HiES$ID]

DEGinput <- cbind(LoES_STAD,HiES_STAD)
DEGgroup <- factor(c(rep("LoES",107),rep("HiES",107)))
DEGlist <- DGEList(counts=DEGinput,genes=rownames(DEGinput),group=DEGgroup)

egSYMBOL <- toTable(org.Hs.egSYMBOL)
m <- match(sub("^.*\\|","",DEGlist$genes$genes), egSYMBOL$gene_id)
DEGlist$genes$Symbol <- egSYMBOL$symbol[m]
DEGlist$genes$genes <- sub("^.*\\|","",DEGlist$genes$genes)

DEGlist <- calcNormFactors(DEGlist)
design <- model.matrix(~0+DEGgroup)
DEGlist <- estimateDisp(DEGlist,design)

fit <- glmFit(DEGlist,design)
tr <- glmTreat(fit,coef=2,contrast=c(1,-1),lfc=1)
topTags(tr)
summary(decideTests(tr))

out <- topTags(tr, n=Inf, adjust.method="BH")
keep <- out$table$FDR <= 0.05 & abs(out$table$logFC) >= 1
selectOut <- out[keep,]
DEGresult <- as.data.frame(selectOut)

write.csv(DEGresult,file="DEG_STAD_Hot-Cold.csv",row.names=F)

#######plot DEG select
Out <- as.data.frame(out)

tiff(file="STAD_DEGselect.tiff",width=10, height=10, units="in", res=500)

plot(Out$logCPM,Out$logFC,main="Differential Expressed Genes from Hot & Cold tumors",xlab="Log-transformed Expression",ylab="Log-transformed Fold Change","cex.main"=2,"cex.lab"=1.5,pch=1,cex=1,col="black")
points(DEGresult[DEGresult$logFC>0,5],DEGresult[DEGresult$logFC>0,3],pch=16,cex=1,col="red")
points(DEGresult[DEGresult$logFC<0,5],DEGresult[DEGresult$logFC<0,3],pch=16,cex=1,col="blue")
abline(h=c(-1, 1), col="black", lty=2, lwd=2)
dev.off()


#############################
go <- goana(tr, geneid=tr$genes$genes )

HotPath <- topGO(go, ont="BP", sort="Up", n=50)
write.csv(HotPath,file="HotPath_STAD_Top50.csv",row.names=F)
ColdPath <- topGO(go, ont="BP", sort="Down", n=50)
write.csv(ColdPath,file="ColdPath_STAD_Top50.csv",row.names=F)

HotPath <- topGO(go, ont="BP", sort="Up", n=20)
write.csv(HotPath,file="HotPath_STAD_Top20.csv",row.names=F)
ColdPath <- topGO(go, ont="BP", sort="Down", n=20)
write.csv(ColdPath,file="ColdPath_STAD_Top20.csv",row.names=F)


#########################Hit count Top20
egSYMBOL <- toTable(org.Hs.egSYMBOL)
egGENENAME <- toTable(org.Hs.egGENENAME)


HotCount <- c(rep(0,nrow(DEGresult)))
for(i in 1:nrow(HotPath)){
x <- org.Hs.egGO2ALLEGS
Rkeys(x) <- rownames(HotPath)[i]
EG <- mappedLkeys(x)
HotCount <- HotCount+as.numeric(DEGresult$genes%in%EG)
}
DEGresult$HotHitCounts <- HotCount

ColdCount <- c(rep(0,nrow(DEGresult)))
for(i in 1:nrow(ColdPath)){
x <- org.Hs.egGO2ALLEGS
Rkeys(x) <- rownames(ColdPath)[i]
EG <- mappedLkeys(x)
ColdCount <- ColdCount+as.numeric(DEGresult$genes%in%EG)
}
DEGresult$ColdHitCounts <- ColdCount

m <- match(DEGresult[,1], egGENENAME$gene_id)
DEGresult$GeneName <- egGENENAME$gene_name[m]

write.csv(DEGresult,file="DEG_STAD_Hot-Cold_new_GO.Hit.Count.csv",row.names=F)


PathFig <- HotPath
PathFig <- PathFig[order(PathFig$P.Up,decreasing=TRUE),]
PathFig$Value <- abs(log(PathFig$P.Up,10))

tiff(file="STAD_HotPath.tiff",width=10, height=10, units="in", res=500)
par(mar = c(5, 20, 5, 10))
barplot(PathFig$Value,horiz=T,main="GO Biological Process:Top20 enrich pathways in higher T-cell responses patients",xlim=c(0,250),xlab="-Log10(p.value)",cex.main=1,cex.lab=1,names.arg=PathFig$Term,col="red",las=2)
dev.off()

PathFig <- ColdPath
PathFig <- PathFig[order(PathFig$P.Down,decreasing=TRUE),]
PathFig$Value <- abs(log(PathFig$P.Down,10))

tiff(file="STAD_ColdPath.tiff",width=10, height=10, units="in", res=500)
par(mar = c(5, 20, 5, 10))
barplot(PathFig$Value,horiz=T,main="GO Biological Process:Top20 enrich pathways in lower T-cell responses patients",xlim=c(0,20),xlab="Log10(p.value)",cex.main=1,cex.lab=1,names.arg=PathFig$Term,col="blue",las=2)
dev.off()


#####################################################PCA
pca <- prcomp(t(TES),scale=T)
plot(pca$x[,1],pca$x[,2],main="STAD",xlab="PC1",ylab="PC2","cex.main"=2,"cex.lab"=1.5,,pch=21,cex=2,col="black",bg=TEScolor,lwd=2)







