################# Fst scan ######################
setwd("/home/aa/2alpine/fstScan/")
library(data.table)
library(stringr)
library(dplyr)
library(SuperExactTest)
treshold = 0.99
tresholdGene = 0.1
genes<-fread("/home/aa/Desktop/references/lyrataV2/genesIDsLyV2.gff",h=T)
genes$ID<-substr(genes$ID,4,12)
pdf(paste("WS1/distr_SNPs",treshold,"_Genes",tresholdGene,".pdf",sep=''),width = 9,height = 10)

pops<-c("SPNLOM","NKMSUN","HIAHIF","HFAHFF","HCADRG","OBIGUN","ZEPSUB")

for (p in pops) { #  p = "ZEPSUB"
print(p)
all<-fread(paste("WS1/",p,"_WS1_MS1_BPM.txt",sep = ""),h=T)

### 1. identify 1% outlier SNPs
#all<-all[ order(all[,13],decreasing = T), ]
outl<-subset(all,all$FstH >= quantile(all$FstH,probs = treshold,na.rm = T))
write.table(x=outl,file = paste('WS1/outSNPs_',p,"_",treshold,'.txt',sep=''),append = F,quote = F,col.names = T,row.names = F)

### 2. annotate to genes 
setkey(genes, scaff, start, end)
annot<-foverlaps(x = outl, y = genes, type="within")

### 3. distr of N outlier SNPs/gene
nOutl<-as.data.frame(table(annot$ID))
#  hist(nOutl$Freq,breaks = 100)
#  summary(nOutl$Freq)
og<-subset(genes,genes$ID %in% nOutl$Var1)
og$nsnps<-nOutl$Freq
og$length<-og$end-og$start
#  hist(og$length,breaks = 100)
#  summary(og$length)
og$density<-as.numeric(og$nsnps)/as.numeric(og$length)
# hist(og$density,breaks = 100)
# summary(og$density)
og<-data.frame(og)

### 4. outliers from there - maybe even 10%?
og<-og[ order(og[,12],decreasing = T), ]
outlG<-og[1:(nrow(og)*tresholdGene),]
outlG<-subset(outlG,nsnps>2) 
print(nrow(outlG))
print(summary(outlG[,10:12]))
par(mfrow=c(2,1))
hist(og$length,breaks = 100,xlim = c(0,max(og$length)),main = p)
hist(outlG$length,breaks = 100,xlim = c(0,max(og$length)))
par(mfrow=c(1,1))
hist(outlG$density,breaks = 100,main = p)
write.table(x=outlG,file = paste('WS1/outGenes_',p,"_SNP",treshold,"_Gene",tresholdGene,'.txt',sep=''),append = F,quote = F,col.names = T,row.names = F)
}
dev.off()

### 5. overlap to see what's sensible
file.remove(paste('outGenes_SNP',treshold,"_Gene",tresholdGene,'.txt',sep=''))
for (p in pops) { #  p = "ZEPSUB"
  o<-read.table(paste('WS1/outGenes_',p,"_SNP",treshold,"_Gene",tresholdGene,'.txt',sep=''),h=T)
  write.table(t(c(p,as.character(o$ID))),paste('outGenes_SNP',treshold,"_Gene",tresholdGene,'.txt',sep=''),quote = F,col.names = F,row.names = F,sep = " ",append = T)
}


### 6. SuperExactTest
d <- scan(paste('outGenes_SNP',treshold,"_Gene",tresholdGene,'.txt',sep=''), what="", sep="\n")
d <- strsplit(d, "[[:space:]]+")
names(d) <- c("SPNLOM","NKMSUN","HIAHIF","HFAHFF","HCADRG","OBIGUN","ZEPSUB")
d <- lapply(d, `[`, -1)
total=34051
res=supertest(d, n=total)
res$overlap.sizes
pdf(paste('outGenes_SNP',treshold,"_Gene",tresholdGene,'.pdf',sep=''),width = 11,height = 11,pointsize = 24)
#png(paste(dat,".png",sep=""),width = 1200,height = 960,pointsize = 24)
plot(res, Layout="landscape", sort.by="size", margin=c(0.5,5,1,2),keep.empty.intersections=F,min.intersection.size = 1, degree=2:6,show.overlap.size = F,color.on="black")
dev.off()
write.csv(summary(res)$Table, file=paste('outGenes_SNP',treshold,"_Gene",tresholdGene,'.csv',sep=''), row.names=FALSE)

  
########### FUNCTIONAL FOLLOW-UP 
#INDEPENDENT FROM THE ABOVE - ASSUMES ONE HAVE A GOOD CANDIDATE LIST

setwd("/home/aa/2alpine/fstScan/")
library(data.table)

### 1. Candidate gene annotation TAIR

dict<-fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
s<-read.csv("outGenes_SNP0.99_Gene0.1.csv",h=T)
ann<-fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR10orth_des_20150927.txt",h=T,quote="")
ann2018<-fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20181231.txt",h=T,quote="")
s11<-droplevels(subset(s,s$Degree ==1))
s<-subset(s,s$Degree ==2)
s1<-toString(s$Elements)
ss<-unlist(strsplit(s1, ", "))
parCand<-subset(dict,dict$AL %in% ss) 
parCand<-parCand[!duplicated(parCand[,c('AL')]),] 
sel<-substr(parCand$AT,1,9)
sumPar<-matrix(nrow = nrow(parCand),ncol = 8,dimnames = list(c(),c("AL","AT","lineages","T15_name","T15_descr","T18_descrShort","T18_descrCurat","T18_descrComp")))
for (i in 1:nrow(parCand)) { # i=1
  g<-unlist(parCand[i,1])
  ann1<-subset(ann, ann$`Version-2` %in% g)
  ann2018_1<-subset(ann2018, ann2018$AL %in% g)
pair<-droplevels(subset(s11,s11$Elements %like% g))
pairc<-paste(pair$Intersections, sep="", collapse=', ')
sumPar[i,1]<-unlist(parCand[i,1])
sumPar[i,2]<-unlist(parCand[i,2])
sumPar[i,3]<-pairc
sumPar[i,4]<-unlist(ann1[1,7])
sumPar[i,5]<-unlist(ann1[1,4])
sumPar[i,6]<-unlist(ann2018_1[1,4])
sumPar[i,7]<-unlist(ann2018_1[1,5])
sumPar[i,8]<-unlist(ann2018_1[1,6])
}

write.table(x=sumPar,file = paste('outGenes_SNP0.99_Gene0.1.ALAT.ann.txt',sep=''),append = F,quote = F,col.names = T,row.names = F,sep='\t')
write.table(x=sumPar[,1:3],file = paste('outGenes_SNP0.99_Gene0.1.ALAT.txt',sep=''),append = F,quote = F,col.names = T,row.names = F,sep='\t')



### 2. topGO
library("biomaRt")
library(topGO)

mart <- biomaRt::useMart(biomart = "plants_mart",dataset = "athaliana_eg_gene", host = 'plants.ensembl.org')
GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id", "go_id"), mart = mart)
GTOGO <- GTOGO[GTOGO$go_id != '',]
geneID2GO <- by(GTOGO$go_id,GTOGO$ensembl_gene_id,function(x) as.character(x))
all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
int.genes <- factor(as.integer(all.genes %in% sel))
names(int.genes) = all.genes
go.obj <- new("topGOdata", ontology='BP', allGenes = int.genes, annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize=10) ## 
resultsFe <- runTest(go.obj, algorithm = "elim", statistic = "fisher") #More conservative
resultsFc <- runTest(go.obj, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(go.obj, classicFisher = resultsFc, elimFisher = resultsFe, orderBy = "elimFisher", ranksOf = "elimFisher", topNodes = 500)

allRes$genes <- sapply(allRes$GO.ID, function(x)
{genes<-genesInTerm(go.obj, x) 
genes[[1]][genes[[1]] %in% sel]})

a<-as.data.frame(subset(allRes, allRes$Annotated>=10 & allRes$Annotated <= 500 & allRes$elimFisher <= 0.05 & allRes$Significant > 1))

a$Genes<-""
for (i in 1:length(a$genes)) {
  a$Genes[i]<-paste(unlist(a$genes[i]), sep="", collapse=', ')}
a<-a[,-8]

write.table(a,file=paste("topgo_outGenes_SNP0.99_Gene0.1_min10max500.txt",sep=""),sep="\t",row.names=F)


### 3. visualize - dotplots
setwd("/home/aa/2alpine/fstScan/")
library(data.table)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(scales)
#1. nacist vsechny pary, gff
SPNLOM<-fread("WS1/SPNLOM_WS1_MS1_BPM.txt")
NKMSUN<-fread("WS1/NKMSUN_WS1_MS1_BPM.txt")
HIAHIF<-fread("WS1/HIAHIF_WS1_MS1_BPM.txt")
HFAHFF<-fread("WS1/HFAHFF_WS1_MS1_BPM.txt")
HCADRG<-fread("WS1/HCADRG_WS1_MS1_BPM.txt")
OBIGUN<-fread("WS1/OBIGUN_WS1_MS1_BPM.txt")
ZEPSUB<-fread("WS1/ZEPSUB_WS1_MS1_BPM.txt")
genes<-fread("/home/aa/Desktop/references/lyrataV2/genesIDsLyV2.gff",h=T)
genes$ID<-as.character(substr(genes$ID,4,12))
cand<-read.table("outGenes_SNP0.99_Gene0.1.ALAT.txt",sep="\t",h=T)
flank<-7500

pdf(paste("Figures/all_",flank,".pdf",sep=''),width = 6,height = 20,pointsize = 17)

for (i in 1:nrow(cand)) { # i=1
gen<-as.character(droplevels(cand[i,1]))
int<-subset(genes,genes$ID %in% gen)
first<-int$start-flank
last<-int$end+flank
SPNLOMgene<-subset(SPNLOM, SPNLOM$end > first & SPNLOM$end<last)
NKMSUNgene<-subset(NKMSUN, NKMSUN$end > first & NKMSUN$end<last)
HIAHIFgene<-subset(HIAHIF, HIAHIF$end > first & HIAHIF$end<last)
HFAHFFgene<-subset(HFAHFF, HFAHFF$end > first & HFAHFF$end<last)
HCADRGgene<-subset(HCADRG, HCADRG$end > first & HCADRG$end<last)
OBIGUNgene<-subset(OBIGUN, OBIGUN$end > first & OBIGUN$end<last)
ZEPSUBgene<-subset(ZEPSUB, ZEPSUB$end > first & ZEPSUB$end<last)
minpos<-min(c(SPNLOMgene$end,NKMSUNgene$end,HIAHIFgene$end,HFAHFFgene$end,HCADRGgene$end,OBIGUNgene$end,ZEPSUBgene$end))
maxpos<-max(c(SPNLOMgene$end,NKMSUNgene$end,HIAHIFgene$end,HFAHFFgene$end,HCADRGgene$end,OBIGUNgene$end,ZEPSUBgene$end))
maxFst<-max(c(SPNLOMgene$FstH,NKMSUNgene$FstH,HIAHIFgene$FstH,HFAHFFgene$FstH,HCADRGgene$FstH,OBIGUNgene$FstH,ZEPSUBgene$FstH))

#pdf(paste("Figures/",gen,"_",flank,".pdf",sep=''),width = 6,height = 20,pointsize = 17)
par(mfrow=c(7,1))
par(mar=c(0.1,2.7,0.1,1), mgp=c(1, 0.5, 0))
par(bg=NA) 
plot(ZEPSUBgene$FstH~ZEPSUBgene$end,ylab = "",col = alpha("grey25", 0.8),ylim=c(0.01,maxFst+0.05),xlab="",pch=19,cex=1,xlim=c(minpos,maxpos),xaxt='n', ann=FALSE)
title(ylab="Fst ZEPSUB", line=1.4, cex.lab=1.2)
if (int$orientation %in% "+") {
  arrows(x0 = int$start,y0 = maxFst+0.04,x1 = int$end,y1 =maxFst+0.04,col = "red",lwd = 4,code = 2)
} else {arrows(x0 = int$end,y0 = maxFst+0.04,x1 = int$start,y1 =maxFst+0.04,col = "red",lwd =4,code=2)} 
plot(OBIGUNgene$FstH~OBIGUNgene$end,ylab = "",col = alpha("grey25", 0.8),ylim=c(0.01,maxFst+0.05),xlab="",pch=19,cex=1,xlim=c(minpos,maxpos),xaxt='n', ann=FALSE)
title(ylab="Fst OBIGUN", line=1.4, cex.lab=1.2)
segments(x0 = int$start,y0 = maxFst+0.04,x1 = int$end,y1 =maxFst+0.04,col = "grey50",lwd = 4)
plot(HCADRGgene$FstH~HCADRGgene$end,ylab = "",col = alpha("grey25", 0.8),ylim=c(0.01,maxFst+0.05),xlab="",pch=19,cex=1,xlim=c(minpos,maxpos),xaxt='n', ann=FALSE)
title(ylab="Fst HCADRG", line=1.4, cex.lab=1.2)
segments(x0 = int$start,y0 = maxFst+0.04,x1 = int$end,y1 =maxFst+0.04,col = "grey50",lwd = 4)
plot(HFAHFFgene$FstH~HFAHFFgene$end,ylab = "",col = alpha("grey25", 0.8),ylim=c(0.01,maxFst+0.05),xlab="",pch=19,cex=1,xlim=c(minpos,maxpos),xaxt='n', ann=FALSE)
title(ylab="Fst HFAHFF", line=1.4, cex.lab=1.2)
segments(x0 = int$start,y0 = maxFst+0.04,x1 = int$end,y1 =maxFst+0.04,col = "grey50",lwd = 4)
plot(HIAHIFgene$FstH~HIAHIFgene$end,ylab = "",col = alpha("grey25", 0.8),ylim=c(0.01,maxFst+0.05),xlab="",pch=19,cex=1,xlim=c(minpos,maxpos),xaxt='n', ann=FALSE)
title(ylab="Fst HIAHIF", line=1.4, cex.lab=1.2)
segments(x0 = int$start,y0 = maxFst+0.04,x1 = int$end,y1 =maxFst+0.04,col = "grey50",lwd = 4)
plot(NKMSUNgene$FstH~NKMSUNgene$end,ylab = "",col = alpha("grey25", 0.8),ylim=c(0.01,maxFst+0.05),xlab="",pch=19,cex=1,xlim=c(minpos,maxpos),xaxt='n', ann=FALSE)
title(ylab="Fst NKMSUN", line=1.4, cex.lab=1.2)
segments(x0 = int$start,y0 = maxFst+0.04,x1 = int$end,y1 =maxFst+0.04,col = "grey50",lwd = 4)
par(mar=c(2.5,2.7,0.1,1))
plot(SPNLOMgene$FstH~SPNLOMgene$end,ylab = "",col = alpha("grey25", 0.8),ylim=c(0.01,maxFst+0.05),xlab="",pch=19,cex=1,xlim=c(minpos,maxpos))
title(ylab="Fst SPNLOM", line=1.4, cex.lab=1.2)
title(xlab=paste(gen,as.character(droplevels(cand[i,2])),as.character(droplevels(cand[i,3])),sep= " - "), line=1.4, cex.lab=1.2)
segments(x0 = int$start,y0 = maxFst+0.04,x1 = int$end,y1 =maxFst+0.04,col = "grey50",lwd = 4)
#dev.off()
}

dev.off()

### 3. visualize sanity check scaffs
NKMSUNs1<-subset(NKMSUN,scaff%in% "scaffold_1")
NKMSUNs2<-subset(NKMSUN,scaff%in% "scaffold_2")
NKMSUNs3<-subset(NKMSUN,scaff%in% "scaffold_3")
NKMSUNs4<-subset(NKMSUN,scaff%in% "scaffold_4")
NKMSUNs5<-subset(NKMSUN,scaff%in% "scaffold_5")
NKMSUNs6<-subset(NKMSUN,scaff%in% "scaffold_6")
NKMSUNs7<-subset(NKMSUN,scaff%in% "scaffold_7")
NKMSUNs8<-subset(NKMSUN,scaff%in% "scaffold_8")
png('NKMSUN.png',width = 2400,height = 3000)
par(mfrow=c(8,1))
plot(NKMSUNs1$FstH~NKMSUNs1$end)
plot(NKMSUNs2$FstH~NKMSUNs2$end)
plot(NKMSUNs3$FstH~NKMSUNs3$end)
plot(NKMSUNs4$FstH~NKMSUNs4$end)
plot(NKMSUNs5$FstH~NKMSUNs5$end)
plot(NKMSUNs6$FstH~NKMSUNs6$end)
plot(NKMSUNs7$FstH~NKMSUNs7$end)
plot(NKMSUNs8$FstH~NKMSUNs8$end)
dev.off()
  
SPNLOMs1<-subset(SPNLOM,scaff%in% "scaffold_1")
SPNLOMs2<-subset(SPNLOM,scaff%in% "scaffold_2")
SPNLOMs3<-subset(SPNLOM,scaff%in% "scaffold_3")
SPNLOMs4<-subset(SPNLOM,scaff%in% "scaffold_4")
SPNLOMs5<-subset(SPNLOM,scaff%in% "scaffold_5")
SPNLOMs6<-subset(SPNLOM,scaff%in% "scaffold_6")
SPNLOMs7<-subset(SPNLOM,scaff%in% "scaffold_7")
SPNLOMs8<-subset(SPNLOM,scaff%in% "scaffold_8")
png('SPNLOM.png',width = 2400,height = 3000)
par(mfrow=c(8,1))
plot(SPNLOMs1$FstH~SPNLOMs1$end)
plot(SPNLOMs2$FstH~SPNLOMs2$end)
plot(SPNLOMs3$FstH~SPNLOMs3$end)
plot(SPNLOMs4$FstH~SPNLOMs4$end)
plot(SPNLOMs5$FstH~SPNLOMs5$end)
plot(SPNLOMs6$FstH~SPNLOMs6$end)
plot(SPNLOMs7$FstH~SPNLOMs7$end)
plot(SPNLOMs8$FstH~SPNLOMs8$end)
dev.off()

HIAHIFs1<-subset(HIAHIF,scaff%in% "scaffold_1")
HIAHIFs2<-subset(HIAHIF,scaff%in% "scaffold_2")
HIAHIFs3<-subset(HIAHIF,scaff%in% "scaffold_3")
HIAHIFs4<-subset(HIAHIF,scaff%in% "scaffold_4")
HIAHIFs5<-subset(HIAHIF,scaff%in% "scaffold_5")
HIAHIFs6<-subset(HIAHIF,scaff%in% "scaffold_6")
HIAHIFs7<-subset(HIAHIF,scaff%in% "scaffold_7")
HIAHIFs8<-subset(HIAHIF,scaff%in% "scaffold_8")
png('HIAHIF.png',width = 2400,height = 3000)
par(mfrow=c(8,1))
plot(HIAHIFs1$FstH~HIAHIFs1$end)
plot(HIAHIFs2$FstH~HIAHIFs2$end)
plot(HIAHIFs3$FstH~HIAHIFs3$end)
plot(HIAHIFs4$FstH~HIAHIFs4$end)
plot(HIAHIFs5$FstH~HIAHIFs5$end)
plot(HIAHIFs6$FstH~HIAHIFs6$end)
plot(HIAHIFs7$FstH~HIAHIFs7$end)
plot(HIAHIFs8$FstH~HIAHIFs8$end)
dev.off()

HFAHFFs1<-subset(HFAHFF,scaff%in% "scaffold_1")
HFAHFFs2<-subset(HFAHFF,scaff%in% "scaffold_2")
HFAHFFs3<-subset(HFAHFF,scaff%in% "scaffold_3")
HFAHFFs4<-subset(HFAHFF,scaff%in% "scaffold_4")
HFAHFFs5<-subset(HFAHFF,scaff%in% "scaffold_5")
HFAHFFs6<-subset(HFAHFF,scaff%in% "scaffold_6")
HFAHFFs7<-subset(HFAHFF,scaff%in% "scaffold_7")
HFAHFFs8<-subset(HFAHFF,scaff%in% "scaffold_8")
png('HFAHFF.png',width = 2400,height = 3000)
par(mfrow=c(8,1))
plot(HFAHFFs1$FstH~HFAHFFs1$end)
plot(HFAHFFs2$FstH~HFAHFFs2$end)
plot(HFAHFFs3$FstH~HFAHFFs3$end)
plot(HFAHFFs4$FstH~HFAHFFs4$end)
plot(HFAHFFs5$FstH~HFAHFFs5$end)
plot(HFAHFFs6$FstH~HFAHFFs6$end)
plot(HFAHFFs7$FstH~HFAHFFs7$end)
plot(HFAHFFs8$FstH~HFAHFFs8$end)
dev.off()

OBIGUNs1<-subset(OBIGUN,scaff%in% "scaffold_1")
OBIGUNs2<-subset(OBIGUN,scaff%in% "scaffold_2")
OBIGUNs3<-subset(OBIGUN,scaff%in% "scaffold_3")
OBIGUNs4<-subset(OBIGUN,scaff%in% "scaffold_4")
OBIGUNs5<-subset(OBIGUN,scaff%in% "scaffold_5")
OBIGUNs6<-subset(OBIGUN,scaff%in% "scaffold_6")
OBIGUNs7<-subset(OBIGUN,scaff%in% "scaffold_7")
OBIGUNs8<-subset(OBIGUN,scaff%in% "scaffold_8")
png('OBIGUN.png',width = 2400,height = 3000)
par(mfrow=c(8,1))
plot(OBIGUNs1$FstH~OBIGUNs1$end)
plot(OBIGUNs2$FstH~OBIGUNs2$end)
plot(OBIGUNs3$FstH~OBIGUNs3$end)
plot(OBIGUNs4$FstH~OBIGUNs4$end)
plot(OBIGUNs5$FstH~OBIGUNs5$end)
plot(OBIGUNs6$FstH~OBIGUNs6$end)
plot(OBIGUNs7$FstH~OBIGUNs7$end)
plot(OBIGUNs8$FstH~OBIGUNs8$end)
dev.off()

HCADRGs1<-subset(HCADRG,scaff%in% "scaffold_1")
HCADRGs2<-subset(HCADRG,scaff%in% "scaffold_2")
HCADRGs3<-subset(HCADRG,scaff%in% "scaffold_3")
HCADRGs4<-subset(HCADRG,scaff%in% "scaffold_4")
HCADRGs5<-subset(HCADRG,scaff%in% "scaffold_5")
HCADRGs6<-subset(HCADRG,scaff%in% "scaffold_6")
HCADRGs7<-subset(HCADRG,scaff%in% "scaffold_7")
HCADRGs8<-subset(HCADRG,scaff%in% "scaffold_8")
png('HCADRG.png',width = 2400,height = 3000)
par(mfrow=c(8,1))
plot(HCADRGs1$FstH~HCADRGs1$end)
plot(HCADRGs2$FstH~HCADRGs2$end)
plot(HCADRGs3$FstH~HCADRGs3$end)
plot(HCADRGs4$FstH~HCADRGs4$end)
plot(HCADRGs5$FstH~HCADRGs5$end)
plot(HCADRGs6$FstH~HCADRGs6$end)
plot(HCADRGs7$FstH~HCADRGs7$end)
plot(HCADRGs8$FstH~HCADRGs8$end)
dev.off()

ZEPSUBs1<-subset(ZEPSUB,scaff%in% "scaffold_1")
ZEPSUBs2<-subset(ZEPSUB,scaff%in% "scaffold_2")
ZEPSUBs3<-subset(ZEPSUB,scaff%in% "scaffold_3")
ZEPSUBs4<-subset(ZEPSUB,scaff%in% "scaffold_4")
ZEPSUBs5<-subset(ZEPSUB,scaff%in% "scaffold_5")
ZEPSUBs6<-subset(ZEPSUB,scaff%in% "scaffold_6")
ZEPSUBs7<-subset(ZEPSUB,scaff%in% "scaffold_7")
ZEPSUBs8<-subset(ZEPSUB,scaff%in% "scaffold_8")
png('ZEPSUB.png',width = 2400,height = 3000)
par(mfrow=c(8,1))
plot(ZEPSUBs1$FstH~ZEPSUBs1$end)
plot(ZEPSUBs2$FstH~ZEPSUBs2$end)
plot(ZEPSUBs3$FstH~ZEPSUBs3$end)
plot(ZEPSUBs4$FstH~ZEPSUBs4$end)
plot(ZEPSUBs5$FstH~ZEPSUBs5$end)
plot(ZEPSUBs6$FstH~ZEPSUBs6$end)
plot(ZEPSUBs7$FstH~ZEPSUBs7$end)
plot(ZEPSUBs8$FstH~ZEPSUBs8$end)
dev.off()

#4. N associated GOs
# I'm using gene list with manual selection but the simple .ann.txt is also fine
library(data.table)
library(ggpubr)
setwd("/home/aa/2alpine/predictors/N_GOs")
go<-fread("file:///home/aa/Desktop/references/lyrataV2/functions/ATH_GO_GOSLIMOct2020.txt",h=T)
cand<-fread("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.ALAT.ann.123.txt",h=T)
go$ATGO<-paste(go$AT,go$GOcode,sep="_")
go1<-go[ order(go[,16]), ]
go2<-go1[!duplicated(go1[,c('ATGO')]),] 
# All GOs, all genes background
genNgos<-as.data.frame(table(go2$AT))
genNgosSel<-subset(genNgos,genNgos$Var1 %in% substr(cand$AT,1,9))
genNgosNonSel<-subset(genNgos,!genNgos$Var1 %in% substr(cand$AT,1,9))
hist(genNgosSel$Freq,breaks = 50)
hist(genNgosNonSel$Freq,breaks = 50)
summary(genNgosSel)
summary(genNgosNonSel)

sel <- cbind("Parallel",genNgosSel)
non <- cbind("Any other",genNgosNonSel)
colnames(non)<-c('"Parallel"',"Var1","Freq")
all<-as.data.frame(rbind(sel,non))
colnames(all)<-c('Gene',"Var1","Freq")
all$Freq<-as.numeric(as.character(all$Freq))
my_comparisons <- list( c("Parallel", "Any other") )
pdf("GOs_all_Allgenes.pdf",width = 4,height = 4,pointsize = 24)
ggviolin(all, x = 'Gene', y = "Freq", fill = "Gene",
         palette = c("red", "grey50"),
         add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = max(all$Freq)-5)
dev.off()
#all GOs, all selected background
dict<-fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
dat<-read.table("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.txt",sep="\t")
datc<-""
for (i in 1:length(dat$V1)) {
 datc<-paste(datc,unlist(dat$V1[i]), sep='', collapse=NULL)}
datc1<-as.data.frame(strsplit(datc, split=" ")[[1]][-1])
datc1$AllSelected<-substr(datc1$`strsplit(datc, split = " ")[[1]][-1]`,1,9)
setDT(datc1)
datc1[, Count := .N, by=AllSelected]
datc2<-datc1[Count==1]
datc3<-subset(dict,dict$AL %in% datc2$AllSelected)
genNgosNonSel<-subset(genNgos,genNgos$Var1 %in% substr(datc3$AT,1,9))
hist(genNgosSel$Freq,breaks = 50)
hist(genNgosNonSel$Freq,breaks = 50)
summary(genNgosSel)
summary(genNgosNonSel)
sel <- cbind("Parallel",genNgosSel)
non <- cbind("Non-parallel selected",genNgosNonSel)
colnames(non)<-c('"Parallel"',"Var1","Freq")
all<-as.data.frame(rbind(sel,non))
colnames(all)<-c('Gene',"Var1","Freq")
all$Freq<-as.numeric(as.character(all$Freq))
my_comparisons <- list( c("Parallel", "Non-parallel selected") )
pdf("GOs_all_Allselected.pdf",width = 4,height = 4,pointsize = 24)
ggviolin(all, x = 'Gene', y = "Freq", fill = "Gene",
         palette = c("red", "grey50"),
         add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = max(all$Freq)-5)
dev.off()

### By category - the loop does not output PDF :(
cat="C" # "F" "P"
print(cat)
print("All")
go3<-subset(go2,go2$category %in% cat)
#all genes background
genNgos<-as.data.frame(table(go3$AT))
genNgosSel<-subset(genNgos,genNgos$Var1 %in% substr(cand$AT,1,9))
genNgosNonSel<-subset(genNgos,!genNgos$Var1 %in% substr(cand$AT,1,9))
#hist(genNgosSel$Freq,breaks = 50)
#hist(genNgosNonSel$Freq,breaks = 50)
print(summary(genNgosSel))
print(summary(genNgosNonSel))
sel <- cbind("Parallel",genNgosSel)
non <- cbind("Any other",genNgosNonSel)
colnames(non)<-c('"Parallel"',"Var1","Freq")
all<-as.data.frame(rbind(sel,non))
colnames(all)<-c('Gene',"Var1","Freq")
all$Freq<-as.numeric(as.character(all$Freq))
my_comparisons <- list( c("Parallel", "Any other") )
pdf(paste("GOs_",cat,"_Allgenes.pdf",sep=""),width = 4,height = 4,pointsize = 24)
ggviolin(all, x = 'Gene', y = "Freq", fill = "Gene",
         palette = c("red", "grey50"),
         add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = max(all$Freq)-0.5)
dev.off()
#all selected background
print("selected")
dict<-fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
dat<-read.table("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.txt",sep="\t")
datc<-""
for (i in 1:length(dat$V1)) {
  datc<-paste(datc,unlist(dat$V1[i]), sep='', collapse=NULL)}
datc1<-as.data.frame(strsplit(datc, split=" ")[[1]][-1])
datc1$AllSelected<-substr(datc1$`strsplit(datc, split = " ")[[1]][-1]`,1,9)
setDT(datc1)
datc1[, Count := .N, by=AllSelected]
datc2<-datc1[Count==1]
datc3<-subset(dict,dict$AL %in% datc2$AllSelected)
genNgosNonSel<-subset(genNgos,genNgos$Var1 %in% substr(datc3$AT,1,9))
#hist(genNgosSel$Freq,breaks = 50)
#hist(genNgosNonSel$Freq,breaks = 50)
print(summary(genNgosSel))
print(summary(genNgosNonSel))
sel <- cbind("Parallel",genNgosSel)
non <- cbind("Non-parallel selected",genNgosNonSel)
colnames(non)<-c('"Parallel"',"Var1","Freq")
all<-as.data.frame(rbind(sel,non))
colnames(all)<-c('Gene',"Var1","Freq")
all$Freq<-as.numeric(as.character(all$Freq))
my_comparisons <- list( c("Parallel", "Non-parallel selected") )
pdf(paste("GOs_",cat,"_Allselected.pdf",sep=""),width = 4,height = 4,pointsize = 24)
ggviolin(all, x = 'Gene', y = "Freq", fill = "Gene",
         palette = c("red", "grey50"),
         add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = max(all$Freq)-0.5)
dev.off()
  
  
###5. N interactors STRING
###TODO: the effect of annotation quality, what are the hyper-connected genes? Parallel pattern? test separately for low/hyper connected genes



library(data.table)
setwd("/home/aa/2alpine/predictors/STRING/")
s1<-read.table("3702.protein.links.v11.0.txt",h=T,comment.char = "")

score = 500
s<-subset(s1,s1$combined_score >= score)
s$x<-substr(s$protein1 ,6,14)
nInt<-as.data.frame(table(s$x))
hist(nInt$Freq,breaks =100)
summary(nInt)

### Parallel
cand<-fread("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.ALAT.ann.123.txt",h=T)
cand<-fread("file:///home/aa/Desktop/vercaSubmission/ATcodes.txt",h=T)

cand$AT<-substr(cand$AT,1,9)
candInt<-subset(nInt,nInt$Var1 %in% cand$AT)
hist(candInt$Freq,breaks =100)
summary(candInt)



### Selected
#To make a list of all selected genes, without parallel ones
dict<-fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
dat<-read.table("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.txt",sep="\t")
datc<-""
for (i in 1:length(dat$V1)) {
  datc<-paste(datc,unlist(dat$V1[i]), sep='', collapse=NULL)}
datc1<-as.data.frame(strsplit(datc, split=" ")[[1]][-1])
datc1$AllSelected<-substr(datc1$`strsplit(datc, split = " ")[[1]][-1]`,1,9)
setDT(datc1)
datc1[, Count := .N, by=AllSelected]
datc2<-datc1[Count==1]
datc3<-subset(dict,dict$AL %in% datc2$AllSelected)
datc3<-subset(datc3,!datc3$AT %in% "nnn")
datc3$AT<-substr(datc3$AT,1,9)
selInt<-subset(nInt,nInt$Var1 %in% datc3$AT)
hist(selInt$Freq,breaks =100)
summary(selInt)

#parallel-all
print(c(mean(candInt$Freq),mean(nInt$Freq)))
print(c(median(candInt$Freq),median(nInt$Freq)))
print(wilcox.test(candInt$Freq,nInt$Freq))
#parallel-selected
print(c(mean(candInt$Freq),mean(selInt$Freq)))
print(c(median(candInt$Freq),median(selInt$Freq)))
print(wilcox.test(candInt$Freq,selInt$Freq))



