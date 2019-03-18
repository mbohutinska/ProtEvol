################# SELECTION SCANS PIPELINE ######################

################# 1. QUARTET FST ######################
### MAKE A SUMMARY MATRIX - 4 pops ###
setwd("/home/aa/alpine/arenosaGenome/selScans/")
#h1<-"VEL"
#h2<-"ZEP"
#f1<-"SUB"
#f2<-"BAB"
#region<-"Tatry2x"
summary_matrix<-function(h1,h2,f1,f2,region){
  #read in all windows, make a list 
  p1<-fread(paste("BPM/",h1,f1,"_WS1000_MS1_BPM.txt",sep=""),h=T)
  p2<-fread(paste("BPM/",h1,f2,"_WS1000_MS1_BPM.txt",sep=""),h=T)
  p3<-fread(paste("BPM/",h2,f1,"_WS1000_MS1_BPM.txt",sep=""),h=T)
  p4<-fread(paste("BPM/",h2,f2,"_WS1000_MS1_BPM.txt",sep=""),h=T)
  n1<-fread(paste("BPM/",h1,h2,"_WS1000_MS1_BPM.txt",sep=""),h=T)
  n2<-fread(paste("BPM/",f1,f2,"_WS1000_MS1_BPM.txt",sep=""),h=T)
  wf1<-fread(paste("WPM/",f1,".WS1.0k_MS1_6ind_WPM.txt",sep=""),h=T)
  wf2<-fread(paste("WPM/",f2,".WS1.0k_MS1_6ind_WPM.txt",sep=""),h=T)
  wh1<-fread(paste("WPM/",h1,".WS1.0k_MS1_6ind_WPM.txt",sep=""),h=T)
  wh2<-fread(paste("WPM/",h2,".WS1.0k_MS1_6ind_WPM.txt",sep=""),h=T)
  a<-rbind(p1[,2:4],p2[,2:4],p3[,2:4],p4[,2:4],n1[,2:4],n2[,2:4])
  a1<-a[ order(a[,1], a[,2]), ]
  a1<-subset(a1,!a1$scaff %in% "Genome")
  a<-a1[!duplicated(a1[,c('scaff','start')]),] 
  cds<-fread("lyrataCDSs.txt",h=T)
  setkey(p1, scaff, start, end)
  setkey(p2, scaff, start, end)
  setkey(p3, scaff, start, end)
  setkey(p4, scaff, start, end)
  setkey(n1, scaff, start, end)
  setkey(n2, scaff, start, end)
  setkey(wf1, scaff, start, end)
  setkey(wf2, scaff, start, end)
  setkey(wh1, scaff, start, end)
  setkey(wh2, scaff, start, end)
  setkey(a, scaff, start, end)
  mp1<-foverlaps(a, p1, type="start")
  mp2<-foverlaps(a, p2, type="start")
  mp3<-foverlaps(a, p3, type="start")
  mp4<-foverlaps(a, p4, type="start")
  mn1<-foverlaps(a, n1, type="start")
  mn2<-foverlaps(a, n2, type="start")
  mwf1<-foverlaps(a, wf1, type="start",mult="first")
  mwf2<-foverlaps(a, wf2, type="start",mult="first")
  mwh1<-foverlaps(a, wh1, type="start",mult="first")
  mwh2<-foverlaps(a, wh2, type="start",mult="first")
  a$p1_Fst<-mp1$FstWC
  a$p2_Fst<-mp2$FstWC
  a$p3_Fst<-mp3$FstWC
  a$p4_Fst<-mp4$FstWC
  a$n1_Fst<-mn1$FstWC
  a$n2_Fst<-mn2$FstWC
  a$p1_mis<-mp1$num_sites
  a$p2_mis<-mp2$num_sites
  a$p3_mis<-mp3$num_sites
  a$p4_mis<-mp4$num_sites
  a$n1_mis<-mn1$num_sites
  a$n2_mis<-mn2$num_sites
  a$p1_snps<-mp1$num_snps
  a$p2_snps<-mp2$num_snps
  a$p3_snps<-mp3$num_snps
  a$p4_snps<-mp4$num_snps
  a$n1_snps<-mn1$num_snps
  a$n2_snps<-mn2$num_snps
  a$wf1Div<-mwf1$Diversity
  a$wf2Div<-mwf2$Diversity
  a$wh1Div<-mwh1$Diversity
  a$wh2Div<-mwh2$Diversity
  write.table(a,append = F,file = paste("QuartetSummary",region,".txt",sep=""),quote = F, sep = "\t",col.names = T,row.names = F)}
#RUN
#source("summary_matrix.R")
summary_matrix('VEL','ZEP','SUB','BAB','Tatry2x')
summary_matrix('TKO','TRT','HRA','SPI','Tatry4x')
summary_matrix('LAC','BAL','DRA','TIS','Fagaras')
summary_matrix('SCH','WIL','ING','KAS','Alps')

### MAKE A SUMMARY MATRIX - 2 pops ###
setwd("/home/aa/alpine/arenosaGenome/selScans/")
h1<-"INE"
f1<-"CAR"
region<-"Rodna"
h1<-"HCA"
f1<-"DRG"
region<-"FagarasH"
h1<-"OBI"
f1<-"GUN"
region<-"AlpsH"
#read in all windows, make a list 
p1<-read.table(paste("BPM/",h1,f1,"_WS1000_MS1_BPM.txt",sep=""),h=T)
a<-rbind(p1[,c(2,3,6,9)])
a1<-a[ order(a[,1], a[,2]), ]
a1<-subset(a1,!a1$scaff %in% "Genome")
a<-a1[!duplicated(a1[,c('scaff','start')]),] 
m<-matrix(nrow = nrow(a), ncol = 6,dimnames =list(c(),c("scaff","start","p1_Fst","p1_mis",'wf1','wh1')))
m[,1]<-as.character(a$scaff)
m[,2]<-a$start
m[,3]<-a$FstWC
m[,4]<-a$num_sites
### WPM - low/high pi in both/all - TODO - foverlap!!!!!!!!!!!!!!!!!!!!!!!!
wf1<-read.table(paste("WPM/",f1,".WS1.0k_MS1_6ind_WPM.txt",sep=""),h=T)
wh1<-read.table(paste("WPM/",h1,".WS1.0k_MS1_6ind_WPM.txt",sep=""),h=T)
con<-c("wf1","wh1")
for (line in 1:nrow(m)){ # line=1
  for (co in con){ # con=1
    iCo<-which(con %in% paste(co))
    v=get(co)
    sub<-subset(v, v[,4] %in% m[line,1] & v[,5] %in% m[line,2])
    if (dim(sub)[1]==0 | dim(sub)[1]>1) {} else 
    {m[line,iCo+4]<-sub$Diversity}}}
write.table(m,append = F,file = paste("QuartetSummary",region,".txt",sep=""),quote = F, sep = "\t",col.names = T,row.names = F)

### RANK THE WINDOWS - 4 pops ###
setwd("/home/aa/alpine/arenosaGenome/selScans/")
library(data.table)
region<-"Alps"
a<-fread(paste("QuartetSummary",region,".txt",sep=''))
mean(a$p1_snps,na.rm = T)
mean(a$p2_snps,na.rm = T)
mean(a$p3_snps,na.rm = T)
mean(a$p4_snps,na.rm = T)
mean(a$n1_snps,na.rm = T)
mean(a$n2_snps,na.rm = T)
a<-a[order(a$p1_Fst,decreasing = T),]
a$rank_p1<-seq(1,nrow(a),1)
a<-a[order(a$p2_Fst,decreasing = T),]
a$rank_p2<-seq(1,nrow(a),1)
a<-a[order(a$p3_Fst,decreasing = T),]
a$rank_p3<-seq(1,nrow(a),1)
a<-a[order(a$p4_Fst,decreasing = T),]
a$rank_p4<-seq(1,nrow(a),1)
a$sumRank<-a$rank_p1+a$rank_p2+a$rank_p3+a$rank_p4
a<-a[order(a$sumRank,decreasing = F),]
n1<-quantile(a$n1_Fst,0.99,na.rm = T)
n2<-quantile(a$n2_Fst,0.99,na.rm = T)
a$out_n1<-0
a$out_n1[a$n1_Fst > n1] <- 1
a$out_n2<-0
a$out_n2[a$n2_Fst > n2] <- 1
#Add annotation
cds<-fread("lyrataCDSs.txt",h=T)
setkey(cds, scaff, start, end)
aaa<-foverlaps(a, cds, type="any")
aaa$start[is.na(aaa$start)]<-99
aaa$start [aaa$start  < aaa$i.start] <-aaa$i.start [aaa$start  < aaa$i.start]
aaa$end[is.na(aaa$end)]<--99
aaa$end [aaa$end  > aaa$i.end] <-aaa$i.end [aaa$end  > aaa$i.end]
aaa$sum<-aaa$end-(aaa$start)
aaa$sum[is.na(aaa$ID)]<-0
aa<-setDT(aaa)[,list(genic=sum(sum)),by=list(Category=paste(scaff,i.start))]
summary(aa$genic)
hist(aa$genic)
table (aa$genic)
a$genic<-aa$genic
a1<-subset(a,a$out_n1 %in% 0 & a$out_n2 %in% 0)
a2<-subset(a1,a1$p1_mis>200 & a1$p2_mis>200 & a1$p3_mis>200 & a1$p4_mis>200)
o1<-a[1:(nrow(a)*0.01),]
pdf(paste("genicVsFst_",region,".pdf",sep=""),width = 12,height = 8)
plot(a2$genic,a2$p1_Fst,ylab = "Fst",xlab="bp within CDSs")
abline(h = min(o1$p1_Fst))
plot(a2$genic,a2$p1_Fst,ylab = "Fst",xlab="bp within CDSs",xlim=c(0,1000))
abline(h = min(o1$p1_Fst))
dev.off()
write.table(a,paste("rankFst",region,".txt",sep=""),quote=F,row.names = F,sep="\t")
write.table(o1,paste("out_1percent_negOut_moreThan200Sites_",region,".txt",sep=''),quote=F,row.names = F,sep="\t")
#overlap with genes
o2<-foverlaps(o1, cds, type="any")
o3<-unique(o2$ID)
write.table(o3,paste("out_1percent_IDs_",region,".txt",sep=''),quote=F,row.names = F,sep="\t")

### RANK THE WINDOWS - 2 pops ###
require(data.table)
setwd("/home/aa/alpine/arenosaGenome/selScans/")
region<-"FagarasH"
a<-fread(paste("QuartetSummary",region,".txt",sep=""),h=T)
#Add annotation
a$end<-a$start+999
cds<-fread("lyrataCDSs.txt",h=T)
setkey(cds, scaff, start, end)
aaa<-foverlaps(a, cds, type="any")
aaa$start[is.na(aaa$start)]<-99
aaa$start [aaa$start  < aaa$i.start] <-aaa$i.start [aaa$start  < aaa$i.start]
aaa$end[is.na(aaa$end)]<--99
aaa$end [aaa$end  > aaa$i.end] <-aaa$i.end [aaa$end  > aaa$i.end]
aaa$sum<-aaa$end-(aaa$start-1)
aaa$sum[is.na(aaa$ID)]<-0
aa<-setDT(aaa)[,list(genic=sum(sum)),by=list(Category=paste(scaff,i.start))]
summary(aa$genic)
hist(aa$genic)
table (aa$genic)
a$genic<-aa$genic
a<-a[order(a$p1_Fst,decreasing = T),]
a2<-subset(a,a$p1_mis>200)
# a2<-subset(a2,a2$genic>0)
o1<-a2[1:(nrow(a2)*0.01),]
pdf(paste("genicVsFst_",region,".pdf",sep=""),width = 12,height = 8)
plot(a2$genic,a2$p1_Fst,ylab = "Fst",xlab="bp within CDSs")
abline(h = min(o1$p1_Fst))
plot(a2$genic,a2$p1_Fst,ylab = "Fst",xlab="bp within CDSs",xlim=c(0,1000))
abline(h = min(o1$p1_Fst))
dev.off()
write.table(a,paste("rankFst",region,".txt",sep=""),quote=F,row.names = F,sep="\t")
write.table(o1,paste("out_1percent_negOut_moreThan200Sites_",region,".txt",sep=''),quote=F,row.names = F,sep="\t")
#overlap with genes
o2<-foverlaps(o1, cds, type="any")
o3<-unique(o2$ID)
write.table(o3,paste("out_1percent_IDs_",region,".txt",sep=''),quote=F,row.names = F,sep="\t")

### PROCESS CANDIDATES 4 pops ###
#fish://holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/arenosa_merged.filtered/snpEff.sh
#splitVCFsann - .table - I want all?? PopKey Tatry2x etc,..
#sed -i -E 's/\t\t/\tna\t/g' TDI.table.repol.txt
library(data.table)
library(stringr)
library(dplyr)
pops<-c('VEL','ZEP','SUB','BAB')
pops<-c('TKO','TRT','HRA','SPI')
pops<-c('LAC','BAL','DRA','TIS')
region<-"Alps"
sh<-"ALParenosa"
nh1<-7
nh2<-8
nf1<-8
nf2<-8
ploidy<-4
pops<-c('SCH','WIL','KAS','ING')
quant<-0.99
setwd("/home/aa/alpine/arenosaGenome/selScans/")
heatmapCand(c('TKO','TRT','HRA','SPI'),"Tatry4x","TTE",8,8,8,8,4,0.99)
heatmapCand(c('VEL','ZEP','SUB','BAB'),"Tatry2x","TDI",8,8,8,8,2,0.99)
heatmapCand(c('SCH','WIL','KAS','ING'),"Alps","ALParenosa",7,8,8,8,4,0.995)
heatmapCand(c('LAC','BAL','DRA','TIS'),"Fagaras","FAGarenosa",8,8,8,8,4,0.99)
heatmapCand2pop(c('INE','CAR'),"Rodna","ROD",7,8,4,0.99)
heatmapCand2pop(c('OBI','GUN'),"AlpsH","ALPhalleri",8,8,2,0.99)
heatmapCand2pop(c('HCA','DRG'),"FagarasH","FAGhalleri",8,8,2,0.99)

heatmapCand<-function(pops,region,sh,nh1,nh2,nf1,nf2,ploidy,quant){
  library(data.table)
  library(stringr)
  library(dplyr)
  s<-fread(paste('ann/',sh,'.table.recode.txt',sep=''),h=F,na.strings = "-9")
  o<-fread(paste("quartetOutliers/out_1percent_negOut_moreThan200Sites_",region,".txt",sep = ""),h=T)
  o$rank<-seq(1,nrow(o),1)
  o1<-o[ order(o[,1], o[,2]), ]
  #1.extract only outlier sites
  colnames(s) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",seq(10,ncol(s),1))
  s$end<-s$start
  setkey(o1, scaff, start, end)
  aaa<-foverlaps(x = s, y = o1, type="within")
  a<-aaa[complete.cases(aaa[ , 2]),]
  #2. calculate AF for each lineage
  i<-which(colnames(a) %in% "10")
  afh1<-a[,i:as.numeric((i-1)+nh1)]
  a$ACh1<-rowSums(afh1,na.rm = T)
  a$NAh1<-apply(is.na(afh1), 1, sum)
  a$ANh1<-(nh1-a$NAh1)*ploidy
  a$AFh1<-a$ACh1/a$ANh1
  afh2<-a[,as.numeric(i+nh1):as.numeric((i-1)+nh1+nh2)]
  a$ACh2<-rowSums(afh2,na.rm = T)
  a$NAh2<-apply(is.na(afh2), 1, sum)
  a$ANh2<-(nh2-a$NAh2)*ploidy
  a$AFh2<-a$ACh2/a$ANh2
  aff1<-a[,as.numeric(i+nh1+nh2):as.numeric((i-1)+nh1+nh2+nf1)]
  a$ACf1<-rowSums(aff1,na.rm = T)
  a$NAf1<-apply(is.na(aff1), 1, sum)
  a$ANf1<-(nf1-a$NAf1)*ploidy
  a$AFf1<-a$ACf1/a$ANf1
  aff2<-a[,as.numeric(i+nh1+nh2+nf1):as.numeric((i-1)+nh1+nh2+nf1+nf2)]
  a$ACf2<-rowSums(aff2,na.rm = T)
  a$NAf2<-apply(is.na(aff2), 1, sum)
  a$ANf2<-(nf2-a$NAf2)*ploidy
  a$AFf2<-a$ACf2/a$ANf2
  s4<-select(a,ID,i.start,ann,aas,AFh1,AFh2,AFf1,AFf2)
  s4$p1<-abs(s4$AFh1-s4$AFf1) 
  s4$p2<-abs(s4$AFh1-s4$AFf2) 
  s4$p3<-abs(s4$AFh2-s4$AFf1) 
  s4$p4<-abs(s4$AFh2-s4$AFf2) 
  s4$aftot<-abs((s4$AFh1+s4$AFh2)/2-(s4$AFf1+s4$AFf2)/2)
    #   s5<-subset(s4,s4$p1>0.8 | s4$p2>0.8 | s4$p3>0.8 | s4$p4>0.8)
  q<-quantile(s4$aftot,probs = quant,na.rm = T)
  s5<-subset(s4,s4$aftot >= q)
  write.table(x=s5,file = paste('outSNPsAll_',region,quant,'.txt',sep=''),append = F,quote = F,col.names = T,row.names = F)
     #} #############alternative end
  #4. plot
  library(gplots)
  library(RColorBrewer)
  my_palette <- colorRampPalette(c("khaki1", "green2", "blue3"))(n = 20)
  df<-as.matrix(s5[,5:8],rownames = paste(s5$ID,s5$ann,sep="         "))
  pdf(paste("figures/heatmap",region,quant,".pdf",sep=""),height = nrow(s5)/10,width = 7)
  heatmap.2(x = df,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,2,4),sepcolor= c("black"),sepwidth = c (0.05),trace="none",ColSideColors = c("orange","orange","green","green"),labCol=pops,colCol= c("orange","orange","green","green"),offsetRow = c(0.05),offsetCol= c(0.05),cexRow = c(1),cexCol = c(2),margins = c(5,25),lhei = c(0.5,30),lwid = c(0.05,0.8))
  dev.off()
  } #############alternative end
  #5. "Adaptability"
  n<-nrow(subset(s,s$ann %in% "missense_variant"))
  t00<-n/nrow(s)
  n<-nrow(subset(s4,s4$ann %in% "missense_variant"))
  al00<-n/nrow(s4)
  q<-quantile(s4$aftot,probs = 0.95,na.rm = T)
  s5<-subset(s4,s4$aftot >= q)
  n<-nrow(subset(s5,s5$ann %in% "missense_variant"))
  al95<-n/nrow(s5)
  q<-quantile(s4$aftot,probs = 0.99,na.rm = T)
  s5<-subset(s4,s4$aftot >= q)
  n<-nrow(subset(s5,s5$ann %in% "missense_variant"))
  al99<-n/nrow(s5)
  q<-quantile(s4$aftot,probs = 0.995,na.rm = T)
  s5<-subset(s4,s4$aftot >= q)
  n<-nrow(subset(s5,s5$ann %in% "missense_variant"))
  al995<-n/nrow(s5)
  q<-quantile(s4$aftot,probs = 0.999,na.rm = T)
  s5<-subset(s4,s4$aftot >= q)
  n<-nrow(subset(s5,s5$ann %in% "missense_variant"))
  al999<-n/nrow(s5)
  q<-quantile(s4$aftot,probs = 0.9995,na.rm = T)
  s5<-subset(s4,s4$aftot >= q)
  n<-nrow(subset(s5,s5$ann %in% "missense_variant"))
  al9995<-n/nrow(s5)
  write.table(x=paste(region,t00,al00,al95,al99,al995,al999,al9995,sep="\t"),file = 'alpha.txt',append = T,quote = F,col.names = F,row.names = F)}

### PROCESS CANDIDATES 2 pops ###
#fish://holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/arenosa_merged.filtered/snpEff.sh
#splitVCFsann - .table - I want all?? PopKey Tatry2x etc,..
#sed -i -E 's/\t\t/\tna\t/g' TDI.table.repol.txt
library(data.table)
library(stringr)
library(dplyr)
pops<-c('INE','CAR')
pops<-c('OBI','GUN')
pops<-c('HCA','DRG')
region<-"AlpsH"
sh<-"ALPhalleri"
nh1<-8
nf1<-8
ploidy<-4
quant<-0.995
setwd("/home/aa/alpine/arenosaGenome/selScans/")
heatmapCand2pop(c('INE','CAR'),"Rodna","ROD",7,8,4,0.995)
heatmapCand2pop(c('OBI','GUN'),"AlpsH","ALPhalleri",8,8,2,0.995)
heatmapCand2pop(c('HCA','DRG'),"FagarasH","FAGhalleri",8,8,2,0.995)

heatmapCand2pop<-function(pops,region,sh,nh1,nf1,ploidy,quant){
  library(data.table)
  library(stringr)
  library(dplyr)
  s<-fread(paste('ann/',sh,'.table.recode.txt',sep=''),h=F,na.strings = "-9")
  o<-fread(paste("out_1percent_negOut_moreThan200Sites_",region,".txt",sep = ""),h=T)
  o$rank<-seq(1,nrow(o),1)
  o1<-o[ order(o[,1], o[,2]), ]
  #1.extract only outlier sites
  colnames(s) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",seq(10,ncol(s),1))
  s$end<-s$start
  setkey(o1, scaff, start, end)
  aaa<-foverlaps(x = s, y = o1, type="within")
  a<-aaa[complete.cases(aaa[ , 2]),]
  #2. calculate AF for each lineage
  i<-which(colnames(a) %in% "10")
  afh1<-a[,i:as.numeric((i-1)+nh1)]
  a$ACh1<-rowSums(afh1,na.rm = T)
  a$NAh1<-apply(is.na(afh1), 1, sum)
  a$ANh1<-(nh1-a$NAh1)*ploidy
  a$AFh1<-a$ACh1/a$ANh1
  aff1<-a[,as.numeric(i+nh1):as.numeric((i-1)+nh1+nf1)]
  a$ACf1<-rowSums(aff1,na.rm = T)
  a$NAf1<-apply(is.na(aff1), 1, sum)
  a$ANf1<-(nf1-a$NAf1)*ploidy
  a$AFf1<-a$ACf1/a$ANf1
  s4<-select(a,ID,i.start,ann,aas,AFh1,AFf1)
  s4$aftot<-abs(s4$AFh1-s4$AFf1) 
  q<-quantile(s4$aftot,probs = quant,na.rm = T)
  s5<-subset(s4,s4$aftot >= q)
  write.table(x=s5,file = paste('outSNPsAll_',region,'.txt',sep=''),append = F,quote = F,col.names = T,row.names = F)
#} #############alternative end
  #4. plot
  library(gplots)
  library(RColorBrewer)
  my_palette <- colorRampPalette(c("khaki1", "green2", "blue3"))(n = 20)
  df<-as.matrix(s5[,5:6],rownames = paste(s5$ID,s5$ann,sep="         "))
  pdf(paste("figures/heatmap",region,quant,".pdf",sep=""),height = nrow(s5)/10,width = 7)
  heatmap.2(x = df,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,1,2),sepcolor= c("black"),sepwidth = c (0.05),trace="none",ColSideColors = c("orange","green"),labCol=pops,colCol= c("orange","green"),offsetRow = c(0.05),offsetCol= c(0.05),cexRow = c(1),cexCol = c(2),margins = c(5,25),lhei = c(0.5,30),lwid = c(0.05,0.8))
  dev.off()
  #5. "Adaptability"
  n<-nrow(subset(s,s$ann %in% "missense_variant"))
  t00<-n/nrow(s)
  n<-nrow(subset(s4,s4$ann %in% "missense_variant"))
  al00<-n/nrow(s4)
  q<-quantile(s4$aftot,probs = 0.95,na.rm = T)
  s5<-subset(s4,s4$aftot >= q)
  n<-nrow(subset(s5,s5$ann %in% "missense_variant"))
  al95<-n/nrow(s5)
  q<-quantile(s4$aftot,probs = 0.99,na.rm = T)
  s5<-subset(s4,s4$aftot >= q)
  n<-nrow(subset(s5,s5$ann %in% "missense_variant"))
  al99<-n/nrow(s5)
  q<-quantile(s4$aftot,probs = 0.995,na.rm = T)
  s5<-subset(s4,s4$aftot >= q)
  n<-nrow(subset(s5,s5$ann %in% "missense_variant"))
  al995<-n/nrow(s5)
  q<-quantile(s4$aftot,probs = 0.999,na.rm = T)
  s5<-subset(s4,s4$aftot >= q)
  n<-nrow(subset(s5,s5$ann %in% "missense_variant"))
  al999<-n/nrow(s5)
  q<-quantile(s4$aftot,probs = 0.9995,na.rm = T)
  s5<-subset(s4,s4$aftot >= q)
  n<-nrow(subset(s5,s5$ann %in% "missense_variant"))
  al9995<-n/nrow(s5)
  write.table(x=paste(region,t00,al00,al95,al99,al995,al999,al9995,sep="\t"),file = 'alpha.txt',append = T,quote = F,col.names = F,row.names = F)}

################# 2. SNP Fst ######################
setwd("/home/aa/alpine/arenosaGenome/selScans/")
library(data.table)
library(stringr)
library(dplyr)
pops<-c('WIL','KAS')
pops<-c('INE','CAR')
pops<-c('BAL','TIS')
pops<-c('TKO','HRA')
pops<-c('ZEP','SUB')
pops<-c('OBI','GUN')
pops<-c('HCA','DRG')
region<-"Alps"
region<-"Rodna"
region<-"Fagaras"
region<-"Tatry4x"
region<-"Tatry2x"
region<-"AlpsH"
region<-"FagarasH"
sh<-"ALParenosa"
sh<-"ROD"
sh<-"FAGarenosa"
sh<-"TTE"
sh<-"TDI"
sh<-"ALPhalleri"
sh<-"FAGhalleri"
quant<-0.999
s<-fread(paste('ann/',sh,'.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
p=paste(pops,sep = "",collapse = "")
o<-fread(paste("WS1/",p,"_WS1_MS1_BPM.txt",sep = ""),h=T)
#1.extract only outlier sites
colnames(s) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",seq(10,ncol(s),1))
s$end<-s$start
o$start<-o$end
setkey(o, scaff, start, end)
aaa<-foverlaps(x = s, y = o, type="within")
s<-""
o<-""
a<-aaa[complete.cases(aaa[ , 2]),]
aaa<-""
#a1<-a[ order(a[,13],decreasing = T), ]
a1<-subset(a,a$FstH >= quantile(a$FstH,probs = quant,na.rm = T))
write.table(x=a1,file = paste('WS1/outSNPs_',region,quant,'.txt',sep=''),append = F,quote = F,col.names = T,row.names = F)
### "Adaptability" ###
n<-nrow(subset(a,a$ann %in% "missense_variant"))
q00<-n/nrow(a)
q<-quantile(a$FstH,probs = 0.90,na.rm = T)
a1<-subset(a,a$FstH >= q)
n<-nrow(subset(a1,a1$ann %in% "missense_variant"))
q90<-n/nrow(a1)
q<-quantile(a$FstH,probs = 0.95,na.rm = T)
a1<-subset(a,a$FstH >= q)
n<-nrow(subset(a1,a1$ann %in% "missense_variant"))
q95<-n/nrow(a1)
q<-quantile(a$FstH,probs = 0.99,na.rm = T)
a1<-subset(a,a$FstH >= q)
n<-nrow(subset(a1,a1$ann %in% "missense_variant"))
q99<-n/nrow(a1)
q<-quantile(a$FstH,probs = 0.995,na.rm = T)
a1<-subset(a,a$FstH >= q)
n<-nrow(subset(a1,a1$ann %in% "missense_variant"))
q995<-n/nrow(a1)
q<-quantile(a$FstH,probs = 0.999,na.rm = T)
a1<-subset(a,a$FstH >= q)
n<-nrow(subset(a1,a1$ann %in% "missense_variant"))
q999<-n/nrow(a1)
q<-quantile(a$FstH,probs = 0.9995,na.rm = T)
a1<-subset(a,a$FstH >= q)
n<-nrow(subset(a1,a1$ann %in% "missense_variant"))
q9995<-n/nrow(a1)
q<-quantile(a$FstH,probs = 0.9999,na.rm = T)
a1<-subset(a,a$FstH >= q)
n<-nrow(subset(a1,a1$ann %in% "missense_variant"))
q9999<-n/nrow(a1)
write.table(x=paste(region,q00,q90,q95,q99,q995,q999,q9995,q9999,sep="\t"),file = 'WS1/alpha.txt',append = T,quote = F,col.names = F,row.names = F)

#1. find overlaps in arenosa
setwd("/home/aa/alpine/arenosaGenome/selScans/")
library(data.table)
library(stringr)
library(dplyr)
a<-fread("WS1/outSNPs_Alps0.999.txt",h=T)
b<-fread("WS1/outSNPs_Fagaras0.999.txt",h=T)
c<-fread("WS1/outSNPs_Rodna0.999.txt",h=T)
d<-fread("WS1/outSNPs_Tatry2x0.999.txt",h=T)
e<-fread("WS1/outSNPs_Tatry4x0.999.txt",h=T)
kk<-rbind(a[,c(20,3)],b[,c(20,3)],c[,c(20,3)],d[,c(20,3)],e[,c(20,3)]) 
kk1<-kk[ order(kk[,1],kk[,2]), ]
kk2<-unique(kk1[duplicated(kk1[,c(1,2)]),])
kk3<-kk1[!duplicated(kk1[,c(1,2)]),]
#2. find overlaps in halleri
f<-fread("WS1/outSNPs_AlpsH0.999.txt",h=T)
g<-fread("WS1/outSNPs_FagarasH0.999.txt",h=T)
hh<-rbind(f[,c(20,3)],g[,c(20,3)]) 
hh1<-hh[ order(hh[,1],hh[,2]), ]
hh2<-unique(hh1[duplicated(hh1[,c(1,2)]),])
hh3<-hh1[!duplicated(hh1[,c(1,2)]),]
#3. between species
tt<-rbind(kk3,hh3)
tt1<-tt[ order(tt[,1],tt[,2]), ]
tt2<-unique(tt1[duplicated(tt1[,c(1,2)]),])

##Gene level - quantile
#1. arenosa
a<-subset(a,!a$ann %in% "intragenic_variant" & !a$ann %in% "downstream_gene_variant" & !a$ann %in% "upstream_gene_variant")
b<-subset(b,!b$ann %in% "intragenic_variant" & !b$ann %in% "downstream_gene_variant" & !b$ann %in% "upstream_gene_variant")
c<-subset(c,!c$ann %in% "intragenic_variant" & !c$ann %in% "downstream_gene_variant" & !c$ann %in% "upstream_gene_variant")
d<-subset(d,!d$ann %in% "intragenic_variant" & !d$ann %in% "downstream_gene_variant" & !d$ann %in% "upstream_gene_variant")
e<-subset(e,!e$ann %in% "intragenic_variant" & !e$ann %in% "downstream_gene_variant" & !e$ann %in% "upstream_gene_variant")
a<-subset(a,!a$ID %like% "-")
b<-subset(b,!b$ID %like% "-")
c<-subset(c,!c$ID %like% "-")
d<-subset(d,!d$ID %like% "-")
e<-subset(e,!e$ID %like% "-")
a$one<-1
b$one<-1
c$one<-1
d$one<-1
e$one<-1
a1<-setDT(a)[,list(perGene=sum(one)),by=list(Category=ID)]
b1<-setDT(b)[,list(perGene=sum(one)),by=list(Category=ID)]
c1<-setDT(c)[,list(perGene=sum(one)),by=list(Category=ID)]
d1<-setDT(d)[,list(perGene=sum(one)),by=list(Category=ID)]
e1<-setDT(e)[,list(perGene=sum(one)),by=list(Category=ID)]

a2<-subset(a1,a1$perGene>=quantile(a1$perGene,0.90))
b2<-subset(b1,b1$perGene>=quantile(b1$perGene,0.90))
c2<-subset(c1,c1$perGene>=quantile(c1$perGene,0.90))
d2<-subset(d1,d1$perGene>=quantile(d1$perGene,0.90))
e2<-subset(e1,e1$perGene>=quantile(e1$perGene,0.90))

a2<-subset(a1,a1$perGene>mean(a1$perGene)) # ALTERNATIVE 1
b2<-subset(b1,b1$perGene>mean(b1$perGene))
c2<-subset(c1,c1$perGene>mean(c1$perGene))
d2<-subset(d1,d1$perGene>mean(d1$perGene))
e2<-subset(e1,e1$perGene>mean(e1$perGene))

a2<-subset(a1,a1$perGene>1) # ALTERNATIVE 2
b2<-subset(b1,b1$perGene>1)
c2<-subset(c1,c1$perGene>1)
d2<-subset(d1,d1$perGene>1)
e2<-subset(e1,e1$perGene>1)
kk<-rbind(a2,b2,c2,d2,e2)
kk$perGene<-1
kk1<-kk[ order(kk[,1]), ]
kk2<-unique(kk1[duplicated(kk1$Category),])
kk3<-kk1[!duplicated(kk1[,c(1,2)]),]
i2<-kk2$Category
write.table(a2$Category,'outlierLists/SNPFst_NarrowGenesMin2SNPs.Alps0.txt',col.names = F,row.names = F,quote = F)
write.table(b2$Category,'outlierLists/SNPFst_NarrowGenesMin2SNPs.Fagaras.txt',col.names = F,row.names = F,quote = F)
write.table(c2$Category,'outlierLists/SNPFst_NarrowGenesMin2SNPs.Rodna.txt',col.names = F,row.names = F,quote = F)
write.table(d2$Category,'outlierLists/SNPFst_NarrowGenesMin2SNPs.Tatry2x.txt',col.names = F,row.names = F,quote = F)
write.table(e2$Category,'outlierLists/SNPFst_NarrowGenesMin2SNPs.Tatry4x.txt',col.names = F,row.names = F,quote = F)

#here is the possibility to plor VennDiagr with a2,b2,...
#2. halleri
f<-subset(f,!f$ann %in% "intragenic_variant" & !f$ann %in% "downstream_gene_variant" & !f$ann %in% "upstream_gene_variant")
g<-subset(g,!g$ann %in% "intragenic_variant" & !g$ann %in% "downstream_gene_variant" & !g$ann %in% "upstream_gene_variant")
f$one<-1
g$one<-1
f<-subset(f,!f$ID %like% "-")
g<-subset(g,!g$ID %like% "-")
f1<-setDT(f)[,list(perGene=sum(one)),by=list(Category=ID)]
g1<-setDT(g)[,list(perGene=sum(one)),by=list(Category=ID)]
f2<-subset(f1,f1$perGene>=quantile(f1$perGene,0.90))
g2<-subset(g1,g1$perGene>=quantile(g1$perGene,0.90))

f2<-subset(f1,f1$perGene>mean(f1$perGene)) #ALTERNATIVE 1
g2<-subset(g1,g1$perGene>mean(g1$perGene))

f2<-subset(f1,f1$perGene>1) #ALTERNATIVE 2
g2<-subset(g1,g1$perGene>1)
hh<-rbind(f2,g2)
hh$perGene<-1
hh1<-hh[ order(hh[,1]), ]
hh2<-unique(hh1[duplicated(hh1$Category),])
hh3<-hh1[!duplicated(hh1[,c(1,2)]),]
#3. between species
tt<-rbind(kk3,hh3)
tt1<-tt[ order(tt[,1],tt[,2]), ]
tt2<-unique(tt1[duplicated(tt1[,c(1,2)]),])
tt3<-subset(tt2,!tt2$Category %like% "-")
i2<-tt3$Category
write.table(f2$Category,'outlierLists/SNPFst_NarrowGenesMin2SNPs.AlpsH.txt',col.names = F,row.names = F,quote = F)
write.table(g2$Category,'outlierLists/SNPFst_NarrowGenesMin2SNPs.FagarasH.txt',col.names = F,row.names = F,quote = F)

write.table(tt3$Category,'outlierLists/SNPFst_arenosa_halleri_NarrowGenesMinAvgSNPs.txt',col.names = F,row.names = F,quote = F)
##find overlaps and annotate
ann<-fread("../../../Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20171231.txt",h=T,quote="")
mm3<-subset(ann,ann$`Version-2` %in% i2)
write.table(mm3,'outlierLists/SNPFst_arenosa_halleri_NarrowGenesMinAvgSNPs_ann.txt',col.names = F,row.names = F,quote = F,sep="\t")

################# 3. SELSCAN ######################
### Kubas results - overlap over SNPs ###
library(data.table)
setwd("/home/aa/alpine/arenosaGenome/selScans/")
#1.find any candidate SNP - my data
a<-fread("WS1/outSNPs_Alps0.999.txt",h=T)
b<-fread("WS1/outSNPs_AlpsH0.999.txt",h=T)
c<-fread("WS1/outSNPs_Fagaras0.999.txt",h=T)
d<-fread("WS1/outSNPs_FagarasH0.999.txt",h=T)
e<-fread("WS1/outSNPs_Rodna0.999.txt",h=T)
f<-fread("WS1/outSNPs_Tatry2x0.999.txt",h=T)
g<-fread("WS1/outSNPs_Tatry4x0.999.txt",h=T)
kk<-rbind(a[,c(1,2,3,4)],b[,c(1,2,3,4)],c[,c(1,2,3,4)],d[,c(1,2,3,4)],e[,c(1,2,3,4)],f[,c(1,2,3,4)],g[,c(1,2,3,4)]) 
kk1<-kk[ order(kk[,1],kk[,2]), ]
kk<-kk1[!duplicated(kk1[,c(1,2)]),] 
#   kk<-fread("outSNPsAll_Tatry2x.txt",h=T)
kk$scaff<-substr(kk$ID,3,3)
kk$end<-kk$i.start
kk$POS<-kk$i.start
#2.extract Kubas outliers
h1<-fread("../../selscan/kubaResults/halleri/scaffold_1.phased.HCA.selscan.window900.nsl.100bins_candindate_list_threshold_based.tsv")
h1$scaff<-1
for (i in 2:8){
  h<-fread(paste("../../selscan/kubaResults/halleri/scaffold_",i,".phased.HCA.selscan.window900.nsl.100bins_candindate_list_threshold_based.tsv",sep=""))
  h$scaff<-i
  h1<-rbind(h1,h)}
l1<-fread("../../selscan/kubaResults/halleri/scaffold_1.phased.DRG.selscan.window900.nsl.100bins_candindate_list_threshold_based.tsv")
l1$scaff<-1
for (i in 2:8){
  l<-fread(paste("../../selscan/kubaResults/halleri/scaffold_",i,".phased.DRG.selscan.window900.nsl.100bins_candindate_list_threshold_based.tsv",sep=""))
  l$scaff<-i
  l1<-rbind(l1,l)}
# 3.Delete those selected in both
al<-rbind(h1,l1)
al1<-al[ order(al[,9],al[,2]), ]
al2<-al1[duplicated(al1[,c(9,2)]),] 
hh<-subset(h1, h1$scaff %in% al2$scaff & !h1$POS %in% al2$POS)
hh$end<-hh$POS
hh$scaff<-as.character(hh$scaff)
# 4. find overlap in mine and Kubas
setkey(kk, scaff, POS,end)
a<-foverlaps(x = hh, y = kk, type="any")
aa<-a[complete.cases(a[ , 2]),]
### Kubas results - overlap over genes ###
#   ar<-fread(paste('ann/ALLarenosa.table.recode.txt',sep=''),h=F,na.strings = "-9")
#   ar<-fread(paste('ann/ALLhalleri.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
#   colnames(ar) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",seq(10,ncol(ar),1))
#   hh$scaff<-paste("scaffold_",hh$scaff,sep='')
#   hh1<-hh[ order(hh[,9],hh[,10]), ]
#   hh1$id<-paste(hh1$scaff,hh1$POS,sep="")
#   ar$id<-paste(ar$scaff,ar$start,sep="")
#   ar1<-subset(ar,paste(ar$scaff,ar$start,sep="") %in% paste(hh1$scaff,hh1$POS,sep=""))
#   write.table(x=ar1,file = paste('outSNPsHCA_Kuba.txt',sep=''),append = F,quote = F,col.names = T,row.names = F)
library(data.table)
setwd("/home/aa/alpine/arenosaGenome/selScans/")
ar1<-read.table('outSNPsHCA_Kuba.txt',h=T)
ar1$one<-1
aa<-setDT(ar1)[,list(perGene=sum(one)),by=list(Category=ID)]
#   hist(aa$perGene)
q<-quantile(aa$perGene,probs = 0.99)
aa1<-subset(aa,aa$perGene>q)
#aa1<-subset(aa,aa$perGene>mean(aa$perGene))
i2<-droplevels(aa1$Category)
i3<-i2
#My genes
kk$one<-1
k<-setDT(kk)[,list(perGene=sum(one)),by=list(Category=ID)]
hist(k$perGene)
q<-quantile(k$perGene,probs = 0.50)
k1<-subset(k,k$perGene>q)
i1<-Reduce(intersect, list(k$Category,aa$Category))
i2<-Reduce(intersect, list(k1$Category,aa1$Category))
i2<-droplevels(aa1$Category)

################# 4. BAYPASS OVERLAP ######################
### OVERLAP BAYPASS WITH FST QUARTETS ###
library(data.table)
setwd("/home/aa/alpine/arenosaGenome/selScans/")
lin<-"Fagaras"
for (lin in c("Fagaras","Tatry4x","Tatry2x",'Alps')) {
b<-fread(paste("../../bayPass/",lin,"/outSNPsAll_",lin,".BPann.txt",sep=''),h=T)
f<-fread(paste("quartetOutliers/outSNPsAll_",lin,"0.99.txt",sep=''),h=T)

f1<-subset(f,!f$ann %in% "intragenic_variant" & !f$ann %in% "downstream_gene_variant" & !f$ann %in% "upstream_gene_variant") 
# FULL
colnames(f)[2]<-"start"
kk<-rbind(b[,c(7,4)],f[,c(1,2)]) 
kk1<-kk[ order(kk[,1],kk[,2]), ]
kk<-kk1[!duplicated(kk1[,c(1,2)]),] #MERGED
kk1<-unique(kk1[duplicated(kk1[,c(1,2)]),]) #OVERLAPPED
# NARROW
colnames(f1)[2]<-"start"
mm<-rbind(b[,c(7,4)],f1[,c(1,2)]) 
mm1<-mm[ order(mm[,1],mm[,2]), ]
mm<-mm1[!duplicated(mm1[,c(1,2)]),] #MERGED
mm1<-unique(mm1[duplicated(mm1[,c(1,2)]),]) #OVERLAPPED
mm1<-subset(mm1,!mm1$ID %like% "-")
#   #1. find outlier genes - outlier n. SNPs, union 
#   kk$one<-1
#   k<-setDT(kk)[,list(perGene=sum(one)),by=list(Category=ID)]
#   k1<-subset(k,k$perGene>1)
#   k2<-subset(k,k$perGene>=quantile(k$perGene,0.99))
#  #2. find outlier genes - overlap 
#  kk1$one<-1 #FULL
#  fu<-setDT(kk1)[,list(perGene=sum(one)),by=list(Category=ID)]
#  fu1<-subset(fu,fu$perGene>1) ############## THREE OPTIONS
#  fu2<-subset(fu,fu$perGene>=quantile(fu$perGene,0.90)) ############## THREE OPTIONS
#  fu3<-subset(fu,fu$perGene>mean(fu$perGene)) ############## THREE OPTIONS
mm1$one<-1 #NARROW
o<-setDT(mm1)[,list(perGene=sum(one)),by=list(Category=ID)]
  o1<-subset(o,o$perGene>1) ############## THREE OPTIONS
#  o2<-subset(o,o$perGene>=quantile(o$perGene,0.90)) ############## THREE OPTIONS
#o3<-subset(o,o$perGene>mean(o$perGene)) ############## THREE OPTIONS
i2<-o1$Category
#write.table(i2,paste("bayPass_quartetFst/bfGenes",lin,sep=""),quote = F,col.names = F,row.names = F)
write.table(kk1,paste("bayPass_quartetFst/bfSNPsOverlap_",lin,sep=""),quote = F,col.names = F,row.names = F)
#write.table(k2$Category,paste("bayPass_quartetFst/bfGenes1percUnion",lin,sep=""),quote = F,col.names = F,row.names = F)
write.table(o1$Category,paste("bayPass_quartetFst/bfGenesOverlap_Min2SNPsNarrowGene_",lin,sep=""),quote = F,col.names = F,row.names = F)
#write.table(fu1$Category,paste("bayPass_quartetFst/bfGenesOverlap_Min2SNPsFullGene_",lin,sep=""),quote = F,col.names = F,row.names = F)
}
### OVERLAP BAYPASS WITH SELSCAN ###
setwd("/home/aa/alpine/arenosaGenome/selScans/")
#1.find any candidate SNP - BP
a<-fread("../../bayPass/Alps/outSNPsAll_Alps.BPann.txt",h=T)
b<-fread("../../bayPass/Fagaras/outSNPsAll_Fagaras.BPann.txt",h=T)
c<-fread("../../bayPass/Tatry2x/outSNPsAll_Tatry2x.BPann.txt",h=T)
d<-fread("../../bayPass/Tatry4x/outSNPsAll_Tatry4x.BPann.txt",h=T)
kk<-rbind(a[,c(1,2,3,4)],b[,c(1,2,3,4)],c[,c(1,2,3,4)],d[,c(1,2,3,4)]) 
kk1<-kk[ order(kk[,1],kk[,4]), ]
kk<-kk1[!duplicated(kk1[,c(1,4)]),] 
#   kk<-fread("../../bayPass/Tatry2x/outSNPsAll_Tatry2x.BPann.txt",h=T)[,1:4]
kk$end<-kk$start
ar1<-fread('outSNPsHIGH_Kuba.txt',h=T)
ar1$end<-ar1$start
# 4. find overlap in mine and Kubas
setkey(kk, scaff, start,end)
a<-foverlaps(x = ar1, y = kk, type="any")
aa<-a[complete.cases(a[ , 2]),]
aa$one<-1
k<-setDT(aa)[,list(perGene=sum(one)),by=list(Category=ID)]
k1<-subset(k,k$perGene>1)
i2<-k1$Category

### OVERLAP BAYPASS STD WITH SNP FST ###
library(data.table)
setwd("/home/aa/alpine/arenosaGenome/selScans/")
bb<-fread("../../bayPass/All/cand.txt",h=F)
colnames(bb)[1:2] <- c("scaff","start")
bb$end<-bb$start
a<-fread("WS1/outSNPs_Alps0.999.txt",h=T)
b<-fread("WS1/outSNPs_Fagaras0.999.txt",h=T)
c<-fread("WS1/outSNPs_Rodna0.999.txt",h=T)
d<-fread("WS1/outSNPs_Tatry2x0.999.txt",h=T)
e<-fread("WS1/outSNPs_Tatry4x0.999.txt",h=T)
kk<-rbind(a[,c(20,3)],b[,c(20,3)],c[,c(20,3)],d[,c(20,3)],e[,c(20,3)]) 
kk1<-kk[ order(kk[,1],kk[,2]), ]
kk2<-unique(kk1[duplicated(kk1[,c(1,2)]),])
kk3<-kk1[!duplicated(kk1[,c(1,2)]),]
kk3$scaff<-paste("scaffold_",substr(kk3$ID,3,3),sep="")
kk3$end<-kk3$start

setkey(kk3, scaff, start,end)
a<-foverlaps(x = bb, y = kk3, type="any")
aa<-a[complete.cases(a[ , 2]),]

### SNP FST OVERLAP ###
"lightgreen" - FagarasH
"yellowgreen" - AlpsH
Alps: "#E0912F"
Fagaras: "#065570"
Rodna: "#068FBD"
Tatry4x: "#960C13"
Tatry2x: "#ED7C7F"

setwd("/home/aa/alpine/arenosaGenome/selScans/WS1")
library(VennDiagram)
library(data.table)
a<-fread("outSNPs_Alps0.999.txt",h=T)
b<-fread("outSNPs_AlpsH0.999.txt",h=T)
c<-fread("outSNPs_Fagaras0.999.txt",h=T)
d<-fread("outSNPs_FagarasH0.999.txt",h=T)
e<-fread("outSNPs_Rodna0.999.txt",h=T)
f<-fread("outSNPs_Tatry2x0.999.txt",h=T)
g<-fread("outSNPs_Tatry4x0.999.txt",h=T)
a$x<-paste(a$ID, a$i.start,sep='.')
b$x<-paste(b$ID, b$i.start,sep='.')
c$x<-paste(c$ID, c$i.start,sep='.')
d$x<-paste(d$ID, d$i.start,sep='.')
e$x<-paste(e$ID, e$i.start,sep='.')
f$x<-paste(f$ID, f$i.start,sep='.')
g$x<-paste(g$ID, g$i.start,sep='.')
length(Reduce(intersect, list(a$x,b$x)))
#SNPs
v<-venn.diagram(x=list("Alps"=a$x,"Fagaras"=c$x,"Rodna"=e$x,"Tatry4x"=g$x,"Tatry2x"=f$x),paste("figures/SNPs_FST_Arenosa_venn0.1percent.tiff",sep=""),fill = c("#E0912F","#065570",'#068FBD','#960C13','#ED7C7F'), lty = "blank",alpha=0.3 , cex=0.7, cat.cex=0.6,main = "Candidate SNPs",height = 2000,width = 2000)
# v<-venn.diagram(x=list("AlpsH"=b$x,"Alps"=a$x,"FagarasH"=d$x,"Fagaras"=c$x,"Tatry2x"=f$x),paste("figures/SNPs_Arenosa_Halleri_venn0.1percent.tiff",sep=""),fill = c("yellowgreen","#E0912F","mediumspringgreen","#065570",'#ED7C7F'), lty = "blank",alpha=0.3 , cex=0.6, cat.cex=0.6,main = "SNPs - 0.1% outliers",height = 2000,width = 2000)
v<-venn.diagram(x=list("AlpsH"=i3,"FagarasH"=i2),paste("figures/selScan_Arenosa_Halleri_1percentGenes.tiff",sep=""),fill = c("yellowgreen","mediumspringgreen"), lty = "blank",alpha=0.3 , cex=0.6, cat.cex=0.6,main = "GENES 1% outliers",height = 2000,width = 2000)
#Genes (run upper part first)
v<-venn.diagram(x=list("Alps"=a2$Category,"Fagaras"=b2$Category,"Rodna"=c2$Category,"Tatry4x"=e2$Category,"Tatry2x"=d2$Category),paste("figures/Genes_FST_Arenosa_venn0.1percent.tiff",sep=""),fill = c("#E0912F","#065570",'#068FBD','#960C13','#ED7C7F'), lty = "blank",alpha=0.3 , cex=0.7, cat.cex=0.6,main = "Candidate genes",height = 2000,width = 2000)
Reduce(intersect, list(a2$Category,b2$Category,c2$Category,e2$Category,d2$Category))
v<-venn.diagram(x=list("AlpsH"=f2$Category,"Alps"=a2$Category,"FagarasH"=g2$Category,"Fagaras"=b2$Category,"Tatry2x"=d2$Category),paste("figures/Genes_Arenosa_Halleri_venn0.1percent.tiff",sep=""),fill = c("yellowgreen","#E0912F","mediumspringgreen","#065570",'#ED7C7F'), lty = "blank",alpha=0.3 , cex=0.6, cat.cex=0.6,main = "Genes - 0.1% outliers",height = 2000,width = 2000)
i2<-Reduce(intersect, list(f2$Category,a2$Category,g2$Category,b2$Category,d2$Category))
### OVERLAP GENE SETS - QF_BP###
library(VennDiagram)
setwd("/home/aa/alpine/arenosaGenome/selScans/")
a<-na.omit(read.table("bayPass_quartetFst/bfGenesOverlap_Min2SNPsNarrowGene_Alps",h=F))
b<-na.omit(read.table("bayPass_quartetFst/bfGenesOverlap_Min2SNPsNarrowGene_Fagaras",h=F))
c<-na.omit(read.table("bayPass_quartetFst/bfGenesOverlap_Min2SNPsNarrowGene_Tatry2x",h=F))
d<-na.omit(read.table("bayPass_quartetFst/bfGenesOverlap_Min2SNPsNarrowGene_Tatry4x",h=F))
a1<-a[ order(a[,1]), ]
b1<-b[ order(b[,1]), ]
c1<-c[ order(c[,1]), ]
d1<-d[ order(d[,1]), ]
length(Reduce(intersect, list(b1,c1,d1)))
v<-venn.diagram(x=list("Alps"=a1,"Fagaras"=b1,"Tatry4x"=d1,"Tatry2x"=c1),paste("bayPass_quartetFst/qf_bf_outlierGenes.tiff",sep=""),fill = c("#E0912F","#065570",'#960C13','#ED7C7F'), lty = "blank",alpha=0.3 , cex=1, cat.cex=0.6,main = "Candidate genes",height = 2000,width = 2000)
##find overlaps and annotate
mm<-rbind(a,b,c,d) 
mm1<-as.data.frame(mm[ order(mm[,1]), ])
mm2<-as.data.frame(mm1[duplicated(mm1[,c(1)]),]) 
mm2$one<-1
mm3<-setDT(mm2)[,list(parallels=sum(one)),by=list(Category=`mm1[duplicated(mm1[, c(1)]), ]`)]
##annotate
ann<-fread("../../../Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20171231.txt",h=T,quote="")
mm3<-subset(ann,ann$AL %in% mm3$Category)
write.table(mm3,paste("bayPass_quartetFst/bfGenesOverlap_ann_Narrow_min2SNPs.txt",sep=""),quote = F,col.names = F,row.names = F,sep="\t")
write.table(mm3,paste("bayPass_quartetFst/bfGenesOverlap_ann_Narrow_minAvgSNPs_Alps.txt",sep=""),quote = F,col.names = F,row.names = F,sep = "\t")
i2<-mm3$`Version-2`
### OVERLAP SNPS SETS - QF_BP###
library(VennDiagram)
a<-na.omit(read.table("bayPass_quartetFst/bfSNPsOverlap_Alps",h=F))
b<-na.omit(read.table("bayPass_quartetFst/bfSNPsOverlap_Fagaras",h=F))
c<-na.omit(read.table("bayPass_quartetFst/bfSNPsOverlap_Tatry2x",h=F))
d<-na.omit(read.table("bayPass_quartetFst/bfSNPsOverlap_Tatry4x",h=F))
a1<-paste(a$V1,a$V2,sep="")
b1<-paste(b$V1,b$V2,sep="")
c1<-paste(c$V1,c$V2,sep="")
d1<-paste(d$V1,d$V2,sep="")
v<-venn.diagram(x=list("Alps"=a1,"Fagaras"=b1,"Tatry4x"=d1,"Tatry2x"=c1),paste("bayPass_quartetFst/qf_bf_outlierSNPs.tiff",sep=""),fill = c("#E0912F","#065570",'#960C13','#ED7C7F'), lty = "blank",alpha=0.3 , cex=0.7, cat.cex=0.6,main = "Candidate SNPs",height = 2000,width = 2000)
##find overlaps
mm<-rbind(a,b,c,d) 
mm1<-as.data.frame(mm[ order(mm[,1],mm[,2]), ])
mm2<-as.data.frame(mm1[duplicated(mm1[,c(1,2)]),]) 
mm2$one<-1
mm3<-setDT(mm2)[,list(parallels=sum(one)),by=list(Category=paste(V1,V2,sep=""))]


### OVERLAP OUTLIER WINDOWS ###
setwd("/home/aa/alpine/arenosaGenome/selScans/")
library(VennDiagram)
library(dplyr)
library(readr)
options(scipen=999)
a<-read.table("out_1percent_negOut_moreThan200Sites_Tatry2x.txt",h=T)
b<-read.table("out_1percent_negOut_moreThan200Sites_Fagaras.txt",h=T)
a1<-as.character(interaction(a$scaff, a$start))
b1<-as.character(interaction(b$scaff, b$start))
length(Reduce(intersect, list(a1,b1)))
#plot
v<-venn.diagram(x=list("Tatry 2x"=a1,"Fagaras"=b1),paste("Tatry2x_Fagaras_venn1percent.tiff",sep=""), lty = "blank", fill = c("red2","blue1"),alpha=0.3 , cex=1, cat.cex=0.6,main = "1% outliers",height = 2000,width = 2000)
#get interval list
aa<-Reduce(intersect, list(a1,b1)) 
aa1<-gsub('.', ':', aa,fixed = T)
aa<-as.data.frame(aa1)
aa$pos<-as.character(as.numeric(substr(aa$aa1,12,30))+15000)
write.table(aa,paste("Tatry2x_Fagaras_1percentOverlaps.intervals",sep=""),quote = F,col.names = F,row.names = F,append = F,sep = "-")


### PROCESS PARALLEL CANDIDATES###
library(data.table)
setwd("/home/aa/alpine/arenosaGenome/selScans/")
#1.find parallel SNPs
a<-fread("outSNPsAll_Alps.txt",h=T)
b<-fread("outSNPsAll_AlpsH.txt",h=T)
c<-fread("outSNPsAll_Fagaras.txt",h=T)
d<-fread("outSNPsAll_FagarasH.txt",h=T)
e<-fread("outSNPsAll_Rodna.txt",h=T)
f<-fread("outSNPsAll_Tatry2x.txt",h=T)
g<-fread("outSNPsAll_Tatry4x.txt",h=T)
kk<-rbind(a[,c(1,2)],b[,c(1,2)],c[,c(1,2)],d[,c(1,2)],e[,c(1,2)],f[,c(1,2)],g[,c(1,2)]) 
kk1<-kk[ order(kk[,1],kk[,2]), ]
kk2<-kk1[duplicated(kk1[,c(1,2)]),] 
kk3<-kk2[!duplicated(kk2[,c(1,2)]),] 
al<-unique(kk3$ID)
#2.extract them from file
ar<-fread(paste('ann/ALLarenosa.table.recode.txt',sep=''),h=F,na.strings = "-9")
ha<-fread(paste('ann/ALLhalleri.table.recode.txt',sep=''),h=F,na.strings = "-9")
colnames(ar) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",seq(10,ncol(ar),1))
colnames(ha) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",seq(10,ncol(ha),1))
ar1<-subset(ar,ar$ID %in% kk3$ID & ar$start %in% kk3$i.start)
ha1<-subset(ha,ha$ID %in% kk3$ID & ha$start %in% kk3$i.start)
ar1$end<-ar1$start
ha1$end<-ha1$start
setkey(ha1, scaff, start,end)
a<-foverlaps(x = ar1, y = ha1, type="any")
write.table(x=a,file = paste('outSNPsAll_parallel.txt',sep=''),append = F,quote = F,col.names = T,row.names = F)
#3. calculate AF for each lineage
i<-which(colnames(a) %in% "10")
afh1<-a[,i:as.numeric((i-1)+8)]
a$ACh1<-rowSums(afh1,na.rm = T)
a$NAh1<-apply(is.na(afh1), 1, sum)
a$ANh1<-(8-a$NAh1)*2
a$OBI<-a$ACh1/a$ANh1
afh2<-a[,as.numeric(i+8):as.numeric((i-1)+16)]
a$ACh2<-rowSums(afh2,na.rm = T)
a$NAh2<-apply(is.na(afh2), 1, sum)
a$ANh2<-(8-a$NAh2)*2
a$GUN<-a$ACh2/a$ANh2
aff1<-a[,as.numeric(i+16):as.numeric((i-1)+24)]
a$ACf1<-rowSums(aff1,na.rm = T)
a$NAf1<-apply(is.na(aff1), 1, sum)
a$ANf1<-(8-a$NAf1)*2
a$HCA<-a$ACf1/a$ANf1
aff2<-a[,as.numeric(i+24):as.numeric((i-1)+32)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*2
a$DRG<-a$ACf2/a$ANf2
i<-which(colnames(a) %in% "i.10")
afh1<-a[,i:as.numeric((i-1)+8)]
a$ACh1<-rowSums(afh1,na.rm = T)
a$NAh1<-apply(is.na(afh1), 1, sum)
a$ANh1<-(8-a$NAh1)*4
a$ING<-a$ACh1/a$ANh1
afh2<-a[,as.numeric(i+8):as.numeric((i-1)+16)]
a$ACh2<-rowSums(afh2,na.rm = T)
a$NAh2<-apply(is.na(afh2), 1, sum)
a$ANh2<-(8-a$NAh2)*2
a$VEL<-a$ACh2/a$ANh2
aff1<-a[,as.numeric(i+16):as.numeric((i-1)+24)]
a$ACf1<-rowSums(aff1,na.rm = T)
a$NAf1<-apply(is.na(aff1), 1, sum)
a$ANf1<-(8-a$NAf1)*2
a$ZEP<-a$ACf1/a$ANf1
aff2<-a[,as.numeric(i+24):as.numeric((i-1)+32)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*2
a$SUB<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+32):as.numeric((i-1)+40)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*2
a$BAB<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+40):as.numeric((i-1)+47)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(7-a$NAf2)*4
a$SCH<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+47):as.numeric((i-1)+55)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$WIL<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+55):as.numeric((i-1)+63)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$KAS<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+63):as.numeric((i-1)+71)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$LAC<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+71):as.numeric((i-1)+79)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$BAL<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+79):as.numeric((i-1)+87)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$DRA<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+87):as.numeric((i-1)+95)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$TIS<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+95):as.numeric((i-1)+102)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(7-a$NAf2)*4
a$INE<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+102):as.numeric((i-1)+110)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$CAR<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+110):as.numeric((i-1)+118)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$TKO<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+118):as.numeric((i-1)+126)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$TRT<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+126):as.numeric((i-1)+134)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$HRA<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+134):as.numeric((i-1)+142)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$SPI<-a$ACf2/a$ANf2
s4<-select(a,i.ID,i.ann,VEL, ZEP, TKO, TRT, LAC, BAL, INE,SCH,WIL,SUB, BAB,HRA,SPI,DRA,TIS,CAR,ING,KAS,OBI,HCA,GUN,DRG)
#4. plot
library(gplots)
library(RColorBrewer)
pops<-c("VEL", 'ZEP', 'TKO', 'TRT', 'LAC', 'BAL', 'INE','SCH','WIL','SUB', 'BAB','HRA','SPI','DRA','TIS','CAR','ING','KAS','OBI','HCA','GUN','DRG')
my_palette <- colorRampPalette(c("khaki1", "green2", "blue3"))(n = 20)
df<-as.matrix(s4[,3:24],rownames = paste(s4$i.ID,s4$i.ann,sep=" - "))
pdf(paste("figures/heatmap_parallels.0.99.pdf",sep=""),height = nrow(s4)/10,width = 14)
heatmap.2(x = df,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,2,4,6,7,9,11,13,15,16,18,19,20,21),sepcolor= c("black","orange","orange","orange","orange","black","green","green","green","green","black","orange","black","green"),sepwidth = c(0.05),trace="none",ColSideColors = c("orange","orange","orange","orange","orange","orange","orange","orange","orange","green","green","green","green","green","green","green","green","green","orange","orange","green","green"),labCol=pops,colCol= c("orange","orange","orange","orange","orange","orange","orange","orange","orange","green","green","green","green","green","green","green","green","green","orange","orange","green","green"),offsetRow = c(0.05),offsetCol= c(0.05),cexRow = c(1),cexCol = c(2),margins = c(5,25),lhei = c(0.5,30),lwid = c(0.05,0.8))
dev.off()


#GO enrichment - TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library(topGO)
library(topGO)
library(ALL)
data(ALL)
data(geneList)
db=read.table("/home/aa/Desktop/references/lyrataV2/mart_Alyrata_V2.txt")
geneID2GO <- readMappings(file = db)
aaa1<-na.omit(aaa)
all <- as.factor(aaa1$ID)
names(all) <- aaa1$p1_Fst
sampleGOdata <- new("topGOdata",description = "Rodna", ontology = "BP", allGenes = all, geneSel = o3, nodeSize = 10, annFUN.file(file = "/home/aa/Desktop/references/lyrataV2/mart_Alyrata_V2.txt",whichOnto = "BP"))
annotationFun = annFUN.file("/home/aa/Desktop/references/lyrataV2/mart_Alyrata_V2.txt")


### PERMUTATION TEST MAJDA ###
setwd("/home/aa/alpine/arenosaGenome/selScans/")
library(dplyr)
#To calculate the number of overlapping candidate windows in real data:
a<-read.table("out_1percent_negOut_moreThan200Sites_Tatry2x.txt",h=T) # candidate list Tatry 2x
b<-read.table("out_1percent_negOut_moreThan200Sites_Fagaras.txt",h=T) # candidate list Fagaras
a1<-as.character(interaction(a$scaff, a$start))
b1<-as.character(interaction(b$scaff, b$start))
p<-length(Reduce(intersect, list(a1,b1))) # number of overlapping windows
a<-read.table("rankFstTatry2x.txt",h=T) # all analysed windows Tatry 2x
b<-read.table("rankFstFagaras.txt",h=T) # all analysed windows Fagaras
aa<-nrow(read.table("out_1percent_negOut_moreThan200Sites_Tatry2x.txt",h=T)) # candidate list Tatry 2x
bb<-nrow(read.table("out_1percent_negOut_moreThan200Sites_Fagaras.txt",h=T)) # candidate list Fagaras
it<-100 # number of random draws
o<-character()
for(i in 1:it){
  a1<-sample_n(a, aa)
  b1<-sample_n(b, bb)
  o <- c(o, length(Reduce(intersect, list(as.character(interaction(a1$scaff, a1$start)),as.character(interaction(b1$scaff, b1$start)))))) # number of overlapping windows between two random draws
}
results<-prop.test(x=as.numeric(o), n = rep (aa+bb,it),p = rep(p/(aa+bb),it)) # two sided test of given proportions

### PERMUTATION TEST BEN LEANEN ###
#there are two functions. One to create the distribution, on another to calculate the pvalue (two sided test)
sampling_prob <- function(n1, n2, N){
  x <- 1:N
  a <- sample(x, n1, replace = F)
  b <- sample(x, n2, replace = F)
  return(sum(a %in% b))
}
pvalue_observed_sampling <- function(distribution, observed){
  pvalue <- 2 * min(sum(observed > distribution) , sum(observed <= distribution)) / length(distribution)
  if(pvalue == 0 ) pvalue <- 1 / length(distribution)
  return(pvalue)
}
distribution <- sapply(1:10000, function(i) sampling_prob(n1=1099,n2=449, N=30000))
hist(distribution)
observed = 32
pvalue_observed_sampling(distribution, observed)

### PERMUTATION TEST ON 4 POPS ###
#there are two functions. One to create the distribution, on another to calculate the pvalue (two sided test)
sampling_prob <- function(n1, n2, n3, n4, N){
  x <- 1:N
  a <- sample(x, n1, replace = F)
  b <- sample(x, n2, replace = F)
  c <- sample(x, n3, replace = F)
  d <- sample(x, n4, replace = F)
  return(sum(a %in% b) + sum(a %in% c) + sum(a %in% d) + sum(b %in% c) + sum(b %in% d) + sum(c %in% d))
}
pvalue_observed_sampling <- function(distribution, observed){
  pvalue <- 2 * min(sum(observed > distribution) , sum(observed <= distribution)) / length(distribution)
  if(pvalue == 0 ) pvalue <- 1 / length(distribution)
  return(pvalue)
}
distribution <- sapply(1:10000, function(i) sampling_prob(n1=277,n2=1109,n3=729,n4=1174, N=11000000))
hist(distribution,nclass = 100)
observed = 116
pvalue_observed_sampling(distribution, observed)

#Proportion:
21/(32+69+55+67)
116/(277+1109+729+1174)



### CHECK FST DATA ###
setwd("/home/aa/alpine/arenosaGenome/selScans/")
region<-"Rodna"
a<-read.table(paste("QuartetSummary",region,".txt",sep=""),h=T)
a$sum1<-(as.numeric(a$wh1)+as.numeric(a$wf1))
a$sum2<-(as.numeric(a$wh1)+as.numeric(a$wf2))
a$sum3<-(as.numeric(a$wh2)+as.numeric(a$wf1))
a$sum4<-(as.numeric(a$wh2)+as.numeric(a$wf2))
pdf("checkFstDataFagaras.pdf",width = 12,height = 8)
plot(a$p1_Fst~a$p1_mis,ylab = "Fst - p1",xlab="number of sites/window",xlim=c(0,14000),ylim=c(-0.08,1))
plot(a$p2_Fst~a$p2_mis,ylab = "Fst - p2",xlab="number of sites/window",xlim=c(0,14000),ylim=c(-0.08,1))
plot(a$p3_Fst~a$p3_mis,ylab = "Fst - p3",xlab="number of sites/window",xlim=c(0,14000),ylim=c(-0.08,1))
plot(a$p4_Fst~a$p4_mis,ylab = "Fst - p4",xlab="number of sites/window",xlim=c(0,14000),ylim=c(-0.08,1))
plot(a$p1_Fst~a$p1_mis,ylab = "Fst - p1",xlab="number of sites/window - xlim=2000",xlim=c(0,2000),ylim=c(-0.08,1))
plot(a$p2_Fst~a$p2_mis,ylab = "Fst - p2",xlab="number of sites/window - xlim=2000",xlim=c(0,2000),ylim=c(-0.08,1))
plot(a$p3_Fst~a$p3_mis,ylab = "Fst - p3",xlab="number of sites/window - xlim=2000",xlim=c(0,2000),ylim=c(-0.08,1))
plot(a$p4_Fst~a$p4_mis,ylab = "Fst - p4",xlab="number of sites/window - xlim=2000",xlim=c(0,2000),ylim=c(-0.08,1))
bpar<-par()
par(mfrow=c(2,1))
par(mar = c(4.1, 4.1, 0.1, 0.1))
plot(a$p1_Fst~a$wh1,ylab = "Fst - p1",xlab="nucleotide diversity in h1",xlim=c(0,0.5),ylim=c(-0.08,1))
plot(a$p1_Fst~a$wf1,ylab = "Fst - p1",xlab="nucleotide diversity in f1",xlim=c(0,0.5),ylim=c(-0.08,1))
plot(a$p2_Fst~a$wh1,ylab = "Fst - p2",xlab="nucleotide diversity in h1",xlim=c(0,0.5),ylim=c(-0.08,1))
plot(a$p2_Fst~a$wf2,ylab = "Fst - p2",xlab="nucleotide diversity in f2",xlim=c(0,0.5),ylim=c(-0.08,1))
plot(a$p3_Fst~a$wh2,ylab = "Fst - p3",xlab="nucleotide diversity in h2",xlim=c(0,0.5),ylim=c(-0.08,1))
plot(a$p3_Fst~a$wf1,ylab = "Fst - p3",xlab="nucleotide diversity in f1",xlim=c(0,0.5),ylim=c(-0.08,1))
plot(a$p4_Fst~a$wh2,ylab = "Fst - p4",xlab="nucleotide diversity in h2",xlim=c(0,0.5),ylim=c(-0.08,1))
plot(a$p4_Fst~a$wf2,ylab = "Fst - p4",xlab="nucleotide diversity in f2",xlim=c(0,0.5),ylim=c(-0.08,1))
par(bpar)
m<- lm(a$p1_Fst~a$sum1)
plot(a$p1_Fst~a$sum1,ylab = "Fst - p1",xlab="sum od nucl. diversity h1+f1")
abline(m, col = "red")
m<- lm(a$p2_Fst~a$sum2)
plot(a$p2_Fst~a$sum2,ylab = "Fst - p2",xlab="sum od nucl. diversity h1+f2")
abline(m, col = "red")
m<- lm(a$p3_Fst~a$sum3)
plot(a$p3_Fst~a$sum3,ylab = "Fst - p3",xlab="sum od nucl. diversity h2+f1")
abline(m, col = "red")
m<- lm(a$p4_Fst~a$sum4)
plot(a$p4_Fst~a$sum4,ylab = "Fst - p4",xlab="sum od nucl. diversity h2+f2")
abline(m, col = "red")
par(mfrow=c(2,1))
par(mar = c(4.1, 4.1, 0.1, 0.1))
plot(a$p1_Fst~a$rh1,ylab = "Fst - p1",xlab="genotypic correlation in h1",xlim=c(0,0.6),ylim=c(-0.08,1))
plot(a$p1_Fst~a$rf1,ylab = "Fst - p1",xlab="genotypic correlation in f1",xlim=c(0,0.6),ylim=c(-0.08,1))
plot(a$p2_Fst~a$rh1,ylab = "Fst - p2",xlab="genotypic correlation in h1",xlim=c(0,0.6),ylim=c(-0.08,1))
plot(a$p2_Fst~a$rf2,ylab = "Fst - p2",xlab="genotypic correlation in f2",xlim=c(0,0.6),ylim=c(-0.08,1))
plot(a$p3_Fst~a$rh2,ylab = "Fst - p3",xlab="genotypic correlation in h2",xlim=c(0,0.6),ylim=c(-0.08,1))
plot(a$p3_Fst~a$rf1,ylab = "Fst - p3",xlab="genotypic correlation in f1",xlim=c(0,0.6),ylim=c(-0.08,1))
plot(a$p4_Fst~a$rh2,ylab = "Fst - p4",xlab="genotypic correlation in h2",xlim=c(0,0.6),ylim=c(-0.08,1))
plot(a$p4_Fst~a$rf2,ylab = "Fst - p4",xlab="genotypic correlation in f2",xlim=c(0,0.6),ylim=c(-0.08,1))
dev.off()
pdf("histFstFagaras.pdf")
hist(a$p1_Fst,breaks = 100,xlim=c(-0.2,1))
hist(a$p2_Fst,breaks = 100,xlim=c(-0.2,1))
hist(a$p3_Fst,breaks = 100,xlim=c(-0.2,1))
hist(a$p4_Fst,breaks = 100,xlim=c(-0.2,1))
hist(a$n1_Fst,breaks = 100,xlim=c(-0.2,1))
hist(a$n2_Fst,breaks = 100,xlim=c(-0.2,1))
dev.off()

### OPTIMAL CUTOFF ###
a<-read.table("rankFstTatry2x.txt",h=T)
b<-read.table("rankFstFagaras.txt",h=T)
a1<-subset(a,a$out_n1 %in% 0 & a$out_n2 %in% 0)
a2<-subset(a1,a1$p1_mis>200 & a1$p2_mis>200 & a1$p3_mis>200 & a1$p4_mis>200)
b1<-subset(b,b$out_n1 %in% 0 & b$out_n2 %in% 0)
b2<-subset(b1,b1$p1_mis>200 & b1$p2_mis>200 & b1$p3_mis>200 & b1$p4_mis>200)
pro<-seq(0.00001,0.05,0.00001)#0.05,0.04,0.03,0.02,
m<-matrix(nrow = length(pro),ncol = 2)
m[,1]<-pro
for(pr in pro){
  ap<-a2[1:(nrow(a2)*pr),]
  bp<-b2[1:(nrow(b2)*pr),]
  a1<-as.character(interaction(ap$scaff, ap$start))
  b1<-as.character(interaction(bp$scaff, bp$start))
  ov<-length(Reduce(intersect, list(a1,b1)))
  tot<-nrow(ap)+nrow(bp)
  iM<-which(pro %in% paste(pr))
  m[iM,2]<- ov/tot}
pdf("Tatry2x_Fagaras_overlapsOverCutoffs.pdf",width = 14,height = 14)
plot(m[,2]~m[,1],xlab = "Cut-off",ylab="Percentage of overlaps")
dev.off()




