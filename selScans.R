
################# selectION SCANS PIPELINE ######################

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
summary_matrix('LAC','BAL','DRA','TIS','FGgaras')
summary_matrix('SCH','WIL','ING','KAS','Alps')

### MAKE A SUMMARY MATRIX - 2 pops ###
setwd("/home/aa/alpine/arenosaGenome/selScans/")
library(data.table)
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
p1<-fread(paste("BPM/",h1,f1,"_WS1000_MS1_BPM.txt",sep=""),h=T)
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

### MAKE A SUMMARY MATRIX - HALLERI/ARENOSA ALL ###
setwd("/home/aa/alpine/arenosaGenome/selScans/")
library(data.table)
h1<-"HIG"
f1<-"LOW"
region<-"Halleri"
region<-"Arenosa"
#read in all windows, make a list 
p1<-fread(paste("BPM/",h1,f1,"_WS1000_MS1_BPM.txt",sep=""),h=T) #h
a<-rbind(p1[,c(2,3,6,9)])
a1<-a[ order(a[,1], a[,2]), ]
a1<-subset(a1,!a1$scaff %in% "Genome")
a<-a1[!duplicated(a1[,c('scaff','start')]),] 
m<-matrix(nrow = nrow(a), ncol = 6,dimnames =list(c(),c("scaff","start","p1_Fst","p1_mis",'wf1','wh1')))
m[,1]<-as.character(a$scaff)
m[,2]<-a$start
m[,3]<-a$FstWC
m[,4]<-a$num_sites
write.table(m,append = F,file = paste("QuartetSummary",region,".txt",sep=""),quote = F, sep = "\t",col.names = T,row.names = F)

### RANK THE WINDOWS - 4 pops ###
setwd("/home/aa/alpine/arenosaGenome/selScans/")
library(data.table)
region<-"Tatry4x"
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
# summary(aa$genic)
# hist(aa$genic)
# table (aa$genic)
# a$genic<-aa$genic
a1<-subset(a,a$out_n1 %in% 0 & a$out_n2 %in% 0)
a2<-subset(a1,a1$p1_mis>200 & a1$p2_mis>200 & a1$p3_mis>200 & a1$p4_mis>200)
#o1<-a[1:(nrow(a)*0.01),]
o1<-a[1:(nrow(a)*0.05),] #
#pdf(paste("genicVsFst_",region,".pdf",sep=""),width = 12,height = 8)
#plot(a2$genic,a2$p1_Fst,ylab = "Fst",xlab="bp within CDSs")
#abline(h = min(o1$p1_Fst))
#plot(a2$genic,a2$p1_Fst,ylab = "Fst",xlab="bp within CDSs",xlim=c(0,1000))
#abline(h = min(o1$p1_Fst))
#dev.off()
write.table(a,paste("rankFst",region,".txt",sep=""),quote=F,row.names = F,sep="\t")
#write.table(o1,paste("out_1percent_negOut_moreThan200Sites_",region,".txt",sep=''),quote=F,row.names = F,sep="\t")
write.table(o1,paste("out_5percent_negOut_moreThan200Sites_",region,".txt",sep=''),quote=F,row.names = F,sep="\t") #
#overlap with genes
o2<-foverlaps(o1, cds, type="any")
o3<-unique(o2$ID)
#write.table(o3,paste("out_1percent_IDs_",region,".txt",sep=''),quote=F,row.names = F,sep="\t")
write.table(o3,paste("out_5percent_IDs_",region,".txt",sep=''),quote=F,row.names = F,sep="\t") #

### RANK THE WINDOWS - 2 pops ###
require(data.table)
setwd("/home/aa/alpine/arenosaGenome/selScans/")
region<-"Arenosa"
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
#summary(aa$genic)
#hist(aa$genic)
#table (aa$genic)
a$genic<-aa$genic
a<-a[order(a$p1_Fst,decreasing = T),]
a2<-subset(a,a$p1_mis>200)
# a2<-subset(a2,a2$genic>0)
#o1<-a2[1:(nrow(a2)*0.01),]
o1<-a2[1:(nrow(a2)*0.05),] #
#pdf(paste("genicVsFst_",region,".pdf",sep=""),width = 12,height = 8)
#plot(a2$genic,a2$p1_Fst,ylab = "Fst",xlab="bp within CDSs")
#abline(h = min(o1$p1_Fst))
#plot(a2$genic,a2$p1_Fst,ylab = "Fst",xlab="bp within CDSs",xlim=c(0,1000))
#abline(h = min(o1$p1_Fst))
#dev.off()
write.table(a,paste("rankFst",region,".txt",sep=""),quote=F,row.names = F,sep="\t")
#write.table(o1,paste("out_1percent_negOut_moreThan200Sites_",region,".txt",sep=''),quote=F,row.names = F,sep="\t")
write.table(o1,paste("out_5percent_negOut_moreThan200Sites_",region,".txt",sep=''),quote=F,row.names = F,sep="\t") #
#overlap with genes
o2<-foverlaps(o1, cds, type="any")
o3<-unique(o2$ID)
#write.table(o3,paste("out_1percent_IDs_",region,".txt",sep=''),quote=F,row.names = F,sep="\t")
write.table(o3,paste("out_5percent_IDs_",region,".txt",sep=''),quote=F,row.names = F,sep="\t") #

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
nf1<-8
nf2<-8
ploidy<-4
pops<-c('SCH','WIL','KAS','ING')
quant<-0.99
setwd("/home/aa/alpine/arenosaGenome/selScans/")
heatmapCand(c('TKO','TRT','HRA','SPI'),"Tatry4x","TTE",8,8,8,8,4,0.95)
heatmapCand(c('VEL','ZEP','SUB','BAB'),"Tatry2x","TDI",8,8,8,8,2,0.95)
heatmapCand(c('SCH','WIL','KAS','ING'),"Alps","ALParenosa",7,8,8,8,4,0.95)
heatmapCand(c('LAC','BAL','DRA','TIS'),"Fagaras","FAGarenosa",8,8,8,8,4,0.95)
heatmapCand2pop(c('INE','CAR'),"Rodna","ROD",7,8,4,0.95)
heatmapCand2pop(c('OBI','GUN'),"AlpsH","ALPhalleri",8,8,2,0.95)
heatmapCand2pop(c('HCA','DRG'),"FagarasH","FAGhalleri",8,8,2,0.95)
heatmapCand2pop(c('HIG','LOW'),"Halleri","ALLhalleri",16,16,2,0.95)


heatmapCand<-function(pops,region,sh,nh1,nh2,nf1,nf2,ploidy,quant){
  library(data.table)
  library(stringr)
  library(dplyr)
  s<-fread(paste('ann/',sh,'.table.recode.txt',sep=''),h=F,na.strings = "-9")
#  o<-fread(paste("quartetOutliers/out_1percent_negOut_moreThan200Sites_",region,".txt",sep = ""),h=T)
  o<-fread(paste("quartetOutliers/out_5percent_negOut_moreThan200Sites_",region,".txt",sep = ""),h=T) #
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
  s4<-dplyr::select(a,ID,i.start,ann,aas,AFh1,AFh2,AFf1,AFf2)
  s4$p1<-abs(s4$AFh1-s4$AFf1) 
  s4$p2<-abs(s4$AFh1-s4$AFf2) 
  s4$p3<-abs(s4$AFh2-s4$AFf1) 
  s4$p4<-abs(s4$AFh2-s4$AFf2) 
  s4$aftot<-abs((s4$AFh1+s4$AFh2)/2-(s4$AFf1+s4$AFf2)/2)
    #   s5<-subset(s4,s4$p1>0.8 | s4$p2>0.8 | s4$p3>0.8 | s4$p4>0.8)
  q<-quantile(s4$aftot,probs = quant,na.rm = T)
  s5<-subset(s4,s4$aftot >= q)
  write.table(x=s5,file = paste('outSNPsAll_',region,quant,'.txt',sep=''),append = F,quote = F,col.names = T,row.names = F)
     } #############alternative end
  #4. plot
  library(gplots)
  library(RColorBrewer)
  my_palette <- colorRampPalette(c("khaki1", "green2", "blue3"))(n = 20)
  df<-as.matrix(s5[,5:8],rownames = paste(s5$ID,s5$ann,sep="         "))
  pdf(paste("figures/heatmap",region,quant,".pdf",sep=""),height = nrow(s5)/10,width = 7)
  heatmap.2(x = df,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,2,4),sepcolor= c("black"),sepwidth = c (0.05),trace="none",ColSideColors = c("orange","orange","green","green"),labCol=pops,colCol= c("orange","orange","green","green"),offsetRow = c(0.05),offsetCol= c(0.05),cexRow = c(1),cexCol = c(2),margins = c(5,25),lhei = c(0.5,30),lwid = c(0.05,0.8))
  dev.off()
#  } #############alternative end
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
region<-"Rodna"
sh<-"ROD"
nh1<-7
nf1<-8
ploidy<-4
quant<-0.99
setwd("/home/aa/alpine/arenosaGenome/selScans/")
heatmapCand2pop(c('INE','CAR'),"Rodna","ROD",7,8,4,0.95)
heatmapCand2pop(c('OBI','GUN'),"AlpsH","ALPhalleri",8,8,2,0.95)
heatmapCand2pop(c('HCA','DRG'),"FagarasH","FAGhalleri",8,8,2,0.95)
heatmapCand2pop(c('HIG','LOW'),"Halleri","ALLhalleri",16,16,2,0.95)


heatmapCand2pop<-function(pops,region,sh,nh1,nf1,ploidy,quant){
  library(data.table)
  library(stringr)
  library(dplyr)
  s<-fread(paste('ann/',sh,'.table.recode.txt',sep=''),h=F,na.strings = "-9")
 # o<-fread(paste("quartetOutliers/out_1percent_negOut_moreThan200Sites_",region,".txt",sep = ""),h=T)
  o<-fread(paste("quartetOutliers/out_5percent_negOut_moreThan200Sites_",region,".txt",sep = ""),h=T) #
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
  s4<-dplyr::select(a,ID,i.start,ann,aas,AFh1,AFf1)
  s4$aftot<-abs(s4$AFh1-s4$AFf1) 
  q<-quantile(s4$aftot,probs = quant,na.rm = T)
  s5<-subset(s4,s4$aftot >= q)
  write.table(x=s5,file = paste('outSNPsAll_',region,quant,'.txt',sep=''),append = F,quote = F,col.names = T,row.names = F)
} #############alternative end
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

### PROCESS CANDIDATES ARENOSA ALL (halleri can be easily done with 2 pop) ### 
library(data.table)
library(stringr)
library(dplyr)
pops<-c('HIG','LOW')
region<-"Arenosa"
sh<-"ALLarenosa"
quant<-0.99
setwd("/home/aa/alpine/arenosaGenome/selScans/")
s<-fread(paste('ann/',sh,'.table.recode.txt2',sep=''),h=F,na.strings = "-9")
o<-fread(paste("quartetOutliers/out_1percent_negOut_moreThan200Sites_",region,".txt",sep = ""),h=T)
o$rank<-seq(1,nrow(o),1)
o1<-o[ order(o[,1], o[,2]), ]
#1.extract only outlier sites
# colnames(s) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",seq(10,ncol(s),1))
s$zero<-rowSums(s[,10:27] == 0,na.rm = T)
colnames(s) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",seq(10,(ncol(s)-2),1),"end","zero")
s<-subset(s,s$zero < 12)
write.table(s,'s.txt')
#  s$end<-s$start
setkey(o1, scaff, start, end)
aaa<-foverlaps(x = s, y = o1, type="within")
a<-aaa[complete.cases(aaa[ , 2]),]
#write.table(a,"overlap.Arenosa.txt")
#2. calculate AF for each lineage
i<-which(colnames(a) %in% "10")
afh1<-a[,c(19,20,23,24,26,27,30,32,33)]
a$ACh1<-rowSums(afh1,na.rm = T)
a$NAh1<-apply(is.na(afh1[,1:2]), 1, sum)
a$NAh2<-apply(is.na(afh1[,3:9]), 1, sum)
a$ANh1<-((2-a$NAh1)*2)+((7-a$NAh2)*4)
a$AFh1<-a$ACh1/a$ANh1
aff1<-a[,c(18,21,22,25,28,29,31,34,35)]
a$ACf1<-rowSums(aff1,na.rm = T)
a$NAf1<-apply(is.na(aff1[,2:3]), 1, sum)
a$NAf2<-apply(is.na(aff1[,c(1,4:9)]), 1, sum)
a$ANf1<-((2-a$NAf1)*2)+((7-a$NAf2)*4)
a$AFf1<-a$ACf1/a$ANf1
s4<-dplyr::select(a,ID,i.start,ann,aas,AFh1,AFf1)
s4$aftot<-abs(s4$AFh1-s4$AFf1) 
s4<-s4[order(s4$aftot,decreasing = T),]
# a2<-subset(a2,a2$genic>0)
s5<-s4[1:824,]
write.table(x=s5,file = paste('quartetOutliers/outSNPsAll_',region,quant,'.txt',sep=''),append = F,quote = F,col.names = T,row.names = F)
#4. plot
library(gplots)
library(RColorBrewer)
my_palette <- colorRampPalette(c("khaki1", "green2", "blue3"))(n = 20)
df<-as.matrix(s5[,5:6],rownames = paste(s5$ID,s5$ann,sep="         "))
pdf(paste("figures/heatmap",region,quant,".pdf",sep=""),height = nrow(s5)/10,width = 7)
heatmap.2(x = df,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,1,2),sepcolor= c("black"),sepwidth = c (0.05),trace="none",ColSideColors = c("orange","green"),labCol=pops,colCol= c("orange","green"),offsetRow = c(0.05),offsetCol= c(0.05),cexRow = c(1),cexCol = c(2),margins = c(5,25),lhei = c(0.5,30),lwid = c(0.05,0.8))
dev.off()

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
write.table(a2$Category,'WS1/outlierLists/SNPFst_NarrowGenesMin2SNPs.Alps.txt',col.names = F,row.names = F,quote = F)
write.table(b2$Category,'WS1/outlierLists/SNPFst_NarrowGenesMin2SNPs.Fagaras.txt',col.names = F,row.names = F,quote = F)
write.table(c2$Category,'WS1/outlierLists/SNPFst_NarrowGenesMin2SNPs.Rodna.txt',col.names = F,row.names = F,quote = F)
write.table(d2$Category,'WS1/outlierLists/SNPFst_NarrowGenesMin2SNPs.Tatry2x.txt',col.names = F,row.names = F,quote = F)
write.table(e2$Category,'WS1/outlierLists/SNPFst_NarrowGenesMin2SNPs.Tatry4x.txt',col.names = F,row.names = F,quote = F)

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
write.table(f2$Category,'WS1/outlierLists/SNPFst_NarrowGenesMin2SNPs.AlpsH.txt',col.names = F,row.names = F,quote = F)
write.table(g2$Category,'WS1/outlierLists/SNPFst_NarrowGenesMin2SNPs.FagarasH.txt',col.names = F,row.names = F,quote = F)

write.table(tt3$Category,'WS1/outlierLists/SNPFst_arenosa_halleri_NarrowGenesMin2SNPs.txt',col.names = F,row.names = F,quote = F)
##find overlaps and annotate
ann<-fread("../../../Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20181231.txt",h=T,quote="")
mm3<-subset(ann,ann$AL %in% i2)
write.table(mm3,'WS1/outlierLists/SNPFst_arenosa_halleri_NarrowGenesMin2SNPs_ann.txt',col.names = F,row.names = F,quote = F,sep="\t")

mm3<-"AL1G10530"
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

################# 4. BAYPASS OVERLAP - WITHIN SPECIES ######################
### OVERLAP BAYPASS WITH FST QUARTETS ###
library(data.table)
setwd("/home/aa/alpine/arenosaGenome/selScans/") ##
write.table(paste("region SNPs_full SNPs_upstream SNPs_genic SNPs_missense SNPs_downstream SNPs_intergenic genes_full.txt genes_upstream genes_genic genes_missense genes_downstream",sep=""),paste("bayPass_quartetFst/BP_Fst_overlaps.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ")

for (lin in c("HalleriA","HalleriF","Tatry2x","Tatry4x","Fagaras","Rodna","Alps")) {# lin<-"Alps" c("Tatry2x","Tatry4x","Fagaras","Rodna","Alps")
b<-fread(paste("../../bayPass/",lin,"/outSNPsAll_",lin,".BPann.0.5.txt",sep=''),h=T) #
#full
f<-fread(paste("quartetOutliers/outSNPsAll_",lin,"0.95.txt",sep=''),h=T) #
#upstream
f1<-subset(f,f$ann %in% "upstream_gene_variant") 
#genic
f2<-subset(f,!f$ann %in% "intragenic_variant" & !f$ann %in% "downstream_gene_variant" & !f$ann %in% "upstream_gene_variant" & !f$ann %in% "intergenic_region") 
#missense
f3<-subset(f,f$ann %in% "missense_variant") 
#downstream
f4<-subset(f,f$ann %in% "downstream_gene_variant") 
#intergenic
f5<-subset(f, f$ID %like% "-") 
#full
colnames(f)[2]<-"start"
k<-rbind(b[,c(7,4)],f[,c(1,2)]) 
k<-k[ order(k[,1],k[,2]), ]
k<-unique(k[duplicated(k[,c(1,2)]),])
#upstream
colnames(f1)[2]<-"start"
k1<-rbind(b[,c(7,4)],f1[,c(1,2)]) 
k1<-k1[ order(k1[,1],k1[,2]), ]
k1<-unique(k1[duplicated(k1[,c(1,2)]),])
#genic
colnames(f2)[2]<-"start"
k2<-rbind(b[,c(7,4)],f2[,c(1,2)]) 
k2<-k2[ order(k2[,1],k2[,2]), ]
k2<-unique(k2[duplicated(k2[,c(1,2)]),])
#missense
colnames(f3)[2]<-"start"
k3<-rbind(b[,c(7,4)],f3[,c(1,2)]) 
k3<-k3[ order(k3[,1],k3[,2]), ]
k3<-unique(k3[duplicated(k3[,c(1,2)]),])
#downstream
colnames(f4)[2]<-"start"
k4<-rbind(b[,c(7,4)],f4[,c(1,2)]) 
k4<-k4[ order(k4[,1],k4[,2]), ]
k4<-unique(k4[duplicated(k4[,c(1,2)]),])
#intergenic
colnames(f5)[2]<-"start"
k5<-rbind(b[,c(7,4)],f5[,c(1,2)]) 
k5<-k5[ order(k5[,1],k5[,2]), ]
k5<-unique(k5[duplicated(k5[,c(1,2)]),])
### assign to gene
k$one<-1
k1$one<-1
k2$one<-1
k3$one<-1
k4$one<-1
o<-setDT(k)[,list(perGene=sum(one)),by=list(Category=ID)]
o1<-setDT(k1)[,list(perGene=sum(one)),by=list(Category=ID)]
o2<-setDT(k2)[,list(perGene=sum(one)),by=list(Category=ID)]
o3<-setDT(k3)[,list(perGene=sum(one)),by=list(Category=ID)]
o4<-setDT(k4)[,list(perGene=sum(one)),by=list(Category=ID)]
#f$one<-1
#f1$one<-1
#f2$one<-1
#f3$one<-1
#f4$one<-1
# o<-setDT(f)[,list(perGene=sum(one)),by=list(Category=ID)]
#o1<-setDT(f1)[,list(perGene=sum(one)),by=list(Category=ID)]
#o2<-setDT(f2)[,list(perGene=sum(one)),by=list(Category=ID)]
#o3<-setDT(f3)[,list(perGene=sum(one)),by=list(Category=ID)]
#o4<-setDT(f4)[,list(perGene=sum(one)),by=list(Category=ID)]
o<-subset(o,o$perGene>5)
o1<-subset(o1,o1$perGene>5)
o2<-subset(o2,o2$perGene>5)
o3<-subset(o3,o3$perGene>5)
o4<-subset(o4,o4$perGene>5)
#export
write.table(k,paste("bayPass_quartetFst/",lin,"/SNPs_full.txt",sep=""),quote = F,col.names = F,row.names = F)
write.table(k1,paste("bayPass_quartetFst/",lin,"/SNPs_upstream.txt",sep=""),quote = F,col.names = F,row.names = F)
write.table(k2,paste("bayPass_quartetFst/",lin,"/SNPs_genic.txt",sep=""),quote = F,col.names = F,row.names = F)
write.table(k3,paste("bayPass_quartetFst/",lin,"/SNPs_missense.txt",sep=""),quote = F,col.names = F,row.names = F)
write.table(k4,paste("bayPass_quartetFst/",lin,"/SNPs_downstream.txt",sep=""),quote = F,col.names = F,row.names = F)
write.table(k5,paste("bayPass_quartetFst/",lin,"/SNPs_intergenic.txt",sep=""),quote = F,col.names = F,row.names = F)
 write.table(o,paste("bayPass_quartetFst/",lin,"/genes_full.txt",sep=""),quote = F,col.names = F,row.names = F)
write.table(o1,paste("bayPass_quartetFst/",lin,"/genes_upstream.txt",sep=""),quote = F,col.names = F,row.names = F)
write.table(o2,paste("bayPass_quartetFst/",lin,"/genes_genic.txt",sep=""),quote = F,col.names = F,row.names = F)
write.table(o3,paste("bayPass_quartetFst/",lin,"/genes_missense.txt",sep=""),quote = F,col.names = F,row.names = F)
write.table(o4,paste("bayPass_quartetFst/",lin,"/genes_downstream.txt",sep=""),quote = F,col.names = F,row.names = F)
write.table(paste(lin,nrow(k),nrow(k1),nrow(k2),nrow(k3),nrow(k4),nrow(k5),nrow(o),nrow(o1),nrow(o2),nrow(o3),nrow(o4)),paste("bayPass_quartetFst/BP_Fst_overlaps.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
write.table(t(c(c(lin),c(o$Category))),paste("bayPass_quartetFst/genes_full.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
write.table(t(c(c(lin),c(o1$Category))),paste("bayPass_quartetFst/genes_upstream.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
write.table(t(c(c(lin),c(o2$Category))),paste("bayPass_quartetFst/genes_genic.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
write.table(t(c(c(lin),c(o3$Category))),paste("bayPass_quartetFst/genes_missense.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
write.table(t(c(c(lin),c(o4$Category))),paste("bayPass_quartetFst/genes_downstream.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  write.table(t(c(c(lin),paste(k$ID,k$start,sep="."))),paste("bayPass_quartetFst/SNPs_full.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
write.table(t(c(c(lin),paste(k1$ID,k1$start,sep="."))),paste("bayPass_quartetFst/SNPs_upstream.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
write.table(t(c(c(lin),paste(k2$ID,k2$start,sep="."))),paste("bayPass_quartetFst/SNPs_genic.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
write.table(t(c(c(lin),paste(k3$ID,k3$start,sep="."))),paste("bayPass_quartetFst/SNPs_missense.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
write.table(t(c(c(lin),paste(k4$ID,k4$start,sep="."))),paste("bayPass_quartetFst/SNPs_downstream.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
write.table(t(c(c(lin),paste(k5$ID,k5$start,sep="."))),paste("bayPass_quartetFst/SNPs_intergenic.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
#write.table(t(c(c(lin),paste(f$ID,f$i.start,sep="."))),paste("bayPass_quartetFst/SNPs_full.txt",sep=""),quote = F,col.names = F,row.names #= F,sep = " ",append = T)
#write.table(t(c(c(lin),paste(f1$ID,f1$i.start,sep="."))),paste("bayPass_quartetFst/SNPs_upstream.txt",sep=""),quote = F,col.names = F,row#.names = F,sep = " ",append = T)
#write.table(t(c(c(lin),paste(f2$ID,f2$i.start,sep="."))),paste("bayPass_quartetFst/SNPs_genic.txt",sep=""),quote = F,col.names = F,row#.names = F,sep = " ",append = T)
#write.table(t(c(c(lin),paste(f3$ID,f3$i.start,sep="."))),paste("bayPass_quartetFst/SNPs_missense.txt",sep=""),quote = F,col.names = F,row#.names = F,sep = " ",append = T)
#write.table(t(c(c(lin),paste(f4$ID,f4$i.start,sep="."))),paste("bayPass_quartetFst/SNPs_downstream.txt",sep=""),quote = F,col.names = F#,row.names = F,sep = " ",append = T)
#write.table(t(c(c(lin),paste(f5$ID,f5$i.start,sep="."))),paste("bayPass_quartetFst/SNPs_intergenic.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
}
################# 4. BAYPASS OVERLAP - ARENOSA  (HALLERI CAN BE DONE USING THE ABOVE) #
### OVERLAP BAYPASS WITH FST QUARTETS ###
library(data.table)
setwd("/home/aa/alpine/arenosaGenome/selScans/")
write.table(paste("region SNPs_full SNPs_upstream SNPs_genic SNPs_missense SNPs_downstream SNPs_intergenic genes_full genes_upstream genes_genic genes_missense genes_downstream",sep=""),paste("bayPass_quartetFst/BP_Fst_overlaps_BetweenSpecies.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ")
lin<-"Arenosa"
  b<-fread("/home/aa/alpine/bayPass/All/candCOREstep2.txt",h=F)
  b$scaff<-substr(b$V1,10,10)
  b$ID<-"AA"
  #full
  f<-fread(paste("quartetOutliers/outSNPsAll_",lin,"0.99.txt",sep=''),h=T)
  f$scaff<-substr(f$ID,3,3)
  #upstream
  f1<-subset(f,f$ann %in% "upstream_gene_variant") 
  #genic
  f2<-subset(f,!f$ann %in% "intragenic_variant" & !f$ann %in% "downstream_gene_variant" & !f$ann %in% "upstream_gene_variant" & !f$ann %in% "intergenic_region") 
  #missense
  f3<-subset(f,f$ann %in% "missense_variant") 
  #downstream
  f4<-subset(f,f$ann %in% "downstream_gene_variant") 
  #intergenic
  f5<-subset(f, f$ID %like% "-") 
  #full
  colnames(f)[2]<-"start"
  colnames(b)[2]<-"start"
  k<-rbind(b[,c(4,2,5)],f[,c(8,2,1)]) 
  k<-k[ order(k[,1],k[,2],k[,3]), ]
  k<-unique(k[duplicated(k[,c(1,2)]),])
  #upstream
  colnames(f1)[2]<-"start"
  k1<-rbind(b[,c(4,2,5)],f1[,c(8,2,1)]) 
  k1<-k1[ order(k1[,1],k1[,2],k1[,3]), ]
  k1<-unique(k1[duplicated(k1[,c(1,2)]),])
  #genic
  colnames(f2)[2]<-"start"
  k2<-rbind(b[,c(4,2,5)],f2[,c(8,2,1)]) 
  k2<-k2[ order(k2[,1],k2[,2],k2[,3]), ]
  k2<-unique(k2[duplicated(k2[,c(1,2)]),])
  #missense
  colnames(f3)[2]<-"start"
  k3<-rbind(b[,c(4,2,5)],f3[,c(8,2,1)]) 
  k3<-k3[ order(k3[,1],k3[,2],k3[,3]), ]
  k3<-unique(k3[duplicated(k3[,c(1,2)]),])
  #downstream
  colnames(f4)[2]<-"start"
  k4<-rbind(b[,c(4,2,5)],f4[,c(8,2,1)]) 
  k4<-k4[ order(k4[,1],k4[,2],k4[,3]), ]
  k4<-unique(k4[duplicated(k4[,c(1,2)]),])
  #intergenic
  colnames(f5)[2]<-"start"
  k5<-rbind(b[,c(4,2,5)],f5[,c(8,2,1)]) 
  k5<-k5[ order(k5[,1],k5[,2],k5[,3]), ]
  k5<-unique(k5[duplicated(k5[,c(1,2)]),])
  ### assign to gene
  k$one<-1
  k1$one<-1
  k2$one<-1
  k3$one<-1
  k4$one<-1
  o<-setDT(k)[,list(perGene=sum(one)),by=list(Category=ID)]
  o1<-setDT(k1)[,list(perGene=sum(one)),by=list(Category=ID)]
  o2<-setDT(k2)[,list(perGene=sum(one)),by=list(Category=ID)]
  o3<-setDT(k3)[,list(perGene=sum(one)),by=list(Category=ID)]
  o4<-setDT(k4)[,list(perGene=sum(one)),by=list(Category=ID)]
  o<-subset(o,o$perGene>1)
  o1<-subset(o1,o1$perGene>1)
  o2<-subset(o2,o2$perGene>1)
  o3<-subset(o3,o3$perGene>1)
  o4<-subset(o4,o4$perGene>1)
  #export
  write.table(k,paste("bayPass_quartetFst/",lin,"/SNPs_full.txt",sep=""),quote = F,col.names = F,row.names = F)
  write.table(k1,paste("bayPass_quartetFst/",lin,"/SNPs_upstream.txt",sep=""),quote = F,col.names = F,row.names = F)
  write.table(k2,paste("bayPass_quartetFst/",lin,"/SNPs_genic.txt",sep=""),quote = F,col.names = F,row.names = F)
  write.table(k3,paste("bayPass_quartetFst/",lin,"/SNPs_missense.txt",sep=""),quote = F,col.names = F,row.names = F)
  write.table(k4,paste("bayPass_quartetFst/",lin,"/SNPs_downstream.txt",sep=""),quote = F,col.names = F,row.names = F)
  write.table(k5,paste("bayPass_quartetFst/",lin,"/SNPs_intergenic.txt",sep=""),quote = F,col.names = F,row.names = F)
  write.table(o,paste("bayPass_quartetFst/",lin,"/genes_full.txt",sep=""),quote = F,col.names = F,row.names = F)
  write.table(o1,paste("bayPass_quartetFst/",lin,"/genes_upstream.txt",sep=""),quote = F,col.names = F,row.names = F)
  write.table(o2,paste("bayPass_quartetFst/",lin,"/genes_genic.txt",sep=""),quote = F,col.names = F,row.names = F)
  write.table(o3,paste("bayPass_quartetFst/",lin,"/genes_missense.txt",sep=""),quote = F,col.names = F,row.names = F)
  write.table(o4,paste("bayPass_quartetFst/",lin,"/genes_downstream.txt",sep=""),quote = F,col.names = F,row.names = F)
  write.table(paste(lin,nrow(k),nrow(k1),nrow(k2),nrow(k3),nrow(k4),nrow(k5),nrow(o),nrow(o1),nrow(o2),nrow(o3),nrow(o4)),paste("bayPass_quartetFst/BP_Fst_overlaps.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  write.table(t(c(c(lin),c(o$Category))),paste("bayPass_quartetFst/genes_full.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  write.table(t(c(c(lin),c(o1$Category))),paste("bayPass_quartetFst/genes_upstream.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  write.table(t(c(c(lin),c(o2$Category))),paste("bayPass_quartetFst/genes_genic.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  write.table(t(c(c(lin),c(o3$Category))),paste("bayPass_quartetFst/genes_missense.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  write.table(t(c(c(lin),c(o4$Category))),paste("bayPass_quartetFst/genes_downstream.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  write.table(t(c(c(lin),paste(k$ID,k$start,sep="."))),paste("bayPass_quartetFst/SNPs_full.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  write.table(t(c(c(lin),paste(k1$ID,k1$start,sep="."))),paste("bayPass_quartetFst/SNPs_upstream.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  write.table(t(c(c(lin),paste(k2$ID,k2$start,sep="."))),paste("bayPass_quartetFst/SNPs_genic.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  write.table(t(c(c(lin),paste(k3$ID,k3$start,sep="."))),paste("bayPass_quartetFst/SNPs_missense.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  write.table(t(c(c(lin),paste(k4$ID,k4$start,sep="."))),paste("bayPass_quartetFst/SNPs_downstream.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  write.table(t(c(c(lin),paste(k5$ID,k5$start,sep="."))),paste("bayPass_quartetFst/SNPs_intergenic.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  
############ OVERLAP SUPERTEST ############
#https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html
library("SuperExactTest")
setwd("/home/aa/alpine/arenosaGenome/selScans/bayPass_quartetFst/SuperTestOutput/")
for (dat in c("SNPs_genic","SNPs_missense","SNPs_upstream")) { # "SNPs_intergenic" dat="SNPs_full" "SNPs_downstream","SNPs_full",
    d <- scan(paste(dat,".txt",sep=""), what="", sep="\n")
    d <- strsplit(d, "[[:space:]]+")
    #  names(d) <- sapply(d, `[[`, 1)
    names(d) <- c("HN","HF","VT","ZT","FG","RD","NT")
    d <- lapply(d, `[`, -1)
    total=11093985
    res=supertest(d, n=total)
    res$overlap.sizes
 #   pdf(paste(dat,".pdf",sep=""),width = 12,height = 8) ###within species comment
    png(paste(dat,".png",sep=""),width = 1200,height = 960,pointsize = 24)
    plot(res, Layout="landscape", sort.by="size", margin=c(0.5,5,1,2),keep.empty.intersections=F,min.intersection.size = 11, degree=2:7)
    #plot(res, sort.by="size", margin=c(2,2,2,2), color.scale.pos=c(0.85,1), legend.pos=c(0.9,0.1),keep.empty.intersections=F)
    #plot(res, sort.by="size", margin=c(2,2,2,2), color.scale.pos=c(0.85,1), legend.pos=c(0.9,0.1),keep.empty.intersections=F)
    dev.off()
    write.csv(summary(res)$Table, file=paste(dat,".csv",sep=""), row.names=FALSE)}

for (dat in c("genes_genic","genes_missense","genes_upstream")) { # dat="genes_genic" # "genes_downstream" "genes_full",
  d <- scan(paste(dat,".txt",sep=""), what="", sep="\n")
  d <- strsplit(d, "[[:space:]]+")
  names(d) <- c("HN","HF","VT","ZT","FG","RD","NT")
  d <- lapply(d, `[`, -1)
  total=30000
  res=supertest(d, n=total)
  res$overlap.sizes
  png(paste(dat,".png",sep=""),width = 1200,height = 960,pointsize = 24)
  plot(res, Layout="landscape", sort.by="size", margin=c(0.5,5,1,2),keep.empty.intersections=F,min.intersection.size = 2, degree=2:7) 
  #plot(res, Layout="landscape", degree=2:7, sort.by="size", margin=c(0.5,5,1,2),keep.empty.intersections=F)
 # plot(res, sort.by="size", margin=c(2,2,2,2), color.scale.pos=c(0.85,1), legend.pos=c(0.9,0.1),keep.empty.intersections=F)
  #plot(res, sort.by="size", margin=c(2,2,2,2), color.scale.pos=c(0.85,1), legend.pos=c(0.9,0.1),keep.empty.intersections=F)
  dev.off()
  write.csv(summary(res)$Table, file=paste(dat,".csv",sep=""), row.names=FALSE)}

############ OVERLAP SUPERTEST ARENOSA HALLERI ############
library("SuperExactTest")
setwd("/home/aa/alpine/arenosaGenome/selScans/bayPass_quartetFst/withinHalleri/") #betweenSpecies
for (dat in c("SNPs_downstream","SNPs_full","SNPs_genic","SNPs_missense","SNPs_upstream")) { # dat="SNPs_genic" # "genes_downstream","SNPs_intergenic"
  d <- scan(paste(dat,".txt",sep=""), what="", sep="\n")
  d <- strsplit(d, "[[:space:]]+")
  #  names(d) <- sapply(d, `[[`, 1)
  names(d) <- c("HN","HF") #c("AH","AA")
  d <- lapply(d, `[`, -1)
  total=6713051
  res=supertest(d, n=total)
  res$overlap.sizes
  write.csv(summary(res)$Table, file=paste(dat,".csv",sep=""), row.names=FALSE)}
for (dat in c("genes_full","genes_genic","genes_missense","genes_upstream")) { # dat="SNPs_genic" # "genes_downstream","SNPs_intergenic"
  d <- scan(paste(dat,".txt",sep=""), what="", sep="\n")
  d <- strsplit(d, "[[:space:]]+")
  #  names(d) <- sapply(d, `[[`, 1)
  names(d) <- c("HN","HF") #c("AH","AA")
  d <- lapply(d, `[`, -1)
  total=30000
  res=supertest(d, n=total)
  res$overlap.sizes
  write.csv(summary(res)$Table, file=paste(dat,".csv",sep=""), row.names=FALSE)}


############ SNP2GO BEFORE START CHANGE ALL ARENOSA TO OTHERS FORMAT !!! ##############
#  install.packages("hash")
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("GenomeInfoDb")
#  biocLite("goProfiles")
#  biocLite("GenomicRanges")
#  install.packages('devtools')
#  install.packages(pkgs="/home/aa/Desktop/programs/SNP2GO_1.0.5.tar.gz", type="source",repos = NULL,Ncpus = 2)
setwd("/home/aa/alpine/arenosaGenome/selScans/SNP2GO/")
library(SNP2GO)
system("cut -f2,4 ../BPM/HIGLOW_WS1000_MS1_BPM.txt | awk 'NR == 1 || NR % 30 == 0' > nonsel.txt") 
system("cut -f2,4 ../BPM/HIGLOWh_WS1000_MS1_BPM.txt | awk 'NR == 1 || NR % 34 == 0' > nonselh.txt") #Halleri

l<- c("Tatry2x","Tatry4x","Fagaras","Rodna","Alps") #  c("HalleriA","HalleriF")
for (lin in l){ # lin = "Alps"  ### "Halleri"
  s<-read.table(paste("../bayPass_quartetFst/",lin,"/SNPs_full.txt",sep=""),h=F) 
  snps<-as.data.frame(paste("scaffold_",substr(s$V1,3,3),sep=""))
  snps$end <- as.numeric(s$V2)
  cand <- GRanges(seqnames=snps[,1],ranges=IRanges(snps[,2],snps[,2]))
  noncand <- read.table("nonsel.txt",h=F)
  noncand[,2] <- as.numeric(noncand[,2])
  noncand <- GRanges(seqnames=noncand[,1],ranges=IRanges(noncand[,2],noncand[,2]))
  y <- snp2go(gtf="/home/aa/Desktop/references/lyrataV2/LyV2.gtf", goFile="/home/aa/Desktop/references/lyrataV2/mart_Alyrata_V2", candidateSNPs=cand, noncandidateSNPs=noncand, FDR=0.05, runs=1000, extension=0, min.regions=1)
  write.table(file=paste("../bayPass_quartetFst/",lin,"/SNPs_full_snp2go.tsv",sep=""),y$enriched,sep="\t",row.names=F)
  }
# revigo - GOs, g, tiny, save 2 r scripts and .csv
# ggsave("revigo-plot.pdf",width = 12,height = 8);
# work further with bayPass_quartetFst/lin/revigo.csv - 0 - not discarded terms

##### SUMMARIZE GO PATHWAYS WITHIN arenosa/halleri #####
setwd("/home/aa/alpine/arenosaGenome/selScans/bayPass_quartetFst/")
write.table(paste("region go_all go_revigo",sep=""),paste("GO_numbers.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ")
for (lin in c("HalleriA","HalleriF","Tatry2x","Tatry4x","Fagaras","Rodna","Alps")) { #lin<-"Tatry2x"  c("Tatry2x","Tatry4x","Fagaras","Rodna","Alps")
b<-read.csv(paste(lin,"/REVIGO.csv",sep=''),h=T)
b1<-subset(b,b$eliminated %in% 0)
write.table(paste(lin,nrow(b),nrow(b1)),paste("GO_numbers.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  write.table(t(c(c(lin),as.character(b$term_ID))),paste("go_all.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  write.table(t(c(c(lin),as.character(b1$term_ID))),paste("go_revigo.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)}
##### SUMMARIZE GO PATHWAYS ARENOSA-HALLERI #####
setwd("/home/aa/alpine/arenosaGenome/selScans/bayPass_quartetFst/")
write.table(paste("region go_all go_revigo",sep=""),paste("betweenSpecies/GO_numbers.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ")
for (lin in c("Arenosa","Halleri")) {   #   lin<-"Tatry2x"
  b<-read.csv(paste(lin,"/REVIGO.csv",sep=''),h=T)
  b1<-subset(b,b$eliminated %in% 0)
  write.table(paste(lin,nrow(b),nrow(b1)),paste("betweenSpecies/GO_numbers.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  write.table(t(c(c(lin),as.character(b$term_ID))),paste("betweenSpecies/go_all.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
  write.table(t(c(c(lin),as.character(b1$term_ID))),paste("betweenSpecies/go_revigo.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)}
############ OVERLAP PATHWAYS SUPERTEST ############
library("SuperExactTest")
setwd("/home/aa/alpine/arenosaGenome/selScans/bayPass_quartetFst/SuperTestOutput/") # withinArenosa/ betweenSpecies
for (dat in c("go_revigo","go_all")) { 
  d <- scan(paste(dat,".txt",sep=""), what="", sep="\n")
  d <- strsplit(d, "[[:space:]]+")
  names(d) <- c("HN","HF","VT","ZT","FG","RD","NT") # c("VT","ZT","FG","RD","NT") c("AA","AH")
  d <- lapply(d, `[`, -1)
  total=6063
  res=supertest(d, n=total)
  res$overlap.sizes
  png(paste(dat,"1.png",sep=""),width = 4800,height = 960,pointsize = 24)
  plot(res, Layout="landscape", sort.by="size", margin=c(0.5,5,1,2),keep.empty.intersections=F,min.intersection.size = 11, degree=2:7,show.overlap.size = T)
#  plot(res, sort.by="size", margin=c(2,2,2,2), color.scale.pos=c(0.85,1), legend.pos=c(0.9,0.1),keep.empty.intersections=F)
  dev.off()
  write.csv(summary(res)$Table, file=paste(dat,".csv",sep=""), row.names=FALSE)}


### CALCULATE PROPORTIONS ###
##SNPs
setwd("/home/aa/alpine/arenosaGenome/selScans/bayPass_quartetFst/")
library(data.table)
names<-c("genes_full","genes_genic","genes_missense","genes_upstream","SNPs_downstream","SNPs_full","SNPs_genic","SNPs_missense","SNPs_upstream","go_all","go_revigo")
pairs<-c("VTZT 0.0774 1.9678","NTZT 0.1240 1.4800","FGRD 0.1336 1.7648","RDZT 0.1409 1.4224","VTNT 0.1420 2.0544","VTRD 0.1422 1.9968","FGZT 0.1428 1.7358","VTFG 0.1550 2.3102","NTFG 0.1553 1.8223","NTRD 0.1644 1.5090","HFHN 0.2040 1.8488","RDHF 0.3951 1.5420","FGHF 0.4028 1.6670","NTHF 0.4060 1.6689","VTHF 0.4088 2.0948","RDHN 0.4178 1.6669","VTHN 0.4280 2.2198","FGHN 0.4289 1.7920","NTHN 0.4328 1.7938","ZTHF 0.4339 1.5106","ZTHN 0.4579 1.6356") ### divergence
#mat<-matrix(nrow = 10,ncol = 9)
for (name in names) { # name<-"SNPs_full"
  s<-read.csv(paste(name,".csv",sep=""),h=T) 
#sa<-read.csv(paste("withinArenosa/",name,".csv",sep=""),h=T)
#sh<-read.csv(paste("withinHalleri/",name,".csv",sep=""),h=T)
#sb<-read.csv(paste("betweenSpecies/",name,".csv",sep=""),h=T)
#s<-rbind(sa,sh,sb)
for (pair in pairs) { # pair="VTFG"
  p1<-substr(pair,1,2)
  p2<-substr(pair,3,4)
t1<-(subset(s,s$Intersections %in% p1))[1,3]
t2<-(subset(s,s$Intersections %in% p2))[1,3]
par<-subset(s,s$Intersections %like% p1 & s$Intersections %like% p2 & s$Degree %in% "2")[1,3]
prop<-par/(t1+t2)
y<-which(names %in% name)
x<-which(pairs %in% pair)
#mat[x,y]<-prop
#write.table(paste(pair,t1,t2,par,prop),paste(name,".proportions.txt",sep=""),quote = F,col.names = F,row.names = F,sep = " ",append = T)
write.table(paste(pair,name,prop),paste("all.ggplot2.proportions.txt",sep=""),quote = F,row.names = F,sep = " ",append = T,col.names = F)
}}
#write.table(mat,"all.proportions.txt",col.names = names,row.names = pairs,quote =F)




###  ANALYZE THE EFFECT OF DIVERGENCE ON PARALLELS ###
library(ade4)
setwd("/home/aa/alpine/arenosaGenome/selScans/bayPass_quartetFst/SuperTestOutput/")
#dfsthh <- read.table("fst.dist.txt",h=T)[1,2]
dfst <- as.dist(read.table("fst.dist.txt",h=T),diag = T,upper = T)
#dfst <- as.dist(read.table("fst.dist.txt",h=T)[3:7,3:7],diag = T,upper = T)
data<-c("SNPs_missense","SNPs_upstream") #c("genes_genic","SNPs_genic","go_all")
for (dat in data) { # dat<-"SNPs_genic"
d <- scan(paste(dat,".txt",sep=""), what="", sep="\n")
d <- strsplit(d, "[[:space:]]+")
names(d) <- c("HN","HF","VT","ZT","FG","RD","NT") # c("VT","ZT","FG","RD","NT") c("AA","AH")
d <- lapply(d, `[`, -1)
my_genes <- unique(unlist(d))
my_function <- function(x){
  is.element(my_genes,x)
}
df <- as.data.frame(lapply(d, my_function))
rownames(df) <- my_genes
df[df=="TRUE"]<- 1
#   df$sum <- apply(df,1,sum)
#   i2<-rownames(subset(df,df$sum > 2))
df1<-t(df)
#dist.eucl <- dist(df, method = "binary")
#round(as.matrix(dist.eucl)[1:7, 1:7], 4)
#for (i in 1:10) {
#  d <- dist.binary(df, method = i)
#  write.table(round(as.matrix(d)[1:7, 1:7], 4),paste(attr(d, "method"),".txt",sep=""),quote = F)}
#for (i in 1:10) {
#  d <- dist.binary(df, method = i)
#  write.table(round(as.matrix(d)[1:7, 1:7], 4),"all.txt",quote = F,append = T)}

dgen <- dist.binary(df1, method = 1, diag = T, upper = T)
#dgen<-as.matrix(dgen)[3:7,3:7]
#dgenhh<-as.matrix(dgen)[1,2]
#dgen<-as.dist(dgen)

#Mantel test
dist <- mantel.randtest(dfst,dgen)
dist
pdf (paste("corr_dgen_dfst_",dat,".pdf",sep=""), width=14, height=7)
plot(dfst,dgen, pch=20,cex=.5)
abline(lm(dgen~dfst))
dev.off()
#   #dbRDA #works but probably nonsense
#   library(vegan)
#   setwd("/home/aa/alpine/arenosaGenome/selScans/bayPass_quartetFst/SuperTestOutput/")
#   df<-as.data.frame(df1)
#   fst <- read.table("fst.dist.txt",h=T)
#   rankindex(fst, df, indices = c("euc", "man", "gow","bra", "kul"),stepacross= FALSE, method =   "spearman")
#  dbRDA=capscale(formula = dgen ~ HN+HF+VT+ZT+FG+RD+NT, data = fst)
#  plot(dbRDA)
#  anova(dbRDA)
#Proctustes
library(vegan)
proc <- procrustes(X = dfst, Y = dgen)
summary(proc)
pdf (paste("Procrustes_",dat,".pdf",sep=""), width=14, height=7)
par(mfrow=c(1,2))
plot(proc,type = "text")
plot(proc, kind=2,xaxt="n")
axis(side = 1, labels=c("HN","HF","VT","ZT","FG","RD","NT"),at = seq(1,7,1))
dev.off()
residuals(proc)
fitted(proc, truemean = TRUE)
protest(X = dfst, Y = dgen, scores = "sites")
#PCoA
library(ape)
a<-pcoa(dgen, correction="none", rn=NULL)
b<-pcoa(dfst, correction="none", rn=NULL)
# The proportion of variances explained is in its element Relative_eig, total "variance" is returned in element trace
 a$values
# a$trace
# a$vectors
# b$values
# b$trace
# b$vectors
pdf (paste("PCoA_dgen_dfst_",dat,".pdf",sep=""), width=14, height=7)
par(mfrow=c(1,2))
biplot(a, Y=NULL, plot.axes = c(1,2), dir.axis1=1, dir.axis2=1, rn=NULL, main="dgen")
biplot(a, Y=NULL, plot.axes = c(1,3), dir.axis1=1, dir.axis2=1, rn=NULL, main=NULL)
par(mfrow=c(1,2))
biplot(b, Y=NULL, plot.axes = c(1,2), dir.axis1=1, dir.axis2=1, rn=NULL, main="dfst")
biplot(b, Y=NULL, plot.axes = c(1,3), dir.axis1=1, dir.axis2=1, rn=NULL, main=NULL)
dev.off()
# Model II regression - on PCA axis 1
library(lmodel2)
res11 <- lmodel2(a$vectors[,1] ~ b$vectors[,1], nperm=999)
#res11
pdf (paste("model.II.regr_dgen_dfst_PCAaxis1_",dat,".pdf",sep=""), width=7, height=7)
plot(res11) 
dev.off()
#res22 <- lmodel2(a$vectors[,2] ~ b$vectors[,2], nperm=999)
#res22
#plot(res22) 

# Model II regression - on distances, with jackknife
library(bootstrap)
df<-c(dfst)
dg<-c(dgen)
#df<-c(dfst,dfsthh)
#dg<-c(dgen,dgenhh)
dftot<-as.data.frame(cbind(df,dg))
res <- lmodel2(dftot$dg ~ dftot$df, nperm=999)
#res
#res$regression.results$Slope[1]
#res[[3]][1,3]
#res$regression.results$`P-perm (1-tailed)`[1]
xdata <-dftot
n <- 11
theta <- function(x,xdata){ lmodel2(xdata[x,2]~xdata[x,1])[[3]][1,3] }
results <- jackknife(1:n,theta,xdata)
mean<-mean(results$jack.values)
error <- qnorm(0.975)*sd(results$jack.values)/sqrt(11)
left <- mean(results$jack.values)-error
right <- mean(results$jack.values)+error
#print(res)
print(dat)
print(mean)
print(left)
print(right)

pdf (paste("model.II.regr_dgen_dfst_",dat,".pdf",sep=""), width=7, height=7,pointsize = 24)
plot(res,xlab="Divergence (Fst)",ylab="distance between parallels",lty=20,main = NA,ylim=c(0.98,1)) 
dev.off()
}


### PLOT AND TEST THE DIVERGENCE - cannot use, points in the regression are not independent ###
setwd("/home/aa/alpine/arenosaGenome/selScans/bayPass_quartetFst/")
all<-read.table("all.ggplot2.proportions.txt",h=F)
colnames(all)<- c("contr","div","adapt","pos","prop")
all$type<-substr(all$pos,1,3)
all$prop1<-(all$prop*100)+1  #can I do this?
library(ggplot2)
library(lme4)
library(ggpubr)
library(data.table)

#1. test the effect of divergence and all possible positions on parallelism.. DOES NOT MAKE MUCH SENSE 
# ggplot(all, aes(x=div, y=prop, color=as.factor(pos))) + geom_point() + geom_smooth(method=lm, se=FALSE, # fullrange=TRUE) + labs(title="All datasets", x="Divergence (neutral Fst)", y = "Proportion of parallels"# )
# alla<-subset(all,!all$pos %like% "go")
# hist(log(alla$prop1))
# hist(alla$prop)
# lm<-aov(log(prop1)~log(div)*as.factor(pos),data=alla)
# summary(lm)
# lma<-aov(log(prop1)~log(adapt)*as.factor(pos),data=alla)
# summary(lma)
# lm1<-aov(prop~as.factor(pos)*div,data=all)
# summary(lm1)
# lm2<-aov(prop~as.factor(pos)*adapt,data=all)
# summary(lm2)
# # glm
# fit <- glm(prop ~ div * pos,data=all,family=binomial())
# summary(fit)
# lmAll1<-glmer(prop ~ div + (type|pos), data = alla,family = "quasibinomial")
# summary(lmAll1)
# lmAll1<-glmer(prop ~ pos + (type|pos), data = alla,family = "binomial")
# summary(lmAll1)


###### 2. test the effect of CODING VS. REGULATORY position ###### 
all1<-subset(all,all$pos %like% "SNP" & !all$pos %like% "downstream" & !all$pos %like% "full" & !all$pos %like% "genic")
all1w<-subset(all1, !contr %in% "VTZT" )
all1ww<-subset(all1,!contr %like% "H" & !contr %in% "VTZT" | contr %in% "HFHN")
all1wt<-subset(all1,!contr %like% "H" | contr %in% "HFHN")
all1be<-subset(all1,contr %like% "H" & !contr %in% "HFHN")

pdf("coding_cisReg.pdf",height = 7,width = 7,pointsize = 12)
ggplot(all1, aes(x=div, y=prop,color=factor(pos, labels=c("Missense","Upstream")), shape=factor(pos,labels=c("Missense","Upstream")))) + geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + labs(title="Proportion of parallels based on their position in the genome - all incl. VTZT", x="Divergence (neutral Fst)", y = "Proportion of parallels",color="",shape="")
ggplot(all1w, aes(x=div, y=prop,color=factor(pos, labels=c("Missense","Upstream")), shape=factor(pos,labels=c("Missense","Upstream")))) + geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + labs(title="Proportion of parallels based on their position in the genome - all without VTZT", x="Divergence (neutral Fst)", y = "Proportion of parallels",color="",shape="")
ggplot(all1wt, aes(x=div, y=prop,color=factor(pos, labels=c("Missense","Upstream")), shape=factor(pos,labels=c("Missense","Upstream")))) + geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + labs(title="Proportion of parallels based on their position in the genome - within incl. VTZT", x="Divergence (neutral Fst)", y = "Proportion of parallels",color="",shape="")
ggplot(all1ww, aes(x=div, y=prop,color=factor(pos, labels=c("Missense","Upstream")), shape=factor(pos,labels=c("Missense","Upstream")))) + geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + labs(title="Proportion of parallels based on their position in the genome - within without VTZT", x="Divergence (neutral Fst)", y = "Proportion of parallels",color="",shape="")
ggplot(all1be, aes(x=div, y=prop,color=factor(pos, labels=c("Missense","Upstream")), shape=factor(pos,labels=c("Missense","Upstream")))) + geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + labs(title="Proportion of parallels based on their position in the genome - between", x="Divergence (neutral Fst)", y = "Proportion of parallels",color="",shape="")
dev.off()

miw<-subset(all1wt,all1wt$pos %like% "missense")
upw<-subset(all1wt,all1wt$pos %like% "upstream")
mib<-subset(all1be,all1be$pos %like% "missense")
upb<-subset(all1be,all1be$pos %like% "upstream")
pdf("coding_cisReg_boxplot.pdf",width = 8,height = 10,pointsize = 24)
boxplot(cbind(miw[,5],upw[,5],mib[,5],upb[,5]),names=c("missensePops","upstreamPops","missenseSpecies","upstreamSpecies"),col = c(rep("gold",2),rep("forestgreen",2)),las=2)
dev.off()
wilcox.test(x = miw[,5],y = upw[,5],paired = T)
wilcox.test(x = mib[,5],y = upb[,5],paired = T)

pdf("coding_cisReg_plot.pdf",width = 14,height = 11)
par(mfrow=c(2,3))
plot(aov(prop~as.factor(pos)*div,data=all1w),main = "linear, without VTZT")
hist(all1w$prop)
hist(all1w$div)
plot(aov(log(prop1)~as.factor(pos)*log(div),data=all1w),main = "logaritmic, without VTZT")
hist(log(all1w$prop1))
hist(log(all1w$div))
plot(aov(sqrt(prop1)~as.factor(pos)*sqrt(div),data=all1w),main = "sqrt, without VTZT")
hist(sqrt(all1w$prop1))
hist(sqrt(all1w$div))
dev.off()

lm4<-aov(prop~as.factor(pos)*div,data=all1w)
lm4<-lm(prop~as.factor(pos)*div,data=all1w)
lm4<-aov(log(prop1)~as.factor(pos)*log(div),data=all1w) ##BETTER LOGARITMIC, BUT STILL NOT OPTIMAL..
summary(lm4)
af4<- anova(lm4)
af4ss <- af4$"Sum Sq"
print(cbind(af4,PctExp=af4ss/sum(af4ss)*100))

###### 3. parallelism by pathway, by gene or by SNP and does it change across the divergence continuum?
all2f<-subset(all, all$pos %like% "SNPs_full" | all$pos %like% "genes_full"| all$pos %like% "go_all") 
#all2g<-subset(all, all$pos %like% "SNPs_genic" | all$pos %like% "genes_genic"| all$pos %like% "go_all")
all2w <-subset(all2f,!contr %in% "VTZT" )
all2ww<-subset(all2f,!contr %like% "H" & !contr %in% "VTZT" | contr %in% "HFHN")
all2wt<-subset(all2f,!contr %like% "H" | contr %in% "HFHN")
all2be<-subset(all2f,contr %like% "H" & !contr %in% "HFHN")
#boxplot
sn<-subset(all2f,type %like% "SNP")
ge<-subset(all2f,type %like% "gen")
pa<-subset(all2f,pos %like% "go_")
pdf("snp_gene_funcion.pdf",height = 7,width = 7,pointsize = 24)
boxplot(cbind(log(sn[,7]),log(ge[,7]),log(pa[,7])),names=c("by SNP","by gene","by function"),ylab="log (Proportion of parallels)",col="grey")
boxplot(cbind(sn[,5],ge[,5],pa[,5]),names=c("by SNP","by gene","by function"),ylab="Proportion of parallels",col="grey")
dev.off()
wilcox.test(x = sn[,5],y = ge[,5],paired = T)
wilcox.test(x = pa[,5],y = ge[,5],paired = T)
wilcox.test(x = sn[,5],y = pa[,5],paired = T)
# across divergence
pdf("div_snp_gene_function.pdf",height = 7,width = 13,pointsize = 6)
a<-ggplot(all2f, aes(x=div, y=prop,color=factor(type, labels=c("By gene","By function","By SNP")), shape=factor(type, labels=c("By gene","By function","By SNP")))) + geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + labs(title="Proportion of parallels by SNP, by gene and by function - full data", x="Divergence (neutral Fst)", y = "Proportion of parallels",color="",shape="")
b<-ggplot(all2f, aes(x=div, y=log(prop1),color=factor(type, labels=c("By gene","By function","By SNP")), shape=factor(type, labels=c("By gene","By function","By SNP")))) + geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + labs(title="Proportion of parallels by SNP, by gene and by function - full data", x="Divergence (neutral Fst)", y = "log (Proportion of parallels)",color="",shape="")
ggarrange(a, b,labels = c("A", "B"), ncol = 2, nrow = 1)
c<-ggplot(all2w, aes(x=div, y=prop,color=factor(type, labels=c("By gene","By function","By SNP")), shape=factor(type, labels=c("By gene","By function","By SNP")))) + geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + labs(title="Proportion of parallels by SNP, by gene and by function - without VTZT", x="Divergence (neutral Fst)", y = "Proportion of parallels",color="",shape="")
d<-ggplot(all2w, aes(x=div, y=log(prop1),color=factor(type, labels=c("By gene","By function","By SNP")), shape=factor(type, labels=c("By gene","By function","By SNP")))) + geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + labs(title="Proportion of parallels by SNP, by gene and by function - without VTZT", x="Divergence (neutral Fst)", y = "log (Proportion of parallels)",color="",shape="")
ggarrange(c, d,labels = c("C", "D"), ncol = 2, nrow = 1)
dev.off()
#test
hist(all2a$prop1)
hist(log(all2a$prop1))
pdf("div_snp_gene_function_plot.pdf",width = 11,height = 11)
par(mfrow=c(2,3))
plot(aov(prop~as.factor(type)*div,data=all2))
plot(aov(log(prop1)~as.factor(pos)*log(div),data=all2a))
dev.off()

pdf("div_snp_gene_function_plot.pdf",width = 14, height = 11)
par(mfrow=c(2,3))
plot(aov(prop~as.factor(pos)*div,data=all2w),main = "linear, without VTZT")
hist(all2w$prop)
hist(all2w$div)
plot(aov(log(prop1)~as.factor(pos)*log(div),data=all2w),main = "logaritmic, without VTZT")
hist(log(all2w$prop1))
hist(log(all2w$div))
plot(aov(sqrt(prop1)~as.factor(pos)*sqrt(div),data=all2w),main = "sqrt, without VTZT")
hist(sqrt(all2w$prop1))
hist(sqrt(all2w$div))
dev.off()

lm5<-aov(log(prop1)~as.factor(pos)*div,data=all2w)
summary(lm5)
af5<- anova(lm5)
af5ss <- af5$"Sum Sq"
print(cbind(af5,PctExp=af5ss/sum(af5ss)*100))

lm5<-aov(prop~as.factor(type)*div,data=all2w)
summary(lm5)
af5<- anova(lm5)
af5ss <- af5$"Sum Sq"
print(cbind(af5,PctExp=af5ss/sum(af5ss)*100))


#4. is divergence or adaptability more causal for parallelism? 
allb<-subset(all,all$pos %in% "SNPs_full" | all$pos %in% "genes_genic")
#rownames(allb) <- paste(allb[,1],substr(allb[,6],1,1),sep="_")
a<- ggplot(allb, aes(x=div, y=prop, color=factor(pos, labels=c("Genes","SNPs")), shape=factor(pos, labels=c("Genes","SNPs")))) + geom_point(size = 3) + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + labs(x="Divergence (neutral Fst)", y = "Proportion of parallels",color="",shape="")
# + scale_x_continuous(breaks = allb$div,labels = allb$contr)
# + geom_text(label=rownames(allb)
b<- ggplot(allb, aes(x=adapt, y=prop, color=factor(pos, labels=c("Genes","SNPs")), shape=factor(pos,labels=c("Genes","SNPs")))) + geom_point(size = 3) + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + labs(x="Evolvability (missense SNPs proportion in candidates vs. all sites)", y = "Proportion of parallels",color="",shape="")
pdf("div_adapt.pdf",height = 7,width = 12,pointsize = 6)
ggarrange(a, b,labels = c("A", "B"), ncol = 2, nrow = 1)
dev.off()
#better not logaritmic
hist(allb$prop)
pdf("div_adapt_plot.pdf",width = 11,height = 11)
par(mfrow=c(2,2))
plot(aov(prop~div*adapt,data=allb))
dev.off()
lm3<-aov(prop~div*adapt+pos,data=allb)
lm3<-glm(prop~div*adapt+pos,data=allb,family = "binomial") ##WHY THIS IS SO INSIGNIFICANT??
summary(lm3)
af <- anova(lm3)
afss <- af$"Sum Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))
#lm3<-aov(log(prop1)~log(div)*log(adapt),data=allb)
#plot(aov(log(prop1)~log(div)*log(adapt),data=allb))



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
FagarasH: '#90EE90'
AlpsH: "#9acd32"
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
v<-venn.diagram(x=list("Alps"=a1,"Fagaras"=b1,"Tatry4x"=d1,"Tatry2x"=c1),paste("bayPass_quartetFst/qf_bf_outlierGenes.tiff",sep=""),fill = c("#E0912F","#065570",'#960C13','#ED7C7F'), lty = "blank",alpha=0.3 , cex=1.4, cat.cex=1.2,main = "Candidate genes",height = 2000,width = 2000)
##find overlaps and annotate
mm<-rbind(a,b,c,d) 
mm1<-as.data.frame(mm[ order(mm[,1]), ])
mm2<-as.data.frame(mm1[duplicated(mm1[,c(1)]),]) 
mm2$one<-1
mm3<-setDT(mm2)[,list(parallels=sum(one)),by=list(Category=`mm1[duplicated(mm1[, c(1)]), ]`)]
##annotate
ann<-fread("../../../Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20181231.txt",h=T,quote="")
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
v<-venn.diagram(x=list("Alps"=a1,"Fagaras"=b1,"Tatry4x"=d1,"Tatry2x"=c1),paste("bayPass_quartetFst/qf_bf_outlierSNPs.tiff",sep=""),fill = c("#E0912F","#065570",'#960C13','#ED7C7F'), lty = "blank",alpha=0.3 , cex=1.4, cat.cex=1.2,main = "Candidate SNPs",height = 2000,width = 2000)
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
s4<-dplyr::select(a,i.ID,i.ann,VEL, ZEP, TKO, TRT, LAC, BAL, INE,SCH,WIL,SUB, BAB,HRA,SPI,DRA,TIS,CAR,ING,KAS,OBI,HCA,GUN,DRG)
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
setwd("/home/aa/alpine/arenosaGenome/selScans/WS1/")
library(topGO)
library(ALL)
library(data.table)
data(ALL)
data(geneList)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.At.tair.db", version = "3.9")
dict<-fread("../../../../Desktop/references/lyrataV2/functions/ALATdict.txt")
a<-fread("outSNPs_Alps0.99.txt")
al<-subset(dict,dict$AL %in% a$ID & !dict$AT %in% "nnn")
al$ATname<-substr(al$AT,1,9)
all <- as.factor(rep("all",nrow(al)))
names(all) <-substr(al$AT,1,9) 
s <- fread("outlierLists/SNPFst_NarrowGenesMin2SNPs.Alps.txt",h=F)
se<-subset(dict,dict$AL %in% s$V1 & !dict$AT %in% "nnn")
sel<-substr(se$AT,1,9)
sel <- as.factor(rep("sel",nrow(se)))
names(sel) <-substr(se$AT,1,9) 

tot<-as.factor(c(rep("1",nrow(se)),rep("0",nrow(al))))
names(tot) <-c(substr(se$AT,1,9),substr(al$AT,1,9))


sampleGOdata <- new("topGOdata",description = "Alps", ontology = "BP", allGenes = tot, nodeSize = 4, annotationFun = annFUN.db, affyLib="org.At.tair.db")

=annFUN.db("org.At.tair.db")
geneSel = sel,
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


### OLD BAYPASS OVERLAP
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


### STEP-BY-STEP TUTORIAL SUPEREXACTTEST
#read file as list
d <- scan("genes_genic.txt", what="", sep="\n")
d <- strsplit(d, "[[:space:]]+")
names(d) <- sapply(d, `[[`, 1)
d <- lapply(d, `[`, -1)
str(d)
(length.gene.sets=sapply(d,length))
total=30000
#expected overlap size
(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))
#probability density distribution of the possible intersection sizes (insert 0:number in the smallest dataset)
(p=sapply(0:32,function(i) dpsets(i, length.gene.sets, n=total)))
(common.genes=intersect(d[[1]], d[[2]], d[[3]],d[[4]],d[[5]]))
#observed intersection
(num.observed.overlap=length(common.genes))
#fold enrichment
(FE=num.observed.overlap/num.expcted.overlap)
#probability density of the observed intersection size
dpsets(num.observed.overlap, length.gene.sets, n=total)
#probability of observing xx or more intersection genes
cpsets(num.observed.overlap-1, length.gene.sets, n=total, lower.tail=FALSE)
### Analyzing all possible intersections
res=supertest(d, n=total)
res$overlap.sizes
#sort.by P value, set, size or degree
#pdf("genes_genic1.pdf")
#plot(res, sort.by="size", margin=c(2,2,2,2), color.scale.pos=c(0.85,1), legend.pos=c(0.9,0.15),keep.empty.intersection#s=F)
#dev.off()
pdf("genes_genic.pdf",width = 9,height = 9)
plot(res, Layout="landscape", degree=2:4, sort.by="size", margin=c(0.5,5,1,2),color.on='forestgreen',keep.empty.intersections=F)
dev.off()
write.csv(summary(res)$Table, file="genes_genic.csv", row.names=FALSE)



########## DMC #############
###1. Generate allele freq table from selected regions
#Works fine but slow locally, better Metacentrum
library(data.table)
library(dplyr)
setwd("/home/aa/alpine/arenosaGenome/selScans/")
ar<-fread(paste('ann/ALLarenosa.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 2)
ha<-fread(paste('ann/ALLhalleri.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
colnames(ar) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",seq(10,ncol(ar),1))
colnames(ha) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",seq(10,ncol(ha),1))

genes <- c("AL1G48930","AL3G42350")
p1<-"BAL"
p2<-"TIS"
p3<-"ZEP"
p4<-"SUB"
dir.create(paste("../../dmc/",p1,p2,p3,p4,sep=""))
 #TODO
#1.extract genes from file
for (id in genes){ #   id = "AL1G48930"
dir.create(paste("../../dmc/",p1,p2,p3,p4,"/",id,sep=""))
setwd(paste("/home/aa/alpine/dmc/",p1,p2,p3,p4,"/",id,sep=""))
###zacatek loopu
ar3<-subset(ar,ar$ID %in% id)
# a<-subset(ar,ar$scaff %in% "scaffold_1" | scaff %in% "scaffold_2")
mid<-min(ar3$start)+round(abs(ar3$start[1]-ar3$start[nrow(ar3)])/2,0)
left<-mid-25000
right<-mid+25000
ar3<-subset(ar,scaff %in% ar3$scaff[1] & start>=left & start<=right)
ha3<-subset(ha,scaff %in% ar3$scaff[1] & start>=left & start<=right) ###halleri
#2.combine together
ar3$end<-ar3$start
ha3$end<-ha3$start
setkey(ha3, scaff, start,end)
a<-foverlaps(x = ar3, y = ha3, type="any")
setkey(ar3, scaff, start,end)
b<-foverlaps(x = ha3, y = ar3, type="any")
bb<-cbind(b[,1],b[,154:195],b[,2:153])
colnames(bb)<-colnames(a)
cc<-rbind(a,bb)
for (i in  1:nrow(cc)){ # i=1
  if (is.na(cc$start[i]))
  {cc[i,4]<-cc[i,46]
  } else {}}
for (i in  1:nrow(cc)){ # i=1
  if (is.na(cc$ID[i]))
  {cc[i,7]<-cc[i,49]
  } else {}}
for (i in  1:nrow(cc)){ # i=1
  if (is.na(cc$ann[i]))
  {cc[i,8]<-cc[i,50]
  } else {}}
for (i in  1:nrow(cc)){ # i=1
  if (is.na(cc$aas[i]))
  {cc[i,9]<-cc[i,51]
  } else {}}
cc<-cc[ order(cc[,1], cc[,4]), ]
cc<-cc[!duplicated(cc[,c('scaff','start')]),]
a<-cc
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
s1<-dplyr::select(a,ID,ann,aas,VEL, ZEP, TKO, TRT, LAC, BAL, INE,SCH,WIL,SUB, BAB,HRA,SPI,DRA,TIS,CAR,ING,KAS,OBI,HCA,GUN,DRG,start)
# sdmc<-dplyr::select(a, scaff,start,BAL, TIS,ZEP,SUB)
for (i in  1:nrow(s1)){ # i=1
  if (is.na(s1[i,4]) & is.na(s1[i,5]) & is.na(s1[i,6]) & is.na(s1[i,7]) & is.na(s1[i,8]) & is.na(s1[i,9]) & is.na(s1[i,10]) & is.na(s1[i,11]) & is.na(s1[i,12]) & is.na(s1[i,13]))
  {s1[i,4:21]<-0
  } else {}}
for (i in  1:nrow(s1)){ # i=1
  if (is.na(s1[i,22]) & is.na(s1[i,23]) & is.na(s1[i,24]) & is.na(s1[i,25]))
  {s1[i,22:25]<-0
  } else {}}
#export
  sdmc<-dplyr::select(s1,start,ID,ann,aas,p1,p2,p3,p4)
  ss<-dplyr::select(s1,p1,p2,p3,p4)
  sdmc$tot<-apply(ss,1,sum)
  sdmc<-subset(sdmc,!sdmc$tot==0)
  sdmc<-subset(sdmc,!sdmc$tot==4)
  write.table(sdmc[,1:8],paste("dmcInput.txt",sep=""),col.names = F,row.names = F,quote = F)
  setwd("../")
  }

### 2. Generate neutral allele frequencies
library(data.table)
library(dplyr)
setwd("/home/aa/alpine/arenosaGenome/neutralomePM/")
#1.extract genes from file
ar3<-fread(paste('ALL.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
ha3<-fread(paste('ALH.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
colnames(ar3) <- c("pop","ploidy","scaff","start","AN","DP",seq(7,ncol(ar3),1))
colnames(ha3) <- c("pop","ploidy","scaff","start","AN","DP",seq(7,ncol(ha3),1))
#2.combine together
ar3$end<-ar3$start
ha3$end<-ha3$start
setkey(ha3, scaff, start,end)
a<-foverlaps(x = ar3, y = ha3, type="any")
a<-subset(a,!is.na(a$pop))
#3. calculate AF for each lineage
i<-which(colnames(a) %in% "7")
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
i<-which(colnames(a) %in% "i.7")
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
s1<-dplyr::select(a,VEL, ZEP, TKO, TRT, LAC, BAL, INE, SCH, WIL, SUB, BAB,HRA,SPI,DRA,TIS,CAR,ING,KAS,OBI,HCA,GUN,DRG,start)
#export
p1<-"LAC"
p2<-"TIS"
p3<-"HCA"
p4<-"DRG"
sdmc<-dplyr::select(s1,start,p1,p2,p3,p4)
ss<-dplyr::select(s1,p1,p2,p3,p4)
sdmc$tot<-apply(ss,1,sum)
sdmc<-subset(sdmc,!sdmc$tot==0)
sdmc<-subset(sdmc,!sdmc$tot==4)
setwd("/home/aa/alpine/dmc/neutralData/")
write.table(sdmc[,2:5],paste(p1,p2,p3,p4,".txt",sep=""),col.names = F,row.names = F,quote = F)



## PREPARE DATA
p1<-"ZEP"
p2<-"SUB"
p3<-"OBI"
p4<-"GUN"
gene<-"AL1G28280"
sampleSizes = c(rep(16, 2),rep(16,2))
rec = 3.7e-08 
Ne = 800000
numPops = 4
selPops = c(1, 3)
numBins = 1000
sels = c(1e-3)
times_sv = c(1e4,1e5,1e6) ## minimum time is split of population pairs
times_sv.source = c(0)  ## maximum time is split of population pairs
times_stagSweeps = c(0,50, 500, 1000,5000) ## maximum time is split of sisters
gs = c(1/(2*Ne))
migs = c(10^-(seq(5, 1, by = -2)), 0.5, 1)
sources = selPops

setwd(paste("/home/aa/alpine/dmc/analysis_2/",p1,p2,p3,p4,"/",gene,sep=""))
a<-read.table("dmcInput.txt",h=F,sep = " ")
aa<-subset(a,a$V2 %in% gene)
left<-min(aa$V1)
right<-max(aa$V1)
se<-as.matrix(t(a[,5:8]))
positions = a$V1
selSite = seq(left-4000, right+4000, length.out = 6)
## GENERATE NEUTRAL MATRIX - F
allFreqs<-as.matrix(t(read.table(paste("../neutralData/",p1,p2,p3,p4,".txt",sep=""))))
neutralF_filename = paste("../neutralData/",p1,p2,p3,p4,"_neutralF",sep="")
source("../calcNeutralF.R")
#plot(F_estimate)
png(paste("../neutralData/",p1,p2,p3,p4,"_neutralF.png",sep=""))
heatmap(F_estimate)
dev.off()
## GENERATE MATRICES WITH selectION - F(s)
# 1 - independent mutations, 2 - migration, 3 - standing variation from source pop, 4 - sv, 5 - stagered sweeps
F_estimate = readRDS(paste("../neutralData/",p1,p2,p3,p4,"_neutralF.RDS",sep=""))
source("../genSelMatrices_individualModes.R")
# model 1
FOmegas_ind = lapply(sels, function(sel) {
  calcFOmegas_indSweeps(sel)
})
saveRDS(FOmegas_ind, "FOmegas_ind.RDS")
# model 2
FOmegas_mig = lapply(sels ,function(sel) {
  lapply(migs, function(mig) {
    lapply(sources, function(my.source) {
      calcFOmegas_mig(sel, mig, my.source)
    })
  })
})
saveRDS(FOmegas_mig, "FOmegas_mig.RDS")
# model 3
FOmegas_sv.source = lapply(sels, function(sel) {
  lapply(gs, function(g) {
    lapply(times_sv.source, function(time) {
      lapply(sources, function(my.source) {
        calcFOmegas_stdVar.source(sel, g, time, my.source)
      })
    })
  })
})
saveRDS(FOmegas_sv.source, "FOmegas_sv_source.RDS")
# model 4
FOmegas_sv = lapply(sels, function(sel) {
  lapply(gs, function(g) {
    lapply(times_sv, function(time) {
        calcFOmegas_stdVar(sel, g, time)
    })
  })
})
saveRDS(FOmegas_sv, "FOmegas_sv.RDS")
# model 5
FOmegas_mig.stagSweeps = lapply(sels, function(sel) {
  lapply(gs, function(g) {
    lapply(times_stagSweeps, function(time_stagSweeps) { ### time
      lapply(sources, function(my.source) {
        calcFOmegas_mig.stagSweeps(sel, g, time_stagSweeps, my.source) ### time
      })
    })
  })
})
saveRDS(FOmegas_mig.stagSweeps, "FOmegas_mig_stagSweeps.RDS")

## GENERATE INVERSES AND DETERMINANTS FOR F(s) MATRICES
library("MASS")
## Neutral model
M = numPops
Tmatrix = matrix(data = rep(-1 / M, (M - 1) * M), nrow = M - 1, ncol = M)
diag(Tmatrix) = (M - 1) / M 
sampleErrorMatrix = diag(1/sampleSizes, nrow = numPops, ncol = numPops)
det_FOmegas_neutral = det(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
saveRDS(det_FOmegas_neutral, "det_FOmegas_neutral.RDS")
inv_FOmegas_neutral = ginv(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
saveRDS(inv_FOmegas_neutral, "inv_FOmegas_neutral.RDS")
## Model 1
FOmegas_ind = readRDS("FOmegas_ind.RDS")
det_FOmegas_ind = lapply(FOmegas_ind, function(sel) {
  lapply(sel, function(dist) {
    det(dist)
  })
})
saveRDS(det_FOmegas_ind, "det_FOmegas_ind.RDS")
inv_FOmegas_ind = lapply(FOmegas_ind, function(sel) {
  lapply(sel, function(dist) {
    ginv(dist)
  })
})
saveRDS(inv_FOmegas_ind, "inv_FOmegas_ind.RDS")
## Model 2
FOmegas_mig = readRDS("FOmegas_mig.RDS")
det_FOmegas_mig = lapply(FOmegas_mig, function(sel) {
  lapply(sel, function(mig) {
    lapply(mig, function(source) {
      lapply(source, function(dist) {
        det(dist)
      })
    })
  })
})
saveRDS(det_FOmegas_mig, "det_FOmegas_mig.RDS")
inv_FOmegas_mig = lapply(FOmegas_mig, function(sel) {
  lapply(sel, function(mig) {
    lapply(mig, function(source) {
      lapply(source, function(dist) {
        ginv(dist)
      })
    })
  })
})
saveRDS(inv_FOmegas_mig, "inv_FOmegas_mig.RDS")
## Model 3
FOmegas_sv.source = readRDS("FOmegas_sv_source.RDS")
det_FOmegas_sv.source = lapply(FOmegas_sv.source, function(sel) {
  lapply(sel, function(g) {
    lapply(g, function(time) {
      lapply(time, function(my.source) {
        lapply(my.source, function(dist) {
          det(dist)
        })
      })
    })
  })
})
saveRDS(det_FOmegas_sv.source, "det_FOmegas_sv_source.RDS")
inv_FOmegas_sv.source = lapply(FOmegas_sv.source, function(sel) {
  lapply(sel, function(g) {
    lapply(g, function(time) {
      lapply(time, function(my.source) {
        lapply(my.source, function(dist) {
          ginv(dist)
        })
      })
    })
  })
})
saveRDS(inv_FOmegas_sv.source, "inv_FOmegas_sv_source.RDS")
## Model 4
FOmegas_sv = readRDS("FOmegas_sv.RDS")
det_FOmegas_sv = lapply(FOmegas_sv, function(sel) {
  lapply(sel, function(g) {
    lapply(g, function(time) {
      lapply(time, function(dist) {
        det(dist)
      })
    })
  })
})
saveRDS(det_FOmegas_sv, "det_FOmegas_sv.RDS")
inv_FOmegas_sv = lapply(FOmegas_sv, function(sel) {
  lapply(sel, function(g) {
    lapply(g, function(time) {
      lapply(time, function(dist) {
          ginv(dist)
      })
    })
  })
})
saveRDS(inv_FOmegas_sv, "inv_FOmegas_sv.RDS")
## Model 5 
FOmegas_mig.stagSweeps = readRDS("FOmegas_mig_stagSweeps.RDS")
det_FOmegas_mig.stagSweeps = lapply(FOmegas_mig.stagSweeps, function(sel) {
  lapply(sel, function(g) {
    lapply(g, function(time) {
      lapply(time, function(my.source) {
        lapply(my.source, function(dist) {
          det(dist)
        })
      })
    })
  })
})
saveRDS(det_FOmegas_mig.stagSweeps, "det_FOmegas_mig_stagSweeps.RDS")
inv_FOmegas_mig.stagSweeps = lapply(FOmegas_mig.stagSweeps, function(sel) {
  lapply(sel, function(g) {
    lapply(g, function(time) {
      lapply(time, function(my.source) {
        lapply(my.source, function(dist) {
          ginv(dist)
        })
      })
    })
  })
})
saveRDS(inv_FOmegas_mig.stagSweeps, "inv_FOmegas_mig_stagSweeps.RDS")

## CALCULATE COMPOSITE LIKELIHOODS
freqs_notRand = se
randFreqs = apply(freqs_notRand, 2, function(my.freqs) {
  if(runif(1) < 0.5) {
    my.freqs = 1 - my.freqs
  }
  my.freqs
})
saveRDS(randFreqs, "selectedRegionAlleleFreqsRand.RDS")
freqs = readRDS("selectedRegionAlleleFreqsRand.RDS")
source("../calcCompositeLike.R")
## Neutral model
det_FOmegas_neutral = readRDS("det_FOmegas_neutral.RDS")
inv_FOmegas_neutral = readRDS("inv_FOmegas_neutral.RDS")
compLikelihood_neutral = lapply(1 : length(selSite), function(j) { # j=1
  calcCompLikelihood_neutral(j, det_FOmegas_neutral, inv_FOmegas_neutral)
})
saveRDS(compLikelihood_neutral, "compLikelihood_neutral.RDS")
## Model 1
det_FOmegas_ind = readRDS("det_FOmegas_ind.RDS")
inv_FOmegas_ind = readRDS("inv_FOmegas_ind.RDS")
compLikelihood_ind = lapply(1 : length(selSite), function(j) {
  lapply(1 : length(sels), function(sel) calcCompLikelihood_1par(j, det_FOmegas_ind,
                                                                 inv_FOmegas_ind, sel))
})
saveRDS(compLikelihood_ind, "compLikelihood_ind.RDS")
## Model 2
det_FOmegas_mig = readRDS("det_FOmegas_mig.RDS")
inv_FOmegas_mig = readRDS("inv_FOmegas_mig.RDS")
compLikelihood_mig = lapply(1 : length(selSite), function(j) {
  lapply(1 : length(sels), function(sel) {
    lapply(1 : length(migs), function(mig) {
      lapply(1 : length(sources), function(my.source) {
        calcCompLikelihood_3par(j, det_FOmegas_mig, inv_FOmegas_mig, sel, mig,
                                my.source)
      })
    })
  })
})
saveRDS(compLikelihood_mig, "compLikelihood_mig.RDS")
## Model 3
det_FOmegas_sv.source = readRDS("det_FOmegas_sv_source.RDS")
inv_FOmegas_sv.source = readRDS("inv_FOmegas_sv_source.RDS")
compLikelihood_sv.source = lapply(1 : length(selSite), function(j) {
  lapply(1 : length(sels), function(sel) {
    lapply(1 : length(gs), function(g) {
      lapply(1 : length(times_sv.source), function(t) {
        lapply(1: length(sources), function(my.source) {
          calcCompLikelihood_4par(j, det_FOmegas_sv.source, inv_FOmegas_sv.source, sel, g, t,
                                  my.source)
        })
      })
    })
  })
})
saveRDS(compLikelihood_sv.source, "compLikelihood_sv_source.RDS")

## Model 4
det_FOmegas_sv = readRDS("det_FOmegas_sv.RDS")
inv_FOmegas_sv = readRDS("inv_FOmegas_sv.RDS")
compLikelihood_sv = lapply(1 : length(selSite), function(j) {
  lapply(1 : length(sels), function(sel) {
    lapply(1 : length(gs), function(g) {
      lapply(1 : length(times_sv), function(t) {
        calcCompLikelihood_3par(j, det_FOmegas_sv, inv_FOmegas_sv, sel, g, t)
      })
    })
  })
})
saveRDS(compLikelihood_sv, "compLikelihood_sv.RDS")

## Model 5
det_FOmegas_mig.stagSweeps = readRDS("det_FOmegas_mig_stagSweeps.RDS")
inv_FOmegas_mig.stagSweeps = readRDS("inv_FOmegas_mig_stagSweeps.RDS")
compLikelihood_mig.stagSweeps = lapply(1 : length(selSite), function(j) {
  lapply(1 : length(sels), function(sel) {
    lapply(1 : length(gs), function(g) {
      lapply(1 : length(times_stagSweeps), function(t) {
        lapply(1: length(sources), function(my.source) {
          calcCompLikelihood_4par(j, det_FOmegas_mig.stagSweeps, inv_FOmegas_mig.stagSweeps, sel, g, t,my.source)
        })
      })
    })
  })
})
saveRDS(compLikelihood_mig.stagSweeps, "compLikelihood_mig_stagSweeps.RDS")

###GET MCL ESTIMATES FOR PARAMETERS
source("../getMCLE.R")
## Model 1
one<-getMCLEind(compLikelihood_ind, selSite, sels)
## Model 2
two<-getMCLEmig(compLikelihood_mig, selSite, sels, migs, sources)
## Model 3
three<-getMCLEsv_source(compLikelihood_sv, selSite, sels, gs, times_sv.source, sources)
## Model 4
four<-getMCLEsv(compLikelihood_sv, selSite, sels, gs, times_sv)
## Model 5
five<-getMCLEmig_stagSweeps(compLikelihood_mig.stagSweeps, selSite, sels, gs, times_stagSweeps, sources)

write.table(as.data.frame(one),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(two),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(three),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(four),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(five),"maxParam.txt",append=T,col.names = T,row.names =F)
## PLOT
#read in composite likelihood files and calculate max for all proposed selected sites
compLikelihood_neutral = readRDS("compLikelihood_neutral.RDS")
compLikelihood_neutral_site = sapply(1 : length(selSite), function(i) {
  max(unlist(compLikelihood_neutral[[i]]))
})
compLikelihood_ind = readRDS("compLikelihood_ind.RDS")
compLikelihood_ind_site = sapply(1 : length(selSite), function(i) {
  max(unlist(compLikelihood_ind[[i]]))
})
compLikelihood_mig = readRDS("compLikelihood_mig.RDS")
compLikelihood_mig_site = sapply(1 : length(selSite), function(i) {
  max(unlist(compLikelihood_mig[[i]]))
})
compLikelihood_sv.source = readRDS("compLikelihood_sv_source.RDS")
compLikelihood_sv.source_site = sapply(1 : length(selSite), function(i) {  
  max(unlist(compLikelihood_sv.source[[i]]))
})
compLikelihood_sv = readRDS("compLikelihood_sv.RDS")
compLikelihood_sv_site = sapply(1 : length(selSite), function(i) {  
  max(unlist(compLikelihood_sv[[i]]))
})
compLikelihood_mig.stagSweeps = readRDS("compLikelihood_mig_stagSweeps.RDS")
compLikelihood_mig.stagSweeps_site = sapply(1 : length(selSite), function(i) {  
  max(unlist(compLikelihood_mig.stagSweeps[[i]]))
})

plot_range = range(c((compLikelihood_ind_site - compLikelihood_neutral_site),
                     (compLikelihood_mig_site - compLikelihood_neutral_site), 
                     (compLikelihood_sv.source_site - compLikelihood_neutral_site),
                     (compLikelihood_sv_site - compLikelihood_neutral_site),
                     (compLikelihood_mig.stagSweeps_site - compLikelihood_neutral_site)))
pdf("compLikel.pdf",width = 9,height = 9,pointsize = 15)
plot(selSite, compLikelihood_ind_site - compLikelihood_neutral_site, type = "b",
     ylim = c(plot_range[1] - 50, plot_range[2] + 50),
     xlab = "Proposed position selected site",
     ylab = "Composite log-likelihood (model - neutral)",main=gene)
lines(selSite, compLikelihood_mig_site - compLikelihood_neutral_site, col = "red",
      type = "b")
lines(selSite, compLikelihood_sv.source_site - compLikelihood_neutral_site, col = "green3",
      type = "b")
lines(selSite, compLikelihood_sv_site - compLikelihood_neutral_site, col = "blue",
      type = "b")
lines(selSite, compLikelihood_mig.stagSweeps_site - compLikelihood_neutral_site, col = "orange2",
      type = "b")
legend("topright", col = c("black", "red","green3", "blue","orange2"),
       lty = 1, sapply(1 : 5, function(i) paste("Model", i)),
       cex = 0.5)
segments(x0 = left,y0 = plot_range[2]+40,x1 = right,y1 = plot_range[2]+40,col = "grey50",lwd = 4)
abline(h = 0,col="grey90",lty=2)
dev.off()

##########
setwd("../")







### PLOT LIKELIHOOD SURFACES
# Time
# Model 3
mcle_sv = getMCLEsv_source(compLikelihood_sv, selSite, sels, gs, times, sources)
compLike_sv_byTime = lapply(1 : length(times), function(time) {
  sapply(1: length(sels), function(sel) {
    sapply(1 : length(gs), function(g) {
      compLikelihood_sv[[which(selSite==mcle_sv[1])]][[sel]][[g]][[time]]
    })
  })
})
profileLike_time_sv = sapply(1: length(compLike_sv_byTime), function(i) {
  max(unlist(compLike_sv_byTime[[i]]))
})
plot(times, profileLike_time_sv, type = "b", xlab = "Time",
     ylab = "Profile composite log-likelihood", main = "Figure 2: Model 3")
abline(v = mcle_sv[4], lty = 2, col = "red")
# selection coefficient
## Model 1
mcle_ind = getMCLEind(compLikelihood_ind, selSite, sels)
profileLike_sel_ind = sapply(1: length(sels), function(i) {
  max(unlist(compLikelihood_ind[[mcle_ind[1]]][[i]]))
})
## Model 2
mcle_mig = getMCLEmig(compLikelihood_mig, selSite, sels, migs, sources)
profileLike_sel_mig = sapply(1: length(sels), function(i) {
  max(unlist(compLikelihood_mig[[mcle_ind[1]]][[i]]))
})
## Model 3
profileLike_sel_sv = sapply(1: length(sels), function(i) {
  max(unlist(compLikelihood_sv[[mcle_sv[1]]][[i]]))
})
pdf("fghf/profileSelCoef.pdf",width = 21,height = 7,pointsize = 18)
par(mfrow = c(1, 3))
plot(sels, profileLike_sel_ind, type = "b", xlab = "Sel",
     ylab = "Profile composite log-likelihood", main = "Model 1")
abline(v = mcle_ind[2], lty = 2, col = "blue")
plot(sels, profileLike_sel_mig, type = "b", xlab = "Sel",
     ylab = "Profile composite log-likelihood", main = "Model 2")
abline(v = mcle_sv[2], lty = 2, col = "blue")
plot(sels, profileLike_sel_sv, type = "b", xlab = "Sel",
     ylab = "Profile composite log-likelihood", main = "Model 3")
abline(v = mcle_sv[2], lty = 2, col = "blue")
dev.off()

### FGSEA
#source("https://bioconductor.org/biocLite.R")
biocLite("fgsea")
library(fgsea)



