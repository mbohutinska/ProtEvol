###WPM
#Stat diversity, pi, thetaw, TajD, 
setwd("/home/aa/alpine/arenosaGenome/4dPM/Until6Jan2020/")
library(data.table)
l<-c('VEL','ZEP','SUB','BAB','SCH','WIL','ING','KAS','LAC','BAL','DRA','TIS','INE','CAR','TKO','TRT','DRA','SPI') #'VEL','ZEP','SUB','BAB','SCH','WIL','ING','KAS','LAC','BAL','DRA','TIS','INE','CAR','TKO','TRT','DRA','SPI',
msum<-matrix(nrow = sum(table(l)), ncol = 16,dimnames =list(c(l),c("ploidy","sampind","ws","num_snps","num_sites","num_singletons","Diversity", "sd_Diversity","ThetaW","sd_ThetaW","Pi","sd_Pi","TajimasD","sd_TajimasD",'lowerCI','upperCI')))
for (lin in l){
  # lin="OBI"
  inp1<-fread(paste(lin,".WS1.0k_MS1_6ind_WPM.txt",sep=""),header=T)
  inp<-head(inp1,-1) 
  inp<-subset(inp,inp$num_sites > 50)
  inp<-subset(inp,inp$num_snps > 15) #HALLERI 2X - 30
  quantile(inp$D,0.975)
  indlin<-which(l %in% paste(lin))
  msum[as.numeric(indlin),1:2]<-as.character(inp[1,2:3])
  msum[as.numeric(indlin),3]<-as.character(inp[1,7])
  msum[as.numeric(indlin),4]<-mean(inp$num_snps)
  msum[as.numeric(indlin),5]<-mean (inp$num_sites)
  msum[as.numeric(indlin),6]<-mean (inp$num_singletons)
  msum[as.numeric(indlin),7]<-mean (inp$Diversity)
  msum[as.numeric(indlin),8]<-sd (inp$Diversity)
  msum[as.numeric(indlin),9]<-mean(inp$ThetaW)
  msum[as.numeric(indlin),10]<-sd (inp$ThetaW)
  msum[as.numeric(indlin),11]<-mean(inp$Pi)
  msum[as.numeric(indlin),12]<-sd (inp$Pi)
  msum[as.numeric(indlin),13]<-mean(inp$D)
  msum[as.numeric(indlin),14]<-sd(inp$D)
  msum[as.numeric(indlin),15]<-quantile(inp$D,0.025)
  msum[as.numeric(indlin),16]<-quantile(inp$D,0.975)
  
}
write.table(msum,append = F,file ="WPM.mean.sd.50_15.txt",quote = F, sep = "\t",col.names = T,row.names = T)

###BPM
setwd("/home/aa/alpine/arenosaGenome/4dPM/Until6Jan2020/")
con<-c('TISSPI','TISHRA','TISGUN','TISDRG','TISCAR','SUBTIS','SUBSPI','SUBKAS','SUBING','SUBHRA','SUBGUN','SUBDRG','SUBDRA','SUBCAR','SUBBAB','SPIGUN','SPIDRG','KASTIS','KASSPI','KASHRA','KASGUN','KASDRG','KASDRA','KASCAR','INGTIS','INGSPI','INGKAS','INGHRA','INGGUN','INGDRG','INGDRA','INGCAR','HRASPI','HRAGUN','HRADRG','GUNDRG','DRATIS','DRASPI','DRAHRA','DRAGUN','DRADRG','DRACAR','CARSPI','CARHRA','CARGUN','CARDRG','BABTIS','BABSPI','BABKAS','BABING','BABHRA','BABGUN','BABDRG','BABDRA','BABCAR')
msum<-matrix(nrow = sum(table(con)), ncol = 12,dimnames =list(c(con),c("ws","num_snps","Rho","sd_Rho","Fst","sd_Fst", "dxy","sd_dxy","AFD","sd_AFD","FixedDiff","sd_FixedDiff")))
for (contr in  con){
  # contr="HCADRG"
  indlin<-which(con %in% paste(contr))
  inp1<-fread(paste(contr,"_WS50000_MS100_BPM.txt",sep=""),header=T)
  inp<-head(inp1,-1) 
  msum[as.numeric(indlin),1 ]<-as.character(inp[1,5])
  msum[as.numeric(indlin),2 ]<-mean(inp$num_snps)
  msum[as.numeric(indlin),3 ]<-mean (inp$Rho)
  msum[as.numeric(indlin),4 ]<-sd (inp$Rho)
  msum[as.numeric(indlin),5 ]<-mean(inp$FstH)
  msum[as.numeric(indlin),6 ]<-sd (inp$FstH)
  msum[as.numeric(indlin),7 ]<-mean(inp$dxy)
  msum[as.numeric(indlin),8 ]<-sd (inp$dxy)
  msum[as.numeric(indlin),9 ]<-mean(inp$AFD)
  msum[as.numeric(indlin),10]<-sd(inp$AFD)
  msum[as.numeric(indlin),11]<-mean(inp$FixedDiff)
  msum[as.numeric(indlin),12]<-sd(inp$FixedDiff)
}
write.table(msum,append = F,file ="BPM.50000.100.mean.sd.V2.txt",quote = F, sep = "\t",col.names = T,row.names = T)


aa<-read.table("BPM.50000.100.mean.sd.V2.txt",h=T)
cor.test(as.numeric(aa$Fst),as.numeric(aa$dxy))
###Get genome values
#cat GUN.WS50.0k_MS10_7ind_WPM.txt | head -n1 >> genome.WS50.0k_MS10_7ind_WPM.txt 
#cat GUN.WS50.0k_MS10_7ind_WPM.txt | grep "Genome" >> genome.WS50.0k_MS10_7ind_WPM.txt
#cat OBI.WS50.0k_MS10_7ind_WPM.txt | grep "Genome" >> genome.WS50.0k_MS10_7ind_WPM.txt
#cat DRG.WS50.0k_MS10_7ind_WPM.txt | grep "Genome" >> genome.WS50.0k_MS10_7ind_WPM.txt
#cat HCA.WS50.0k_MS10_7ind_WPM.txt | grep "Genome" >> genome.WS50.0k_MS10_7ind_WPM.txt
#
#
#cat HCADRG_WS50000_MS10_BPM.txt | head -n1 >> genome_WS50000_MS10_BPM.txt
#cat HCADRG_WS50000_MS10_BPM.txt | grep "Genome" >> genome_WS50000_MS10_BPM.txt
#cat OBIDRG_WS50000_MS10_BPM.txt | grep "Genome" >> genome_WS50000_MS10_BPM.txt
#cat OBIHCA_WS50000_MS10_BPM.txt | grep "Genome" >> genome_WS50000_MS10_BPM.txt
#cat GUNHCA_WS50000_MS10_BPM.txt | grep "Genome" >> genome_WS50000_MS10_BPM.txt
#cat GUNDRG_WS50000_MS10_BPM.txt | grep "Genome" >> genome_WS50000_MS10_BPM.txt
#cat OBIGUN_WS50000_MS10_BPM.txt | grep "Genome" >> genome_WS50000_MS10_BPM.txt

################################################################ 30th of July ARENOSA
setwd("/home/aa/alpine/arenosaGenome/4dPM/")
#full
pop<-c('VEL','ZEP','TRD','SUB','BAB','SCH','WIL','ING','GUL','KAS','KOS','HOC','LAC','BAL','DRA','TIS','INE','CAR','TKO','TRT','HRA','SPI','ZAP','KAM')
con<-c("SCHTRT","TRDSCH","SCHZAP","TRDTRT","TRTZAP","SCHDRA","DRATRT","SCHHOC","SCHKAS","VELSCH","SCHGUL","HOCTRT","KASTRT","SCHLAC","SCHSPI","GULTRT","SCHKAM","SCHBAL","LACTRT","ZEPSCH","SCHTKO","SCHKOS","SCHING","TRDZAP","KOSTRT","SCHTIS","VELTRT","SCHWIL","TRDDRA","BALTRT","INGTRT","SUBSCH","SCHCAR","WILTRT","BABSCH","TRTSPI","TRTKAM","CARTRT","SCHINE","ZEPTRT","INETRT","TISTRT","SCHHRA","DRAZAP","TRDKAS","TRDHOC","TRDGUL","BABTRT","SUBTRT","TRDLAC","TKOTRT","KASZAP","GULZAP","TRDKOS","TRDSPI","TRDBAL","LACZAP","TRDING","TRDWIL","TRTHRA","TRDCAR","TRDTIS","HOCZAP","VELTRD","TRDKAM","HOCDRA","KASDRA","VELDRA","GULDRA","TRDINE","ZEPTRD","KOSZAP","DRASPI","VELZAP","LACDRA","BALZAP","INGZAP","SPIZAP","DRAKAM","ZAPKAM","KOSDRA","CARZAP","TRDTKO","ZEPDRA","WILZAP","VELKAS","INGDRA","ZEPZAP","TRDBAB","VELGUL","DRATKO","TISZAP","HOCLAC","VELHOC","KASLAC","WILDRA","KASSPI","TRDSUB","INEZAP","GULLAC","VELKOS","TRDHRA","SUBDRA","VELLAC","GULSPI","BABDRA","HOCSPI","KOSLAC","DRAINE","DRAHRA","GULKAS","BALDRA","LACSPI","BABZAP","DRACAR","SUBZAP","VELSPI","KOSSPI","KASHOC","KASKOS","VELING","ZEPGUL","ZEPKAS","TKOZAP","KASTKO","GULKAM","INGLAC","ZEPKOS","VELKAM","SUBKAS","VELBAL","GULKOS","ZEPHOC","GULTKO","BABKAS","LACTKO","HOCTKO","VELCAR","HOCKAM","DRATIS","SUBGUL","INGGUL","ZEPLAC","GULHOC","GULTIS","HRAZAP","VELWIL","VELTIS","LACKAM","INGKAS","INGSPI","BABGUL","KOSTKO","KASKAM","GULINE","KASTIS","KOSKAM","VELINE","GULBAL","INGHOC","BALSPI","SUBKOS","KOSBAL","KOSTIS","KOSHOC","WILKAS","ZEPSPI","CARSPI","WILLAC","VELSUB","KASBAL","SUBLAC","WILSPI","TISSPI","KASCAR","KOSCAR","WILGUL","BABKOS","BABLAC","LACCAR","GULCAR","GULHRA","BABHOC","KASINE","SUBHOC","KASHRA","VELBAB","INESPI","HOCTIS","SPIKAM","HOCINE","SUBSPI","KOSINE","ZEPING","HOCBAL","HOCHRA","WILHOC","LACTIS","INGKOS","LACHRA","HOCCAR","VELTKO","WILKOS","LACINE","INGTKO","TKOSPI","KOSHRA","LACBAL","INGKAM","BABSPI","INGTIS","VELHRA","CARTKO","BABING","SUBING","INGBAL","BALTKO","INGCAR","HRASPI","WILTKO","ZEPCAR","TISTKO","INETKO","ZEPWIL","TKOKAM","BABBAL","ZEPBAL","ZEPTIS","ZEPKAM","ZEPTKO","INGINE","BABCAR","INGHRA","CARKAM","ZEPINE","BABWIL","SUBBAL","BALKAM","WILKAM","BABTIS","SUBCAR","ZEPBAB","WILING","BABKAM","SUBTKO","BABTKO","WILBAL","TISKAM","BABINE","WILCAR","ZEPSUB","BALCAR","SUBKAM","INEKAM","WILTIS","SUBWIL","VELZEP","CARHRA","SUBINE","SUBTIS","BALHRA","INEHRA","TISHRA","TISCAR","WILHRA","WILINE","HRAKAM","ZEPHRA","BALTIS","BALINE","SUBHRA","TISINE","BABHRA","TKOHRA","INECAR","SUBBAB")
#subsampled
pop<-c('VEL','ZEP','SUB','BAB','SCH','WIL','ING','KAS','LAC','BAL','DRA','TIS','INE','CAR','TKO','TRT','HRA','SPI')
con<-c('ZEPWIL','ZEPTRT','ZEPTKO','ZEPTIS','ZEPSUB','ZEPSPI','ZEPSCH','ZEPLAC','ZEPKAS','ZEPING','ZEPINE','ZEPHRA','ZEPDRA','ZEPCAR','ZEPBAL','ZEPBAB','WILTRT','WILTKO','WILTIS','WILSPI','WILLAC','WILKAS','WILING','WILINE','WILHRA','WILDRA','WILCAR','WILBAL','VELZEP','VELWIL','VELTRT','VELTKO','VELTIS','VELSUB','VELSPI','VELSCH','VELLAC','VELKAS','VELING','VELINE','VELHRA','VELDRA','VELCAR','VELBAL','VELBAB','TRTSPI','TRTHRA','TKOTRT','TKOSPI','TKOHRA','TISTRT','TISTKO','TISSPI','TISINE','TISHRA','TISCAR','SUBWIL','SUBTRT','SUBTKO','SUBTIS','SUBSPI','SUBSCH','SUBLAC','SUBKAS','SUBING','SUBINE','SUBHRA','SUBDRA','SUBCAR','SUBBAL','SUBBAB','SCHWIL','SCHTRT','SCHTKO','SCHTIS','SCHSPI','SCHLAC','SCHKAS','SCHING','SCHINE','SCHHRA','SCHDRA','SCHCAR','SCHBAL','LACTRT','LACTKO','LACTIS','LACSPI','LACINE','LACHRA','LACDRA','LACCAR','LACBAL','KASTRT','KASTKO','KASTIS','KASSPI','KASLAC','KASINE','KASHRA','KASDRA','KASCAR','KASBAL','INGTRT','INGTKO','INGTIS','INGSPI','INGLAC','INGKAS','INGINE','INGHRA','INGDRA','INGCAR','INGBAL','INETRT','INETKO','INESPI','INEHRA','INECAR','HRASPI','DRATRT','DRATKO','DRATIS','DRASPI','DRAINE','DRAHRA','DRACAR','CARTRT','CARTKO','CARSPI','CARHRA','BALTRT','BALTKO','BALTIS','BALSPI','BALINE','BALHRA','BALDRA','BALCAR','BABWIL','BABTRT','BABTKO','BABTIS','BABSPI','BABSCH','BABLAC','BABKAS','BABING','BABINE','BABHRA','BABDRA','BABCAR','BABBAL')
msum<-matrix(nrow = sum(table(pop)), ncol = sum(table(pop)),dimnames =list(c(pop),c(pop)))
stats<-c("FstH","Rho") #"FstWC","dxy",
for(s in stats){
  for(contr in  con){ # contr="HRAKAM"
    p1<-substr(contr,1,3)
    p2<-substr(contr,4,7)
    p1i<-which(pop %in% paste(p1))
    p2i<-which(pop %in% paste(p2))
    d<-read.table(paste(contr,"_WS50000_MS1_BPM.txt",sep=""),h=T)
    ds<-which(colnames(d) %in% paste(s))
    msum[p1i,p2i]<-d[nrow(d),ds]
    msum[p2i,p1i]<-d[nrow(d),ds]
  }
  write.table(msum,paste("genome_",s,"_WS50000_MS1_BPM.txt",sep=""),quote = F)}

#neutralome
setwd("/home/aa/alpine/arenosaGenome/4dPM/")
pop<-c(l<-c('SUB','BAB','ING','KAS','DRA','TIS','HRA','SPI','CAR','GUN','DRG'))

con<-c('TISSPI','TISHRA','TISGUN','TISDRG','TISCAR','SUBTIS','SUBSPI','SUBKAS','SUBING','SUBHRA','SUBGUN','SUBDRG','SUBDRA','SUBCAR','SUBBAB','SPIGUN','SPIDRG','KASTIS','KASSPI','KASHRA','KASGUN','KASDRG','KASDRA','KASCAR','INGTIS','INGSPI','INGKAS','INGHRA','INGGUN','INGDRG','INGDRA','INGCAR','HRASPI','HRAGUN','HRADRG','GUNDRG','DRATIS','DRASPI','DRAHRA','DRAGUN','DRADRG','DRACAR','CARSPI','CARHRA','CARGUN','CARDRG','BABTIS','BABSPI','BABKAS','BABING','BABHRA','BABGUN','BABDRG','BABDRA','BABCAR')
msum<-matrix(nrow = sum(table(pop)), ncol = sum(table(pop)),dimnames =list(c(pop),c(pop)))
stats<-c("FstH") #"FstWC","dxy",,"Rho"
for(s in stats){
  for(contr in  con){ # contr="HRAKAM"
    p1<-substr(contr,1,3)
    p2<-substr(contr,4,7)
    p1i<-which(pop %in% paste(p1))
    p2i<-which(pop %in% paste(p2))
    d<-read.table(paste(contr,"_WS50000_MS100_BPM.txt",sep=""),h=T)
    ds<-which(colnames(d) %in% paste(s))
    msum[p1i,p2i]<-d[nrow(d),ds]
    msum[p2i,p1i]<-d[nrow(d),ds]
  }
  write.table(msum,paste("genome_",s,"_WS50000_MS100_BPM.txt",sep=""),quote = F)}

###Correlate Fsts
a<-read.table("genome_FstH_WS50000_MS1_BPM_inclMean.txt",h=T)
a[1,1]<-NA
a[4,4]<-NA
a[7,7]<-NA
a[10,10]<-NA
a[13,13]<-NA
a[16,16]<-NA
a[21,21]<-NA
a[24,24]<-NA
#foothill
a1<-as.matrix(a[c(4,10,16,20,24,27,29),c(4,10,16,20,24,27,29)])
#alpine
a2<-as.matrix(a[c(1,7,13,19,21,28,30),c(1,7,13,19,21,28,30)])
cor(c(a1), c(a2),use = "complete.obs")
plot(c(a1)~ c(a2))
pdf("divergenceCorr.pdf")
plot(c(a1)~ c(a2),ylab = "Foothill", xlab = "Alpine")
dev.off()
text(c(a1)~ c(a2), labels = paste(row.names(a1)), pos = 4)
write.table(a1,"divergenceFoothill.txt",quote = F)



###AFS
setwd("/home/aa/alpine/arenosaGenome/AFS4d0d/")
l<-c('VEL','ZEP','TRD','SUB','BAB','SCH','WIL','ING','GUL','KAS','KOS','HOC','LAC','BAL','DRA','TIS','INE','CAR','TKO','TRT','HRA','SPI','ZAP','KAM')
for (lin in l){ # lin="VEL"
a1<-read.table(paste(lin,"_0d_AFS.txt",sep="")) #0d
b1<-read.table(paste(lin,"_4d_AFS.txt",sep=""))
as<-sum(a1[2:(max(a1$V3)+1),4])
bs<-sum(b1[2:(max(b1$V3)+1),4])
s1<-(as+bs)/2
a1$V5<-a1$V4/as
b1$V5<-b1$V4/bs
zero_d<-a1$V5[1:max(a1$V3)+1]
four_d<-b1$V5[1:max(a1$V3)+1]
test<-rbind(zero_d,four_d)
png(paste(lin,"_SFSunfoldedSubsampl1_0d4d.png"),width = 1300,height = 850,pointsize = 20)
barplot(test,beside=T,names.arg = seq(1,(max(a1$V3)),1), axes = T, main = paste(lin),col = c("orange","blue"),ylim = c(0,0.35),las=2)
legend("topright", c("0dg","4dg"),fill = c("orange","blue"),bg = NA) #0dg #x = 0, y = (b1[(max(b1$V3)+1),5])
dev.off()}
##To have it in nice table:
a2<-as.data.frame(t(a1)[c(3,5),0:max(a1$V3)+1], row.names = c("AC","0d"))
b2<-as.data.frame(t(b1)[c(3,5),0:max(b1$V3)+1], row.names = c("AC","4d"))
c1<-rbind(a2[1:2,],b2[2,])
colnames(c1) <- as.character(unlist(c1[1,]))
c1 = c1[-1, ]

###Venn diagram
# install.packages('VennDiagram')
library(VennDiagram)
setwd("/home/aa/alpine/halleriGenome/selectionScans/")

#0.999
a<-read.table("HCADRG_WS380_MS10_BPM_0.999tile_OutAnnot.csv",sep = "=", header = F,skip = 1)
b<-read.table("OBIGUN_WS380_MS10_BPM_0.999tile_OutAnnot.csv",sep = "=", header = F,skip = 1)
v<-venn.diagram(x=list("Alps"=unlist(b$V3),"Fagaras"=unlist(a$V3)),"venn0_999ws380.tiff", lty = "blank", fill = c("blue1","red2"),alpha=0.3 , main = "Overlap of selection candidates",cex=1, cat.cex=1,sub = "Arabidopsis halleri",height = 2000,width = 2500)

intersect(a$V3, b$V3)

###Plot stat around windows
i3<- subset (i2, i2$POS %in% 9810000:10550000)




  ################
  con<-c("WCATET", "SECWCA", "SECTET", "PANWCA", "PANTET", "PANSEC", "PANDIN", "PANBAL", "DINWCA", "DINTET", "DINSEC", "DINBAL", "CROWCA", "CROTET", "CROSEC", "CROPAN", "CRODIN", "CROBAL", "BALWCA", "BALTET", "BALSEC")
  msum<-matrix(nrow = sum(table(con)), ncol = 16,dimnames =list(c(con),c("Smin","S1Q","Smedian","Smean","S3Q","Smax","Ssd", "Nmin","N1Q","Nmedian","Nmean","N3Q","Nmax","Nsd","Squantile95","Squantile99")))
  #al <- fread("../ALcodes.txt",header = F)
  #mat <- matrix(nrow = nrow(al), ncol = sum(table(con)),dimnames = list(c(al$V1),c(con))) 
  for (contr in  con){
    # contr="DINTET"
    indlin<-which(con %in% paste(contr))
    inp<-fread(paste(contr,"_WS1_MS1_BPM.txt",sep=""),header=T)
    s<-subset(inp,inp$AAS %in% "synonymous_variant")
    n<-subset(inp,inp$AAS %in% "missense_variant")
    sq <- quantile(abs(s$AFD),probs=c(.99))
    msum[as.numeric(indlin),1:6]<-summary (abs(s$AFD))
    msum[as.numeric(indlin),7]<-sd(abs(s$AFD))
    msum[as.numeric(indlin),8:13]<-summary (abs(n$AFD))
    msum[as.numeric(indlin),14]<-sd(abs(n$AFD))
    msum[as.numeric(indlin),15] <-quantile(abs(s$AFD),probs=c(.95))
    msum[as.numeric(indlin),16]<-sq
    
  }
  write.table(msum,append = F,file ="AFSStatsPerContrast.txt",quote = F, sep = "\t",col.names = T,row.names = T)
  
  
############## ALPINE 2
  ###AFS
  setwd("/home/aa/2alpine/PopStructure/AFS/")
  l<-c('SUB','ZEP','LOM','NKM','SPN','SUN','DRG','GUN','HCA','OBI','HFA','HFF','HIA','HIF')
  pdf(paste("SFSunfoldedSubsampl1_4d.pdf"),width = 6,height = 6,pointsize = 14)
  
  for (lin in l){ # lin="DRG"
    # a1<-read.table(paste(lin,"_0d_AFS.txt",sep="")) #0d
    b1<-read.table(paste(lin,"_4d_AFS.txt",sep=""))
    #   as<-sum(a1[2:(max(a1$V3)+1),4])
    bs<-sum(b1[2:(max(b1$V3)+1),4])
    s1<-bs #(as+bs)/2
    #a1$V5<-a1$V4/as
    b1$V5<-b1$V4/bs
    #zero_d<-a1$V5[1:max(a1$V3)+1]
    four_d<-b1$V5[1:max(b1$V3)+1]
    test<-rbind(four_d)
    barplot(test,names.arg = seq(1,(max(b1$V3)),1), axes = T, main = paste(lin),col = c("blue"),ylim = c(0,0.9),las=2)
    # legend("topright", c("0dg","4dg"),fill = c("orange","blue"),bg = NA) #0dg #x = 0, y = (b1[(max(b1$V3)+1),5])
  } 
  dev.off()

  
  ###BPM
  setwd("/home/aa/2alpine/PopStructure/BPM/")
  
  pop<-c('SUB','ZEP','LOM','NKM','SPN','SUN','DRG','GUN','HCA','OBI','HFA','HFF','HIA','HIF')
  con<-c('DRGGUN','DRGHCA','DRGHFA','DRGHFF','DRGHIA','DRGHIF','DRGLOM','DRGNKM','DRGOBI','DRGSPN','DRGSUB','DRGSUN','DRGZEP','GUNHCA','GUNHFA','GUNHFF','GUNHIA','GUNHIF','GUNLOM','GUNNKM','GUNOBI','GUNSPN','GUNSUB','GUNSUN','GUNZEP','HCAHFA','HCAHFF','HCAHIA','HCAHIF','HCALOM','HCANKM','HCAOBI','HCASPN','HCASUB','HCASUN','HCAZEP','HFAHFF','HFAHIA','HFAHIF','HFALOM','HFANKM','HFAOBI','HFASPN','HFASUB','HFASUN','HFAZEP','HFFHIA','HFFHIF','HFFLOM','HFFNKM','HFFOBI','HFFSPN','HFFSUB','HFFSUN','HFFZEP','HIAHIF','HIALOM','HIANKM','HIAOBI','HIASPN','HIASUB','HIASUN','HIAZEP','HIFLOM','HIFNKM','HIFOBI','HIFSPN','HIFSUB','HIFSUN','HIFZEP','LOMNKM','LOMOBI','LOMSPN','LOMSUB','LOMSUN','LOMZEP','NKMOBI','NKMSPN','NKMSUB','NKMSUN','NKMZEP','OBISPN','OBISUB','OBISUN','OBIZEP','SPNSUB','SPNSUN','SPNZEP','SUBSUN','SUBZEP','SUNZEP')
  msum<-matrix(nrow = sum(table(pop)), ncol = sum(table(pop)),dimnames =list(c(pop),c(pop)))
  stats<-c("FstH","dxy") #"FstWC","dxy",
  for(s in stats){
    for(contr in  con){ # contr="HRAKAM"
      p1<-substr(contr,1,3)
      p2<-substr(contr,4,7)
      p1i<-which(pop %in% paste(p1))
      p2i<-which(pop %in% paste(p2))
      d<-read.table(paste(contr,"_WS50000_MS1_BPM.txt",sep=""),h=T)
      ds<-which(colnames(d) %in% paste(s))
      msum[p1i,p2i]<-d[nrow(d),ds]
      msum[p2i,p1i]<-d[nrow(d),ds]
    }
    write.table(msum,paste("genome_",s,"_WS50000_MS1_BPM.txt",sep=""),quote = F)}
  
  
tajD<-c(0.08, 0.09, -0.03, -0.02, -0.15, 0.36, 0.48, 0.47, 0.66, 0.41, 0.01, 0.08, 0.02, 0.11, 0.03, -0.16, 0.01, 0.29, 0.54, 0.50, 0.17, 0.14)
quantile(tajD,0.025)  
quantile(tajD,0.975)  
mean(tajD)
