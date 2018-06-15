########N/S analysis per gene
#  #works on Metacentrum
#  library(data.table)
#  setwd("/home/aa/JICAutumn2016/finalAnalysis29Apr")
#  lineages<-c('CRO','PAN','DIN','BAL','SEC','WCA','DIP','TET')
#  #lineages<-c('TET')
#  for (lin in  lineages){
#    system(command = paste("grep -E 'missense_variant|synonymous_variant' data/",lin,".table.repol.txt > data/",lin,".synNon  .txt",sep=""))
#    inp<-fread(paste("data/",lin,".synNon.txt",sep=""),header=F)
#    inp[inp=="-9"]<-NA
#    for (g in readLines("lists/ALcodesAll.txt"))
#    {# g='AL8G38150'
#      mat <- matrix(nrow = 1, ncol = 6)
#      gene<-subset(x = inp,subset = inp$V7 %in% g)
#      if (nrow(subset(gene,gene$V8 %in% "missense_variant"))==0) {n<-0} else 
#        {n<-sum(subset(gene,gene$V8 %in% "missense_variant")[,10:ncol(gene)],na.rm = T)}
#      if (nrow(subset(gene,gene$V8 %in% "synonymous_variant"))==0) {s<-0} else 
#        {s<-sum(subset(gene,gene$V8 %in% "synonymous_variant")[,10:ncol(gene)],na.rm = T)}
#      if ((n+s)<5) {check<-NA} else 
#        {check<-1}
#      ns<-as.numeric(n)/as.numeric(s)
#      mat[1,1]<-as.character(gene[1,1])
#      mat[1,2]<-g
#      mat[1,3]<-n
#      mat[1,4]<-s
#      mat[1,5]<-ns
#      mat[1,6]<-check
#      write.table(mat,append = T,file = "results/NSperGene.AllLin.Allgenes.txt",quote = F, sep = "\t",col.names = F,row.names =   F)}
#    mat1 <- matrix(nrow = 1, ncol = 6)
#    n<-sum(subset(inp,inp$V8 %in% "missense_variant")[,10:ncol(inp)],na.rm = T)
#    s<-sum(subset(inp,inp$V8 %in% "synonymous_variant")[,10:ncol(inp)],na.rm = T)
#    ns<-as.numeric(n)/as.numeric(s)
#    mat1[1,1]<-as.character(inp[1,1])
#    mat1[1,2]<-"genome"
#    mat1[1,3]<-n
#    mat1[1,4]<-s
#    mat1[1,5]<-ns
#    mat1[1,6]<-"1"
#    write.table(mat1,append = T,file = "results/NSgenomeWide.AllLin.txt",quote = F, sep = "\t",col.names = F,row.names = F)}
#  ###########the end of Metacentrum

###per meiotic genes
lineages<-c('CRO','PAN','DIN','BAL','SEC','WCA','TET')
if (file.exists("results/meiotic.PropTest.AllLin.txt")) file.remove("results/meiotic.PropTest.AllLin.txt")
if (file.exists("results/meiotic.ModeOfSelection.NRP.AllLin.txt")) file.remove("results/meiotic.ModeOfSelection.NRP.AllLin.txt")
ns<-fread("results/NSperGene.AllLin.Allgenes.txt",header=F)
gwns<-read.table("results/NSgenomeWide.AllLin.txt",header=F)
genes<-read.table("lists/ALmeiotic.ALcode.tab.Name.txt",header=F)
mtot<-matrix(nrow = nrow(genes), ncol = sum(table(lineages)),dimnames =list(c(paste(genes$V1,"_",genes$V2,sep="")),c(lineages)))
for (lin in  lineages){
indlin<-which(lineages %in% paste(lin))
nslin<-subset(ns,ns$V1 %in% paste(lin))
gwnslin<-subset(gwns,gwns$V1 %in% paste(lin))
png(filename=paste("histGW.meiotic.N_S.",lin,".png",sep=""), 
    units ="px", 
    width=1300*2, 
    height=875*2, 
    pointsize=16*2) 
freq<-hist(log(nslin$V5),breaks = 300,main = paste("dN/dS of all genes -", lin, sep=" "),xlab = "log(dN/dS)",freq = F,col = "grey")
abline(v=log(gwnslin$V5),col="red",lwd=6)
for (g in genes$V1){
  value<-log(as.numeric(subset(nslin,nslin$V2 %in% paste(g))[1,5]))
  arrows(y0 = max(freq$density),x0 = value,x1 = value,y1 = (max(freq$density)-0.06*(max(freq$density))),length = 0.2,col = "blue",lwd = 3)
  arrows(5,40,5,50, length = 0.15, angle = 30, lty=1,code=1 , lwd=2, col="purple")}
dev.off()
#vyndat meiozni - cely radek + prop.test #spocitat prop.test
mat<-matrix(nrow = nrow(genes),ncol=9)
for (g in genes$V1){
  index<-which(genes$V1 %in% paste(g))
  line<-subset(nslin,nslin$V2 %in% paste(g))
  tot<-line$V3+line$V4
  if (nrow(line)<1) 
   {mat[as.numeric(index),1]<-genes$V2[index]
    mat[as.numeric(index),2]<-lin
    mat[as.numeric(index),3]<-g
    mat[as.numeric(index),4:9]<-"NA"} 
  else if (line$V6 %in% NA)
    {mat[as.numeric(index),1]<-genes$V2[index]
    mat[as.numeric(index),2:7]<-as.matrix(line)
    mat[as.numeric(index),8]<-tot
    mat[as.numeric(index),9]<-"NA"} 
  else
  {p<-prop.test(x = line$V3,n = tot,p = gwnslin$V3/(gwnslin$V3+gwnslin$V4))
  mat[as.numeric(index),1]<-genes$V2[index]
  mat[as.numeric(index),2:7]<-as.matrix(line)
  mat[as.numeric(index),8]<-tot
  mat[as.numeric(index),9]<-p$p.value}
  #make a matrix prot/lineages - only Neg/Pos
  if (mat[as.numeric(index),7] %in% NA)
  {mtot[as.numeric(index),as.numeric(indlin)]<-"NA NA"} else if (mat[as.numeric(index),9] >= 0.05)
  {mtot[as.numeric(index),as.numeric(indlin)]<-paste("R ",formatC(p$p.value, format = "e", digits = 1), sep="")} else if (mat[as.numeric(index),9] < 0.05 & mat[as.numeric(index),6] < gwnslin$V3/(gwnslin$V3+gwnslin$V4))
  {mtot[as.numeric(index),as.numeric(indlin)]<-paste("N ",formatC(p$p.value, format = "e", digits = 1), sep="")} else
  {mtot[as.numeric(index),as.numeric(indlin)]<-paste("P ",formatC(p$p.value, format = "e", digits = 1), sep="")}
}
write.table(mat,append = T,file = "results/meiotic.PropTest.AllLin.txt",quote = F, sep = "\t",col.names = F,row.names = F)}
#write.table(t(as.data.frame(c("gene",rep(lineages,each=2)))),append = F,file = "meiotic.ModeOfSelection.NRP.AllLin.txt",quote = F, sep = "\t",col.names = F,row.names=F)
write.table(mtot,append = T,file = "results/meiotic.ModeOfSelection.NRP.AllLin.txt",quote = F, sep = "\t",col.names = F,row.names =T)

###########MKT
#To minimize the impact of slightly deleterious mutations, it has been proposed to exclude polymorphisms that are below a certain cutoff frequency, such as <8% or <5%
#We have suggested that polymorphisms below 15% be excluded from MK-type analyses when there is evidence that some nonsilent polymorphisms are slightly deleterious.
#removing low-frequency polymorphisms increases the variance of the estimate of a, but this is a price that has to be paid to obtain a less biased estimate. Hopefully, by providing a recommended cutoff frequency, we will remove the temptation to search for the frequency that yields the highest value of a because this is statistically difficult to defend.
library(data.table)
setwd("~/JICAutumn2016/finalAnalysis29Apr/")
con<-c("TETDIP","PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA")
# con<-c("CROPAN")
al <- fread("lists/ALcodes.txt",header = F)
#to filter for missing data: 0.9 = 90% nomissing
mffg<-0.49
#to filter for minimum derived allele frequency (0.15 - recommended)
mdafP<-0.1
##pridat az budu vedet GW treshold!!
mdafD<-0.6
for (contr in  con){
  system(command = paste("head -n 1 data/",contr,"_WS1_MS1_BPM.txt > data/",contr,".synNon_WS1_MS1_BPM.txt",sep=""))
  system(command = paste("grep -E 'missense_variant|synonymous_variant' data/",contr,"_WS1_MS1_BPM.txt >> data/",contr,".synNon_WS1_MS1_BPM.txt",sep=""))
  #  contr<-"DINTET"
  mat <- matrix(nrow = nrow(al), ncol = 10,dimnames = list(c(),c("contr","gene","pN","pS","dN","dS","alfa","NI","DoS","name")))
  p<-fread(paste("head -n -1 data/",contr,".synNon_WS1_MS1_BPM.txt",sep = ""),header = T)
  # Histogram Colored (blue and red)
#  png(filename = paste("results/AFS_", contr,"_",mdafP,"_",mdafD,".png"),width = 1300,height = 850,pointsize = 24)
#  contrN<-subset(x = p, p$ann %in% "missense_variant")
#  contrS<-subset(x = p, p$ann %in% "synonymous_variant")
#  hist(contrN$AC0/contrN$AN0, col=rgb(1,0.7,0,0.6),xlim=c(0,1), xlab="AF",breaks = 100,main = paste("AFS in ",contr," pair"#,sep=""))
#  hist(contrS$AC0/contrS$AN0, col=rgb(0,0,1,0.6), add=T,breaks = 100)
#  legend(legend = c("synonymous", "nonsynonymous"),fill = c("blue", "orange"),x = "topright")
#  abline(v = mdafP,col="red")
#  dev.off()
  #filter missing data and minDAF
  p1<-subset(x = p,subset = p$AN0 > max(p$AN0)*mffg & p$AN1 > max(p$AN1)*mffg & abs(p$AFD) >= mdafP & abs(p$AFD) < mdafD)
  #filter missing data and almost fixed variants 
  d1<-subset(x = p,subset = p$AN0 > max(p$AN0)*mffg & p$AN1 > max(p$AN1)*mffg & abs(p$AFD) >= mdafD)
  
  for (gene in readLines("lists/ALmeiotic.ALcode.tab.Name.txt"))
  {# g='AL4G46460 SHOC'
    g<-substr(gene,1,9)
    name<-substr(gene,11,30)
    index<-which(al$V1 %in% paste(g))
    pN<-nrow(subset(x = p1,subset = p1$ALcode %in% g & p1$ann %in% "missense_variant"))
    pS<-nrow(subset(x = p1,subset = p1$ALcode %in% g & p1$ann %in% "synonymous_variant"))
    dN<-nrow(subset(x = d1,subset = d1$ALcode %in% g & d1$ann %in% "missense_variant"))
    dS<-nrow(subset(x = d1,subset = d1$ALcode %in% g & d1$ann %in% "synonymous_variant"))
    if (pS>0 & dN>0) 
    {alfa=1-(dS*pN)/(dN*pS)
    NI=(pN/pS)/(dN/dS)} else 
    {alfa<-"NA"
    NI<-"NA"}
    if (pS+pN>0 & dS+dN>0) 
    {dos=dN/(dN+dS)-pN/(pN+pS)} else 
    {dos<-"NA"}
    mat[as.numeric(index),1]<-contr
    mat[as.numeric(index),2]<-g
    mat[as.numeric(index),3]<-pN
    mat[as.numeric(index),4]<-pS
    mat[as.numeric(index),5]<-dN
    mat[as.numeric(index),6]<-dS
    mat[as.numeric(index),7]<-alfa #to estimate the fraction (α) of substitutions driven to fixation by positive selection at the functional sites
    mat[as.numeric(index),8]<-NI
    mat[as.numeric(index),9]<-dos #Stoletzki, 2011
    mat[as.numeric(index),10]<-name}
  write.table(mat,append = F,file = paste("results/",contr,"_mafdP_",mdafP,"_mafdD_",mdafD,"_MKT.txt",sep=""),quote = F, sep = "\t",col.names = T,row.names = F)}

###########HighFst AASs###########
#  ###  -- THIS IS INCORPORATED INTO SCANTOOLS -- ###
#  setwd("~/JICAutumn2016/finalAnalysis29Apr/")
#  library(data.table)
#  con<-c("PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA","CROPAN","CRODIN","CROSEC" # ,"CROBAL","CROWCA")
#  con<-c("TETDIP")
#  msum<-matrix(nrow = sum(table(con)), ncol = 16,dimnames =list(c(con),c("Smin","S1Q","Smedian","Smean","S3Q","Smax","Ssd",   "Nmin","N1Q","Nmedian","Nmean","N3Q","Nmax","Nsd","Squantile95","Squantile99")))
#  al <- fread("lists/ALcodesAll.txt",header = F)
#  mat <- matrix(nrow = nrow(al), ncol = sum(table(con)),dimnames = list(c(al$V1),c(con))) 
#  for (contr in  con){
#    indlin<-which(con %in% paste(contr))
#    inp<-fread(paste("head -n -1 data/",contr,"_WS1_MS1_BPM.txt",sep=""),header=T)
#    s<-subset(inp,inp$AAS %in% "synonymous_variant")
#    n<-subset(inp,inp$AAS %in% "missense_variant")
#    sq <- quantile(s$FstH,probs=c(.99))
#    msum[as.numeric(indlin),1:6]<-summary (s$FstH)
#    msum[as.numeric(indlin),7]<-sd(s$FstH)
#    msum[as.numeric(indlin),8:13]<-summary (n$FstH)
#    msum[as.numeric(indlin),14]<-sd(n$FstH)
#    msum[as.numeric(indlin),15] <-quantile(s$FstH,probs=c(.95))
#    msum[as.numeric(indlin),16]<-sq
#    png(filename=paste("results/",contr,"_FstHperNandSHist.png", sep=""), 
#        units = "px",
#        width=1800, 
#        height=1000, 
#        pointsize=16)
#    hist(log(s$FstH), freq = F, breaks = 200, xlab =paste("log(Fst) - ",contr, sep=""), col=rgb(1,0,0,0.5), main =NA,xlim=c(-7  ,0))
#    hist(log(n$FstH), freq=F, breaks = 200, col=rgb(0,0,1,0.5), add=T,xlim=c(-7,0))
#    abline(v=log(quantile(s$FstH,probs=c(.99))),col="red")
#    legend("topleft", c("synonymous","nonsynonymous"), fill=c( "red", "blue"))
#    box()
#    dev.off()
#    for (g in readLines("lists/ALcodesAll.txt"))
#    {# g='AL8G38150'
#      index<-which(al$V1 %in% paste(g))
#      gene<-subset(x = inp,subset = inp$ALcode %in% g & inp$AAS %in% "missense_variant")
#      
#      if (nrow(gene)==0) {high<-0} else 
#      {high<-nrow(subset(gene,gene$FstH >= as.numeric(sq)))}
#      mat[index,indlin]<-high
#    }}
#  write.table(mat,append = F,file = "results/FstHighperGene.AllLin.Allgenes.txt",quote = F, sep = "\t",col.names = T,row.names   = T)
#  write.table(msum,append = F,file ="results/FstStatsPerContrast.txt",quote = F, sep = "\t",col.names = T,row.names = T)
#  ###  -- THE END OF SCANTOOLS PART  -- ###

######TETRAPLOIDS#######
library(data.table)
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/")
###COPY METACENTRUM OUTPUT TO /DATA!!1
system("bash gatherFstperGene.sh")
quant = 0.99 #quant = 0.999
hfall<-read.table(paste("results/FstHighperGene.Allgenes.",quant,".txt", sep=""),header=T,row.names = 1)
colnames(hfall)[1]<-"TETDIP"
fst<-read.table(paste("results/FstStatsPerContrast.",quant,".txt", sep=""),header=T)
colnames(fst)[15]<-"Squantile"
genes<-read.table("lists/ALmeiotic.ALcode.tab.Name.txt",header=F)
genes<-read.table("lists/ALcodesAll.txt",h=F)
con<-c("TETDIP")
if (file.exists(paste("results/meiotic.highFstsitesList.",quant,".txt", sep=""))) file.remove(paste("results/meiotic.highFstsitesList.",quant,".txt", sep=""))
for (contr in  con){
  indCon<-which(con %in% paste(contr))
  png(filename=paste("results/histGW.meiotic.highFst.",contr,quant,".png",sep=""), 
      units ="px", 
      width=1300*2, 
      height=875*2, 
      pointsize=16*2) 
  freq<-hist(subset(hfall[,indCon],hfall[,indCon] > 0),breaks = 100,main = paste("number of high diff AASs in all genes -", contr, sep=" "),xlab = "no. high diff AASs",freq = F,col = "grey") #xlim = c(1,max(log(hfall[,indCon]))
  for (g in genes$V1){#  g="AL6G30890"
    index<-which(rownames(hfall) %in% paste(g))
    value<-hfall[index,indCon]
    if (value>0) 
    {arrows(y0 = max(freq$density),x0 = value+0.02,x1 = value+0.02,y1 = (max(freq$density)-0.12*(max(freq$density))),length =  0.2,col = "blue",lwd = 3)} else 
    {}}
  dev.off()
  mat<-matrix(data = as.matrix(subset(hfall, rownames(hfall) %in% genes$V1)[names(hfall)==con]),nrow = nrow(genes),ncol=length(con),dimnames = list(paste(genes$V1,"_",genes$V2,sep=""),c(con)))
  inp<-fread(paste("head -n -1 data/",contr,".synNon_WS1_MS1_BPM.txt",sep=""),header=T)
  q<-fst$Squantile[indCon]
  n<-subset(inp,inp$ann %in% "missense_variant" & inp$FstH >= q)
  ge<-subset(n,n$ALcode %in% genes$V1)
  ge1<-cbind(ge$scaff,ge$end)
  write.table(cbind(contr,ge1),append = T,file = paste("results/meiotic.highFstsitesList.",quant,".txt", sep=""),quote = F, sep = "\t",col.names = F,row.names = F)
  write.table(cbind(ge$scaff,":",ge$end),append = T,file = paste("results/quantile",quant,".intervals",sep=""),quote = F, sep = "",col.names = F,row.names = F)
  }
d<-read.table(paste("results/meiotic.highFstsitesList.",quant,".txt",sep=""),header = F)
d1<-d[ order(d[,2], d[,3]), ][,c(2,3)] 
d2<-d1[!duplicated(d1[,c('V2','V3')]),]
mat<-matrix(data=-9,nrow = nrow(d2),ncol=(9+sum(table(con))))
for (contr in  con){
  indCon<-which(con %in% paste(contr))
  for (i in 1:nrow(d2)){# i=5
    line<-d2[i,]
  if (nrow(subset(inp,inp$scaff %in% line[,1] & inp$end %in% line[,2]))==0)
    {l2<-c(NA)} else
    {l<-subset(inp,inp$scaff %in% line[,1] & inp$end %in% line[,2])
    l1<-l[,2:10]
    l2<-l$FstH
    mat[i,1:9]<-as.matrix(l1[1,])
    }
    mat[i,indCon+9]<-l2  }}
write.table(mat,append = F,file = paste("results/meiotic.highFstFinalSummary.",quant,".txt",sep=""),quote = F, sep = "\t",col.names = F,row.names = F)  
#####ANNOTATE
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/")
library(data.table)
library(stringr)
library(readr)
a<-read.table("lists/pI_hyd.txt",header = T)
n<-read.table("lists/ALmeiotic.ALcode.tab.Name.txt",h=F)
i<-read.table("pai_Malvidae/Alignment_Identities_allSites.txt",h=F)
d<-read.table(paste("results/meiotic.highFstFinalSummary.",quant,".txt",sep=""),header = F,col.names = c("scaf","V2","V3","V4","V5","alcode","aas","ann","pos","TETDIP"))
d$V2 <- d$V3 <- d$V4 <- d$V5 <- NULL
d$pos<-d$pos+1
d$anc<-substr(x = d$aas,start = 1,stop = 3)
d$der<-str_sub(string = d$aas, start= -3)
indAnc<-which(colnames(d) %in% "anc")
indDer<-which(colnames(d) %in% "der")
d$anc_pI<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,4])
d$der_pI<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,4])
d$delta_pI<-abs(d$anc_pI-d$der_pI)
#hist(d$delta_pI,breaks = 50)
d$anc_hydCA<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,2])
d$der_hydCA<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,2])
d$delta_hydCA<-abs(d$anc_hydCA-d$der_hydCA)
#hist(d$delta_hydCA,breaks = 50)
d$anc_hydCC<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,3])
d$der_hydCC<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,3])
d$delta_hydCC<-abs(d$anc_hydCC-d$der_hydCC)
#hist(d$delta_hydCC,breaks = 50)
d$prot<-apply(X = d,1, function(x) subset(x = n, subset = n$V1 %in% x[2])[1,2])
d$aas_pos<-as.character(parse_number(as.character(d$aas)))
indAl<-which(colnames(d) %in% "alcode")
indAP<-which(colnames(d) %in% "aas_pos")
d$al_ident<-apply(X = d,1, function(x) subset(x = i, subset = i$V6 %in% x[indAl] & as.character(i$V1) %in% as.character(x[indAP]))[1,3])
d$al_hydr<-apply(X = d,1, function(x) subset(x = i, subset = i$V6 %in% x[indAl] & as.character(i$V1) %in% as.character(x[indAP]))[1,4])
d$al_pI<-apply(X = d,1, function(x) subset(x = i, subset = i$V6 %in% x[indAl] & as.character(i$V1) %in% as.character(x[indAP]))[1,5])
indProt<-which(colnames(d) %in% "prot")
indDpI<-which(colnames(d) %in% "delta_pI")
indDhyC<-which(colnames(d) %in% "delta_hydCC")
indDhyA<-which(colnames(d) %in% "delta_hydCA")
indPAI<-which(colnames(d) %in% "al_ident")
indPAH<-which(colnames(d) %in% "al_hydr")
indPApI<-which(colnames(d) %in% "al_pI")
d1<-d[ order(d[,2], d[,5]), ]
d2<-d1[,c(indAl,5,indProt,3,6,indDpI,indDhyA,indDhyC, indPAI,indPAH,indPApI)]
#d2a<-subset(d2, d2$delta_pI!=0)
d3<-as.data.frame(d1[,1])
d3$t<-":"
d3$pos<-d1[,5]
write.table(d1, paste("results/meiotic.highFstFinalSummaryAnn",quant,".txt",sep=""),append = F,quote = F,sep = "\t",row.names = F)
write.table(d2, paste("results/meiotic.highFstFinalSummaryAnn.Subset",quant,".txt",sep=""),append = F,quote = F,sep = "\t",row.names = F)
write.table(d3, paste("results/quantile",quant,".intervals",sep=""),append = F,quote = F,sep = "",row.names = F,col.names = F)



######DIPLOIDS########
library(data.table)
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/")
###COPY METACENTRUM OUTPUT TO /DATA!!1
system("bash gatherFstperGene.sh")
quant = 0.999
hfall<-read.table(paste("results/FstHighperGene.Allgenes.",quant,".txt", sep=""),header=T,row.names = 1)
colnames(hfall)[1]<-"TETDIP"
fst<-read.table(paste("results/FstStatsPerContrast.",quant,".txt", sep=""),header=T)
colnames(fst)[15]<-"Squantile"
genes<-read.table("lists/ALmeiotic.ALcode.tab.Name.txt",header=F)
con<-c("PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA")
if (file.exists(paste("results/meiotic.highFstsitesList.",quant,".txt", sep=""))) file.remove(paste("results/meiotic.highFstsitesList.",quant,".txt", sep=""))
for (contr in  con){ # contr="PANDIN"
  indCon<-which(con %in% paste(contr))
  png(filename=paste("results/histGW.meiotic.highFst.",contr,quant,".png",sep=""), 
      units ="px", 
      width=1300*2, 
      height=875*2, 
      pointsize=16*2) 
  freq<-hist(subset(hfall[,indCon],hfall[,indCon] > 0),breaks = 100,main = paste("number of high diff AASs in all genes -", contr, sep=" "),xlab = "no. high diff AASs",freq = F,col = "grey") #xlim = c(1,max(log(hfall[,indCon]))
  for (g in genes$V1){#  g="AL6G30890"
    index<-which(rownames(hfall) %in% paste(g))
    value<-hfall[index,indCon]
    if (value>0) 
    {arrows(y0 = max(freq$density),x0 = value+0.02,x1 = value+0.02,y1 = (max(freq$density)-0.12*(max(freq$density))),length =  0.2,col = "blue",lwd = 3)} else 
    {}}
  dev.off()
  #vyndat meiozni
  mat<-matrix(data = as.matrix(subset(hfall, rownames(hfall) %in% genes$V1)[names(hfall)==con]),nrow = nrow(genes),ncol=length(con),dimnames = list(paste(genes$V1,"_",genes$V2,sep=""),c(con)))
  inp<-fread(paste("head -n -1 data/",contr,".synNon_WS1_MS1_BPM.txt",sep=""),header=T)
  q<-fst$Squantile[indCon]
  n<-subset(inp,inp$ann %in% "missense_variant" & inp$FstH >= q)
  ge<-subset(n,n$ALcode %in% genes$V1)
  ge1<-cbind(ge$scaff,ge$end)
  write.table(cbind(contr,ge1),append = T,file = paste("results/meiotic.highFstsitesList.",quant,".txt", sep=""),quote = F, sep = "\t",col.names = F,row.names = F)
}
#get the list - DIPLOIDS - TODO:scaffold instead of ALCcode!!!!
con<-c("PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA","CROPAN","CRODIN","CROSEC","CROBAL","CROWCA")
d<-read.table(paste("results/meiotic.highFstsitesList.",quant,".txt",sep=""),header = F)
d1<-d[ order(d[,2], d[,3]), ][,c(2,3)] 
d2<-d1[!duplicated(d1[,c('V2','V3')]),]
mat<-matrix(data=-9,nrow = nrow(d2),ncol=(9+sum(table(con))))
for (contr in  con){
  indCon<-which(con %in% paste(contr))
  inp<-fread(paste("head -n -1 data/",contr,"_WS1_MS1_BPM.txt",sep=""),header=T)
  for (i in 1:nrow(d2)){
    # i=5
    line<-d2[i,]
    if (nrow(subset(inp,inp$scaff %in% line[,1] & inp$end %in% line[,2]))==0)
    {#l1<-(matrix(rep(NA,9),rep(1,9)))
      l2<-c(NA)} else
      {l<-subset(inp,inp$scaff %in% line[,1] & inp$end %in% line[,2])
      l1<-l[,2:10]
      l2<-l$FstH
      mat[i,1:9]<-as.matrix(l1[1,])
      }
    mat[i,indCon+9]<-l2  }}
write.table(mat,append = F,file = paste("results/meiotic.highFstFinalSummary.",quant,".txt",sep=""),quote = F, sep = "\t",col.names = F,row.names = F)  
#####ANNOTATE
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/")
library(data.table)
library(stringr)
library(readr)
a<-read.table("lists/pI_hyd.txt",header = T)
n<-read.table("lists/ALmeiotic.ALcode.tab.Name.txt",h=F)
i<-read.table("pai_Malvidae/Alignment_Identities_allSites.txt",h=F)
d<-read.table(paste("results/meiotic.highFstFinalSummary.",quant,".txt",sep=""),header = F,col.names = c("scaf","V2","V3","V4","V5","alcode","aas","ann","pos","PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA","CROPAN","CRODIN","CROSEC","CROBAL","CROWCA"))
#for high effect:
d$V2 <- d$V3 <- d$V4 <- d$V5 <- NULL
d$pos<-d$pos+1
d$anc<-substr(x = d$aas,start = 1,stop = 3)
d$der<-str_sub(string = d$aas, start= -3)
indAnc<-which(colnames(d) %in% "anc")
indDer<-which(colnames(d) %in% "der")
d$anc_pI<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,4])
d$der_pI<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,4])
d$delta_pI<-abs(d$anc_pI-d$der_pI)
#hist(d$delta_pI,breaks = 50)
d$anc_hydCA<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,2])
d$der_hydCA<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,2])
d$delta_hydCA<-abs(d$anc_hydCA-d$der_hydCA)
#hist(d$delta_hydCA,breaks = 50)
d$anc_hydCC<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,3])
d$der_hydCC<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,3])
d$delta_hydCC<-abs(d$anc_hydCC-d$der_hydCC)
#hist(d$delta_hydCC,breaks = 50)
d$prot<-apply(X = d,1, function(x) subset(x = n, subset = n$V1 %in% x[2])[1,2])
d$aas_pos<-as.character(parse_number(as.character(d$aas)))
indAl<-which(colnames(d) %in% "alcode")
indAP<-which(colnames(d) %in% "aas_pos")
d$al_ident<-apply(X = d,1, function(x) subset(x = i, subset = i$V6 %in% x[indAl] & as.character(i$V1) %in% as.character(x[indAP]))[1,3])
d$al_hydr<-apply(X = d,1, function(x) subset(x = i, subset = i$V6 %in% x[indAl] & as.character(i$V1) %in% as.character(x[indAP]))[1,4])
d$al_pI<-apply(X = d,1, function(x) subset(x = i, subset = i$V6 %in% x[indAl] & as.character(i$V1) %in% as.character(x[indAP]))[1,5])
indProt<-which(colnames(d) %in% "prot")
indDpI<-which(colnames(d) %in% "delta_pI")
indDhyC<-which(colnames(d) %in% "delta_hydCC")
indDhyA<-which(colnames(d) %in% "delta_hydCA")
indPAI<-which(colnames(d) %in% "al_ident")
indPAH<-which(colnames(d) %in% "al_hydr")
indPApI<-which(colnames(d) %in% "al_pI")
d1<-d[ order(d[,2], d[,5]), ]
d2<-d1[,c(indAl,5,indProt,3,6:(5+length(con)),indDpI,indDhyA,indDhyC, indPAI,indPAH,indPApI)]
d3<-as.data.frame(d1[,1])
d3$t<-":"
d3$pos<-d1[,5]
write.table(d1, paste("results/meiotic.highFstFinalSummaryAnn",quant,".txt",sep=""),append = F,quote = F,sep = "\t",row.names = F)
write.table(d2, paste("results/meiotic.highFstFinalSummaryAnn.Subset",quant,".txt",sep=""),append = F,quote = F,sep = "\t",row.names = F)
write.table(d3, paste("results/quantile",quant,".intervals",sep=""),append = F,quote = F,sep = "",row.names = F,col.names = F)
###Get outlier - DIPLOIDS
quant = 0.999
out<-read.table(paste("results/meiotic.highFstFinalSummaryAnn.Subset",quant,".txt",sep=""),header = T)
fst<-read.table(paste("results/FstStatsPerContrast.",quant,".txt", sep=""),header=T)
con<-c("PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA","CROPAN","CRODIN","CROSEC","CROBAL","CROWCA")
mat<-matrix(data=-9,nrow = nrow(out),ncol=(length(con)),dimnames =list(c(), c(con)))
for (contr in  con){# contr="PANDIN"
 indCon<-which(con %in% paste(contr))
 q<-fst$Squantile[indCon]
 mat[,indCon]<-out[,(4+indCon)]/q
}
o<-cbind(out,mat)
write.table(o, paste("results/meiotic.highFstFinalSummaryAnn.Subset.relativeFst",quant,".txt",sep=""),append = F,quote = F,sep = "\t",row.names = F)


#############High effect changes
library(data.table)
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/")
###COPY METACENTRUM OUTPUT TO /DATA!!1
system("bash gatherFstperGene.sh")
quant = 0.95
#hfall<-read.table(paste("results/FstHighperGene.Allgenes.",quant,".txt", sep=""),header=T,row.names = 1)
#colnames(hfall)[1]<-"TETDIP"
fst<-read.table(paste("results/FstStatsPerContrast.",quant,".txt", sep=""),header=T)
colnames(fst)[15]<-"Squantile"
genes<-read.table("lists/ALmeiotic.ALcode.tab.Name.txt",header=F)
con<-c("TETDIP","PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA")
if (file.exists(paste("results/meiotic.highFstsitesList.",quant,".txt", sep=""))) file.remove(paste("results/meiotic.highFstsitesList.",quant,".txt", sep=""))
for (contr in  con){ # contr="PANDIN"
  indCon<-which(con %in% paste(contr))
  #vyndat meiozni
  #mat<-matrix(data = as.matrix(subset(hfall, rownames(hfall) %in% genes$V1)[names(hfall)==con]),nrow = nrow(genes),ncol=length(con),dimnames = list(paste(genes$V1,"_",genes$V2,sep=""),c(con)))
  inp<-fread(paste("head -n -1 data/",contr,"_WS1_MS1_BPM.txt",sep=""),header=T)
  q<-fst$Squantile[indCon]
  k<-subset(inp,inp$ann %in% "splice_acceptor_variant" | inp$ann %in% "splice_donor_variant" | inp$ann %in% "splice_region_variant" | inp$ann %in% "stop_lost" | inp$ann %in% "start_lost" | inp$ann %in% "stop_gained")
  n<-subset(k, k$FstH >= q)
  ge<-subset(n,n$ALcode %in% genes$V1)
  ge1<-cbind(ge$scaff,ge$end)
  write.table(cbind(contr,ge1),append = T,file = paste("results/meiotic.highFstsitesList.",quant,".txt", sep=""),quote = F, sep = "\t",col.names = F,row.names = F)
  write.table(cbind(ge$scaff,":",ge$end),append = T,file = paste("results/quantile",quant,".intervals",sep=""),quote = F, sep = "",col.names = F,row.names = F)
}
con<-c("TETDIP","PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA")
d<-read.table(paste("results/meiotic.highFstsitesList.",quant,".txt",sep=""),header = F)
d1<-d[ order(d[,2], d[,3]), ][,c(2,3)] 
d2<-d1[!duplicated(d1[,c('V2','V3')]),]
mat<-matrix(data=-9,nrow = nrow(d2),ncol=(9+sum(table(con))))
for (contr in  con){
  indCon<-which(con %in% paste(contr))
  inp<-fread(paste("head -n -1 data/",contr,"_WS1_MS1_BPM.txt",sep=""),header=T)
  for (i in 1:nrow(d2)){
    # i=5
    line<-d2[i,]
    if (nrow(subset(inp,inp$scaff %in% line[,1] & inp$end %in% line[,2]))==0)
    {#l1<-(matrix(rep(NA,9),rep(1,9)))
      l2<-c(NA)} else
      {l<-subset(inp,inp$scaff %in% line[,1] & inp$end %in% line[,2])
      l1<-l[,2:10]
      l2<-l$FstH
      mat[i,1:9]<-as.matrix(l1[1,])
      }
    mat[i,indCon+9]<-l2  }}
write.table(mat,append = F,file = paste("results/meiotic.highFstFinalSummary.",quant,".txt",sep=""),quote = F, sep = "\t",col.names = F,row.names = F)  
#####ANNOTATE
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/")
library(data.table)
library(stringr)
library(readr)
a<-read.table("lists/pI_hyd.txt",header = T)
n<-read.table("lists/ALmeiotic.ALcode.tab.Name.txt",h=F)
i<-read.table("pai_Malvidae/Alignment_Identities_allSites.txt",h=F)
#for high effect:
d<-read.table(paste("results/meiotic.highFstFinalSummary.",quant,".txt",sep=""),header = F,col.names = c("scaf","V2","V3","V4","V5","alcode","aas","ann","pos","TETDIP","PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA"))
d$V2 <- d$V3 <- d$V4 <- d$V5 <- NULL
d$pos<-d$pos+1
d$anc<-substr(x = d$aas,start = 1,stop = 3)
d$der<-str_sub(string = d$aas, start= -3)
indAnc<-which(colnames(d) %in% "anc")
indDer<-which(colnames(d) %in% "der")
d$anc_pI<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,4])
d$der_pI<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,4])
d$delta_pI<-abs(d$anc_pI-d$der_pI)
#hist(d$delta_pI,breaks = 50)
d$anc_hydCA<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,2])
d$der_hydCA<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,2])
d$delta_hydCA<-abs(d$anc_hydCA-d$der_hydCA)
#hist(d$delta_hydCA,breaks = 50)
d$anc_hydCC<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,3])
d$der_hydCC<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,3])
d$delta_hydCC<-abs(d$anc_hydCC-d$der_hydCC)
#hist(d$delta_hydCC,breaks = 50)
d$prot<-apply(X = d,1, function(x) subset(x = n, subset = n$V1 %in% x[2])[1,2])
d$aas_pos<-as.character(parse_number(as.character(d$aas)))
indAl<-which(colnames(d) %in% "alcode")
indAP<-which(colnames(d) %in% "aas_pos")
d$al_ident<-apply(X = d,1, function(x) subset(x = i, subset = i$V6 %in% x[indAl] & as.character(i$V1) %in% as.character(x[indAP]))[1,3])
d$al_hydr<-apply(X = d,1, function(x) subset(x = i, subset = i$V6 %in% x[indAl] & as.character(i$V1) %in% as.character(x[indAP]))[1,4])
d$al_pI<-apply(X = d,1, function(x) subset(x = i, subset = i$V6 %in% x[indAl] & as.character(i$V1) %in% as.character(x[indAP]))[1,5])
indProt<-which(colnames(d) %in% "prot")
indDpI<-which(colnames(d) %in% "delta_pI")
indDhyC<-which(colnames(d) %in% "delta_hydCC")
indDhyA<-which(colnames(d) %in% "delta_hydCA")
indPAI<-which(colnames(d) %in% "al_ident")
indPAH<-which(colnames(d) %in% "al_hydr")
indPApI<-which(colnames(d) %in% "al_pI")
d1<-d[ order(d[,2], d[,5]), ]
d2<-d1[,c(indAl,5,indProt,3,4,6:(5+length(con)),indDpI,indDhyA,indDhyC, indPAI,indPAH,indPApI)] #for High effect
d3<-as.data.frame(d1[,1])
d3$t<-":"
d3$pos<-d1[,5]
write.table(d1, paste("results/meiotic.highFstFinalSummaryAnn",quant,".txt",sep=""),append = F,quote = F,sep = "\t",row.names = F)
write.table(d2, paste("results/meiotic.highFstFinalSummaryAnn.Subset",quant,".txt",sep=""),append = F,quote = F,sep = "\t",row.names = F)


###################To get position - AAS dictionary
k<-read.table("lists/Pos_AAS_300AllSynNon.tab.txt", h=T)
k$ALcode<-substr(k$ANN....GENE,1,9)
k$pp<-substr(k$ANN....HGVS_P,3,12)
k$pp1<-gsub(",","",x = k$pp,fixed = T)
k$aas<-gsub(".","",x = k$pp1,fixed = T)
write.table(file = "lists/Pos_AAS_300AllSynNon.tab.txt", x = k[,c(1,2,5,8)], sep = "\t", quote = F,append = F,row.names = F)

###################To get list of identities from Geneious output 
###newer
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/pai_Malvidae/")
library(data.table)
library(stringr)
d<-read.table("listCSV.txt",header = F)
for (g in d$V1){
  dd<-read.table(paste("csvMalvidae/",g,sep=""),h=T,sep = ",")
  gg<-substr(g,1,9)
  dd$V3<-gg
  write.table(dd[1:nrow(dd),],"Alignment_Identities_allSites.txt",append = T,quote = F,row.names = F,col.names = F,sep = "\t")}

###################Haplotype allele frequencies
#get AF per lineage, combine into one matrix, also number anc, der
library (dplyr)
library(stringr)
library(readr)
setwd(dir = "/home/aa/JICAutumn2016/finalAnalysis29Apr/haplotypes/")
pops<-c('HNE', 'KZL', 'SZI', 'BEL', 'BIH', 'FOJ', 'GOR', 'RZA', 'VID', 'MIE', 'PRE', 'HNI', 'SNO', 'TRD', 'VEL', 'DRA', 'LAC', 'TZI', 'BRD', 'WEK', 'CHO', 'RFT', 'SWA', 'BGS', 'GUL', 'HOC', 'KAS', 'KOS', 'SCH', 'SPI', 'TKO', 'TRE', 'TRT', 'ZAP', 'KOW', 'STE', 'TBG')
contr<-c("diploid","tetraploid","other")
colors<-c("red3", "blue2", "grey")
#1.neutral
c1<-c('HNE', 'KZL', 'SZI', 'BEL', 'BIH', 'FOJ', 'GOR', 'RZA', 'VID', 'MIE', 'PRE', 'HNI', 'SNO', 'TRD', 'VEL')
#2.selected
c2<-c('DRA', 'LAC', 'TZI', 'BRD', 'WEK', 'CHO', 'RFT', 'SWA', 'BGS', 'GUL', 'HOC', 'KAS', 'KOS', 'SCH', 'SPI', 'TKO', 'TRE', 'TRT', 'ZAP', 'KOW', 'STE', 'TBG')

nnn<-read.table(paste(pops[1],".table.repol.txt",sep=""),h=F)
mtot<-matrix(nrow = nrow(nnn), ncol = sum(table(pops))+4,dimnames =list(c(),c("scaff","pos","alcode","AAS",pops)))
mtot[,1]<-as.character(nnn[,3])
mtot[,2]<-as.character(nnn[,4])
mtot[,3]<-as.character(nnn[,7])
mtot[,4]<-as.character(nnn[,9])
for (pop in  pops){ # pop="KZL"
  ip<-which(pops %in% paste(pop))
  a<-read.table(paste(pop,".table.repol.txt",sep=""), header =F)
  a[ a == "-9" ] <- NA
  a$na<-apply(is.na(a), 1, sum)
  a$an<-((ncol(a)-10)-a$na)*as.numeric(a[1,2])
  a$af<-rowSums(a[,10:(ncol(a)-2)],na.rm = T) / a$an
  mtot[,ip+4]<-a$af}
write.table(mtot,"AFsperLin.txt",append = F,quote = F,row.names = F,col.names = T)
##repolarize
k<-read.table("AFsperLin.txt",h=T)
k$anc<-substr(x = k$AAS,start = 1,stop = 3)
k$der<-str_sub(string = k$AAS, start= -3)
k$aas_pos<-as.character(parse_number(as.character(k$AAS)))
k1<-select(.data = k,c1)
k2<-select(.data = k,c2)
k$c1<-rowSums(k1,na.rm = T)/ncol(k1)
k$c2<-rowSums(k2,na.rm = T)/ncol(k2)
#e<-0
for (i in  1:nrow(k)){ # i=1
  if (k$c1[i]>k$c2[i])
  {k[i,5:(ncol(k)-10)] <- 1 - k[i,5:(ncol(k)-10)]
  a<-k$anc[i]
  k$anc[i]<-k$der[i]
  k$der[i]<-a
  #e<-e+1
  #print(k[i,1:2])
  } else {}}
#print(e)
write.table(k[1:(ncol(k)-2)],"AFsperLinRepol.txt",append = F,quote = F,row.names = F,col.names = T)
#calculate haplotype allele frequencies
k<-read.table("AFsperLinRepol.txt",h=T)
g<-read.table("PopGeogr.txt",h=T)
for (al in unique(k$alcode)){ # al<-"AL1G10680"
 a1<-subset(k,k$alcode %in% al)
 if (nrow(a1)==1)
 {} else {
 a2<-as.data.frame(t(a1[,5:(ncol(a1)-3)]))
 #neutral
 g$c1<-1-(apply(a2, 1, max))
 #selected
 g$c2<-apply(a2, 1, min)
 #recombined
 g$oth<-abs(1-g$c1-g$c2)
 n<-read.table("ALmeiotic.ALcode.tab.Name.InclZyp.txt")
 name<-as.matrix(subset(n,n$V1 %in% al))[1,2]
 write.table(g,paste(al,"_",name,".hapl.txt",sep=""),quote = F,sep = "\t",row.names = F)}}
#MAP
library(RgoogleMaps)
library(argosfilter)
library(mapplots)
list<-list.files(path = ".", pattern = "*.hapl.txt")
for (l in list){
  df<-read.table(paste(l), header=T)
  lat <- c(43,55.8) #define our map's ylim
  lon <- c(4.5,27) #define our map's xlim
  name<-strsplit(l, split='[.]')[[1]][1]
  prot<-strsplit(name, split='[_]')[[1]][2]
  terrmap<-GetMap.bbox(lon,lat,maptype= "terrain", destfile = "terrain.png",path = paste0("&style=feature:road|element:all|visibility:off","&style=feature:administrative.country|element:geometry.stroke|visibility:off","&style=feature:all|element:labels|visibility:off"), MINIMUMSIZE = T)
  size<-terrmap$size
  pdf(file=paste("figures/",name,".hapl.pdf", sep=""),width = size[1]/72, height = size[2]/72)
  PlotOnStaticMap(terrmap)
  df<-na.omit(df)
  for (i in 1:nrow(df)){
    d1<-distance(df$lat[i], df$lat[i], lon[1], lon[2]) 
    d2<-distance(df$lat[i], df$lat[i], lon[1], df$lon[i])
    a1<-size[1]/d1
    a2<-d2*a1
    a3<- -(size[1]/2)+a2+df$offsetx[i]
    b1<-size[2]/(lat[2]-lat[1])
    b2<-(df$lat[i]-lat[1])*b1
    b3<- -(size[2]/2)+b2+df$offsety[i]
    terrmap<-add.pie(z=c(df$c1[i],df$c2[i],df$oth[i]), x=a3, y=b3, radius=10, col = colors, labels=NA)
    terrmap<-text( x=a3, y=b3, labels = df$pop[i],pos = df$pos[i],col = as.character(df$col[i]), offset = 0.83,cex = 0.7)}
  terrmap<-legend("topleft", inset=.005, title=paste(prot," haplotypes",sep=""),contr,fill = colors, horiz=F, cex=1,bg = 'white')
  dev.off()}


###TODO
##################################################
#for 3 haplotypes
##################################################
library (dplyr)
library(stringr)
library(readr)
setwd(dir = "/home/aa/JICAutumn2016/finalAnalysis29Apr/haplotypes/")
pops<-c('HNE', 'KZL', 'SZI', 'BEL', 'BIH', 'FOJ', 'GOR', 'RZA', 'VID', 'MIE', 'PRE', 'HNI', 'SNO', 'TRD', 'VEL', 'DRA', 'LAC', 'TZI', 'BRD', 'WEK', 'CHO', 'RFT', 'SWA', 'BGS', 'GUL', 'HOC', 'KAS', 'KOS', 'SCH', 'SPI', 'TKO', 'TRE', 'TRT', 'ZAP', 'KOW', 'STE', 'TBG')
contr<-c("pannonian","diploid","tetraploid","other")
colors<-c("orange","red3", "blue2", "grey")
#1.neutral
c1<-c('BEL', 'BIH', 'FOJ', 'GOR', 'RZA', 'VID', 'MIE', 'PRE', 'HNI', 'SNO', 'TRD', 'VEL')
#2.main selected
c2<-c('DRA', 'LAC', 'TZI', 'BRD', 'WEK', 'CHO', 'RFT', 'SWA', 'BGS', 'GUL', 'HOC', 'KAS', 'KOS', 'SCH', 'SPI', 'TKO', 'TRE', 'TRT', 'ZAP', 'KOW', 'STE', 'TBG')
#3. also selected
c3<-c('HNE', 'KZL', 'SZI')


nnn<-read.table(paste(pops[1],".table.repol.txt",sep=""),h=F)
mtot<-matrix(nrow = nrow(nnn), ncol = sum(table(pops))+4,dimnames =list(c(),c("scaff","pos","alcode","AAS",pops)))
mtot[,1]<-as.character(nnn[,3])
mtot[,2]<-as.character(nnn[,4])
mtot[,3]<-as.character(nnn[,7])
mtot[,4]<-as.character(nnn[,9])
for (pop in  pops){ # pop="KZL"
  ip<-which(pops %in% paste(pop))
  a<-read.table(paste(pop,".table.repol.txt",sep=""), header =F)
  a[ a == "-9" ] <- NA
  a$na<-apply(is.na(a), 1, sum)
  a$an<-((ncol(a)-10)-a$na)*as.numeric(a[1,2])
  a$af<-rowSums(a[,10:(ncol(a)-2)],na.rm = T) / a$an
  mtot[,ip+4]<-a$af}
write.table(mtot,"AFsperLin.txt",append = F,quote = F,row.names = F,col.names = T)
##repolarize
k<-read.table("AFsperLin.txt",h=T)
k$anc<-substr(x = k$AAS,start = 1,stop = 3)
k$der<-str_sub(string = k$AAS, start= -3)
k$aas_pos<-as.character(parse_number(as.character(k$AAS)))
k1<-select(.data = k,c1)
k2<-select(.data = k,c2)
k$c1<-rowSums(k1,na.rm = T)/ncol(k1)
k$c2<-rowSums(k2,na.rm = T)/ncol(k2)
#e<-0
for (i in  1:nrow(k)){ # i=1
  if (k$c1[i]>k$c2[i])
  {k[i,5:(ncol(k)-10)] <- 1 - k[i,5:(ncol(k)-10)]
  a<-k$anc[i]
  k$anc[i]<-k$der[i]
  k$der[i]<-a
  #e<-e+1
  #print(k[i,1:2])
  } else {}}
#print(e)
write.table(k[1:(ncol(k)-2)],"AFsperLinRepol.txt",append = F,quote = F,row.names = F,col.names = T)
#calculate haplotype allele frequencies
k<-read.table("AFsperLinRepol.txt",h=T)
g<-read.table("PopGeogr.txt",h=T)
for (al in unique(k$alcode)){ # al<-"AL1G10680"
  a1<-subset(k,k$alcode %in% al)
  if (nrow(a1)==1)
  {} else {
    a2<-as.data.frame(t(a1[,5:(ncol(a1)-3)]))
    #neutral
    g$c1<-1-(apply(a2, 1, max))
    #selected
    g$c2<-apply(a2, 1, min)
    #recombined
    g$oth<-abs(1-g$c1-g$c2)
    n<-read.table("ALmeiotic.ALcode.tab.Name.InclZyp.txt")
    name<-as.matrix(subset(n,n$V1 %in% al))[1,2]
    write.table(g,paste(al,"_",name,".hapl.txt",sep=""),quote = F,sep = "\t",row.names = F)}}
#MAP
library(RgoogleMaps)
library(argosfilter)
library(mapplots)
list<-list.files(path = ".", pattern = "*.hapl.txt")
for (l in list){
  df<-read.table(paste(l), header=T)
  lat <- c(43,55.8) #define our map's ylim
  lon <- c(4.5,27) #define our map's xlim
  name<-strsplit(l, split='[.]')[[1]][1]
  prot<-strsplit(name, split='[_]')[[1]][2]
  terrmap<-GetMap.bbox(lon,lat,maptype= "terrain", destfile = "terrain.png",path = paste0("&style=feature:road|element:all|visibility:off","&style=feature:administrative.country|element:geometry.stroke|visibility:off","&style=feature:all|element:labels|visibility:off"), MINIMUMSIZE = T)
  size<-terrmap$size
  pdf(file=paste("figures/",name,".hapl.pdf", sep=""),width = size[1]/72, height = size[2]/72)
  PlotOnStaticMap(terrmap)
  df<-na.omit(df)
  for (i in 1:nrow(df)){
    d1<-distance(df$lat[i], df$lat[i], lon[1], lon[2]) 
    d2<-distance(df$lat[i], df$lat[i], lon[1], df$lon[i])
    a1<-size[1]/d1
    a2<-d2*a1
    a3<- -(size[1]/2)+a2+df$offsetx[i]
    b1<-size[2]/(lat[2]-lat[1])
    b2<-(df$lat[i]-lat[1])*b1
    b3<- -(size[2]/2)+b2+df$offsety[i]
    terrmap<-add.pie(z=c(df$c1[i],df$c2[i],df$oth[i]), x=a3, y=b3, radius=10, col = colors, labels=NA)
    terrmap<-text( x=a3, y=b3, labels = df$pop[i],pos = df$pos[i],col = as.character(df$col[i]), offset = 0.83,cex = 0.7)}
  terrmap<-legend("topleft", inset=.005, title=paste(prot," haplotypes",sep=""),contr,fill = colors, horiz=F, cex=1,bg = 'white')
  dev.off()}






##3hapl
  terrmap<-add.pie(z=c(df$pan[i],df$dip[i],df$tet[i],df$oth[i]), x=a3, y=b3, radius=10, col = c("orange","red3", "blue2", "grey"), labels=NA)
  terrmap<-text( x=a3, y=b3, labels = df$pop[i],pos = df$pos[i],col = as.character(df$col[i]), offset = 0.83,cex = 0.7)}
terrmap<-legend("topleft", inset=.005, title="PMS1 haplotypes",c("pannonian","diploid","tetraploid","other"),fill = c("orange","red3", "blue2", "grey"), horiz=F, cex=1,bg = 'white')
dev.off()




###BARPLOT
setwd("/home/majda/Desktop/ownCloud/JICAutumn2016/haplotypes/freq")
##############4X
pdf(file="TetraploidHaplotypes.pdf",width = 11.69, height = 8.27,title = "Tetraploid haplotypes",pointsize = 15)
par(mai=c(0,0.15,0,0))
par(mfrow=c(12,1))
plot.new()
prot<-as.data.frame(c('HEI10', 'SMC1','SCC4', 'DYAD', 'PRD3','SDS','ASY1','ASY3','PDS5B','SYN1'))
for(i in 1:nrow(prot)){ 
  name=as.character(paste(as.character(prot[i,1]),'hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[,3:6]
  rownames(d1) <- d0[,2]
  barplot(t(as.matrix(d1)),col = c("orange","red3", "blue2", "grey"),border="NA", las=2,xlab = NA,ylab=NA, axes = F, axisnames = T,space=0)
  title(ylab=title, line=-1.8, cex.lab=1.2)}
dev.off()
#############2X
pdf(file="DiploidHaplotypes.pdf",width = 11.69, height = 8.27,title = "Diploid haplotypes",pointsize = 15)
par(mai=c(0,0.15,0,0))
par(mfrow=c(12,1))
plot.new()
prot<-as.data.frame(c('BRCA2B', 'PDS5E','PDS5D', 'RAD54', 'RAD17','RAD50'))
for(i in 1:nrow(prot)){ 
  name=as.character(paste(as.character(prot[i,1]),'hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[,3:5]
  rownames(d1) <- d0[,2]
  barplot(t(as.matrix(d1)),col = c("darkslateblue", "orange", "grey"),border="NA", las=2,xlab = NA,ylab=NA, axes = F, axisnames = T,space=0)
  title(ylab=title, line=-1.8, cex.lab=1.2)}
prot<-as.data.frame(c('SMC2-1'))
for(i in 1:nrow(prot)){ 
  name=as.character(paste(as.character(prot[i,1]),'hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[,3:5]
  rownames(d1) <- d0[,2]
  barplot(t(as.matrix(d1)),col = c("darkslateblue", "green4", "grey"), border="NA", las=2,xlab = NA,ylab=NA, axes = F, axisnames = T,space=0)
  title(ylab=title, line=-1.8, cex.lab=1.2)}
prot<-as.data.frame(c('PMS1','SMC6B'))
for(i in 1:nrow(prot)){ 
  name=as.character(paste(as.character(prot[i,1]),'hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[,3:5]
  rownames(d1) <- d0[,2]
  barplot(t(as.matrix(d1)),col = c("darkslateblue", "steelblue1", "grey"),border="NA", las=2,xlab = NA,ylab=NA, axes = F, axisnames = T,space=0)
  title(ylab=title, line=-1.8, cex.lab=1.2)}
prot<-as.data.frame(c('SMG7'))
for(i in 1:nrow(prot)){ 
  name=as.character(paste(as.character(prot[i,1]),'hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[,3:5]
  rownames(d1) <- d0[,2]
  barplot(t(as.matrix(d1)),col = c("black","darkslateblue", "grey"),border="NA", las=2,xlab = NA,ylab=NA, axes = F, axisnames = T,space=0)
  title(ylab=title, line=-1.8, cex.lab=1.2)}
dev.off()

#########ALL
pdf(file="AllHaplotypes.pdf",width = 11.69, height = 8.27*2,title = "All haplotypes",pointsize = 18)
par(mai=c(0,0.15,0,0))
par(mfrow=c(23,1))
plot.new()
prot<-as.data.frame(c('HEI10', 'SMC1','SCC4', 'DYAD', 'PRD3','SDS','ASY1','ASY3','PDS5B','SYN1'))
for(i in 1:nrow(prot)){ 
  name=as.character(paste(as.character(prot[i,1]),'hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[,3:6]
  rownames(d1) <- d0[,2]
  barplot(t(as.matrix(d1)),col = c("orange","red3", "blue2", "grey"),border="NA", las=2,xlab = NA,ylab=NA, axes = F, axisnames = T,space=0)
  title(ylab=title, line=-1.8, cex.lab=1.2)}
prot<-as.data.frame(c('BRCA2B', 'PDS5E','PDS5D', 'RAD54','RAD50'))
for(i in 1:nrow(prot)){ 
  name=as.character(paste(as.character(prot[i,1]),'hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[,3:5]
  rownames(d1) <- d0[,2]
  barplot(t(as.matrix(d1)),col = c("darkslateblue", "orange", "grey"),border="NA", las=2,xlab = NA,ylab=NA, axes = F, axisnames = T,space=0)
  title(ylab=title, line=-1.8, cex.lab=1.2)}
prot<-as.data.frame(c('SMC2-1'))
for(i in 1:nrow(prot)){ 
  name=as.character(paste(as.character(prot[i,1]),'hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[,3:5]
  rownames(d1) <- d0[,2]
  barplot(t(as.matrix(d1)),col = c("darkslateblue", "green4", "grey"), border="NA", las=2,xlab = NA,ylab=NA, axes = F, axisnames = T,space=0)
  title(ylab=title, line=-1.8, cex.lab=1.2)}
prot<-as.data.frame(c('PMS1','SMC6B'))
for(i in 1:nrow(prot)){ 
  name=as.character(paste(as.character(prot[i,1]),'hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[,3:5]
  rownames(d1) <- d0[,2]
  barplot(t(as.matrix(d1)),col = c("darkslateblue", "steelblue1", "grey"),border="NA", las=2,xlab = NA,ylab=NA, axes = F, axisnames = T,space=0)
  title(ylab=title, line=-1.8, cex.lab=1.2)}
prot<-as.data.frame(c('SMG7'))
for(i in 1:nrow(prot)){ 
  name=as.character(paste(as.character(prot[i,1]),'hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[,3:5]
  rownames(d1) <- d0[,2]
  barplot(t(as.matrix(d1)),col = c("black","darkslateblue", "grey"),border="NA", las=2,xlab = NA,ylab=NA, axes = F, axisnames = T,space=0)
  title(ylab=title, line=-1.8, cex.lab=1.2)}
dev.off()
#

###Overlapping datasets
library(stringr)
library(VennDiagram)
setwd("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/")
#inF<-"FstSubsampled/TETDIP/AFDmin0.1/"
inF<-"FAAD_overlaps/FstH/TETDIP/"
inG<-"FAAD_overlaps/Grantham/TETDIP/"
inS<-"FAAD_overlaps/SIFT/TETDIP/"
mat<-matrix(nrow = 7,ncol = 2,dimnames = list(c("G","S","F","FSG","GS","GF","SF"),c("0.1%","1%")))
g1<-read.table(paste(inG,"top_0.01_outliers_grantham_dip_tet.2.5.scores",sep=""))
ss<-read.table(paste(inS,"2.5_outliers_0.01.scores.sorted.reformat.txt",sep=""),sep = "\t",fill = T)
s1<-as.data.frame(do.call(rbind, str_split(ss$V1, ' ')))
f1<-read.table(paste(inF,"quantile0.99.intervals",sep=""),sep = ":",h=F)
g<-as.character(interaction(g1$V1, g1$V2))
s<-as.character(interaction(s1$V1, s1$V2))
f<-as.character(interaction(f1$V1, f1$V2))
mat[1,1]<-length(g)
mat[2,1]<-length(s)
mat[3,1]<-length(f)
mat[4,1]<-length(Reduce(intersect, list(g,s,f)))
mat[5,1]<-length(Reduce(intersect, list(g,s)))
mat[6,1]<-length(Reduce(intersect, list(g,f)))
mat[7,1]<-length(Reduce(intersect, list(s,f)))

mat[1,2]<-length(g)
mat[2,2]<-length(s)
mat[3,2]<-length(f)
mat[4,2]<-length(Reduce(intersect, list(g,s,f)))
mat[5,2]<-length(Reduce(intersect, list(g,s)))
mat[6,2]<-length(Reduce(intersect, list(g,f)))
mat[7,2]<-length(Reduce(intersect, list(s,f)))
write.table(Reduce(intersect, list(g,s,f)), "FAAD_overlaps/overlaps0.99.genomeWide.FSG.intervals",quote = F,col.names = F,row.names = F)
v<-venn.diagram(x=list("Fst_H"=f,"FAAD_S"=s,"FAAD_G"=g),"FAAD_overlaps/venn0_999.tiff", lty = "blank", fill = c("blue1","red2","green"),alpha=0.3 , main = "Overlaps of different tests",cex=1, cat.cex=1,sub = "0.1% outliers",height = 2000,width = 2500)
v<-venn.diagram(x=list("Fst_H"=f,"FAAD_S"=s,"FAAD_G"=g),"FAAD_overlaps/venn0_99.tiff", lty = "blank", fill = c("blue1","red2","green"),alpha=0.3 , main = "Overlaps of different tests",cex=1, cat.cex=1,sub = "1% outliers",height = 2000,width = 2500)
#combinations
v<-venn.diagram(x=list("Fst_H 1%"=f,"FAAD_S 0.1%"=s,"FAAD_G 0.1%"=g),"FAAD_overlaps/venn99_999_999.tiff", lty = "blank", fill = c("blue1","red2","green"),alpha=0.3 , main = "Overlaps of different tests",cex=1, cat.cex=1,sub = "1% outliers",height = 2000,width = 2500)
v<-venn.diagram(x=list("Fst_H 1%"=f,"FAAD_S 1%"=s,"FAAD_G 0.1%"=g),"FAAD_overlaps/venn99_99_999.tiff", lty = "blank", fill = c("blue1","red2","green"),alpha=0.3 , main = "Overlaps of different tests",cex=1, cat.cex=1,sub = "1% outliers",height = 2000,width = 2500)
write.table(mat,file = "FAAD_overlaps/overlaps.txt",quote = F,sep = "\t",row.names = T,col.names = T)
#ANNOTATE overlaps
#setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/")
library(data.table)
library(stringr)
library(readr)
quant=0.99
a<-read.table("../lists/pI_hyd.txt",header = T)
n<-read.table("../lists/ALmeiotic.ALcode.tab.Name.txt",h=F)
i<-read.table("../pai_Malvidae/Alignment_Identities_allSites.txt",h=F)
dd<-read.table(paste("meiotic.highFstFinalSummary.",quant,".txt",sep=""),header = F,col.names = c("scaf","AN1","AN2","AC1","AC2","alcode","aas","ann","pos","TETDIP"))
o<-read.table("FAAD_overlaps/overlaps0.99.genomeWide.FSG.intervals",h = F,sep = ".")
dd$pos<-dd$pos+1
#get only overlaps
d<-subset(dd,dd$pos %in% o$V2 & dd$scaf %in% o$V1)
#d$V2 <- d$V3 <- d$V4 <- d$V5 <- NULL
d$anc<-substr(x = d$aas,start = 1,stop = 3)
d$der<-str_sub(string = d$aas, start= -3)
indAnc<-which(colnames(d) %in% "anc")
indDer<-which(colnames(d) %in% "der")
d$anc_pI<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,4])
d$der_pI<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,4])
d$delta_pI<-abs(d$anc_pI-d$der_pI)
#hist(d$delta_pI,breaks = 50)
d$anc_hydCA<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,2])
d$der_hydCA<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,2])
d$delta_hydCA<-abs(d$anc_hydCA-d$der_hydCA)
#hist(d$delta_hydCA,breaks = 50)
d$anc_hydCC<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,3])
d$der_hydCC<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,3])
d$delta_hydCC<-abs(d$anc_hydCC-d$der_hydCC)
#hist(d$delta_hydCC,breaks = 50)
indAl<-which(colnames(d) %in% "alcode")
d$prot<-apply(X = d,1, function(x) subset(x = n, subset = n$V1 %in% x[indAl])[1,2])
d$aas_pos<-as.character(parse_number(as.character(d$aas)))
indAP<-which(colnames(d) %in% "aas_pos")
d$al_ident<-apply(X = d,1, function(x) subset(x = i, subset = i$V6 %in% x[indAl] & as.character(i$V1) %in% as.character(x[indAP]))[1,3])
d$al_hydr<-apply(X = d,1, function(x) subset(x = i, subset = i$V6 %in% x[indAl] & as.character(i$V1) %in% as.character(x[indAP]))[1,4])
d$al_pI<-apply(X = d,1, function(x) subset(x = i, subset = i$V6 %in% x[indAl] & as.character(i$V1) %in% as.character(x[indAP]))[1,5])
g2<-g1[ order(g1[,1], g1[,2]), ]
d$grantham<-apply(X = d,1, function(x) subset(x = g2, subset = g2$V1 %in% x[1] & as.numeric(g2$V2) %in% as.numeric(x[9]))[1,8])
d$FAAD_dip<-apply(X = d,1, function(x) subset(x = g2, subset = g2$V1 %in% x[1] & as.numeric(g2$V2) %in% as.numeric(x[9]))[1,9])
d$FAAD_tet<-apply(X = d,1, function(x) subset(x = g2, subset = g2$V1 %in% x[1] & as.numeric(g2$V2) %in% as.numeric(x[9]))[1,10])
d$delta_FAAD<-apply(X = d,1, function(x) subset(x = g2, subset = g2$V1 %in% x[1] & as.numeric(g2$V2) %in% as.numeric(x[9]))[1,11])
d$SIFT<-apply(X = d,1, function(x) subset(x = s1, subset = s1$V1 %in% x[1] & as.numeric(as.character(s1$V2)) %in% as.numeric(x[9]))[1,8])
d$SIFTann<-apply(X = d,1, function(x) subset(x = s1, subset = s1$V1 %in% x[1] & as.numeric(as.character(s1$V2)) %in% as.numeric(x[9]))[1,12])
colnames(d)[10]<-"Fst_H"
write.table(d, paste("FAAD_overlaps/overlaps.FinalSummary.Ann",quant,".txt",sep=""),append = F,quote = F,sep = "\t",row.names = F)
#





























#MKT with CRO outgroup - for between species comparison - useful to campare arenosa/lyrata?
#To minimize the impact of slightly deleterious mutations, it has been proposed to exclude polymorphisms that are below a certain cutoff frequency, such as <8% or <5%
#We have suggested that polymorphisms below 15% be excluded from MK-type analyses when there is evidence that some nonsilent polymorphisms are slightly deleterious.
#removing low-frequency polymorphisms increases the variance of the estimate of a, but this is a price that has to be paid to obtain a less biased estimate. Hopefully, by providing a recommended cutoff frequency, we will remove the temptation to search for the frequency that yields the highest value of a because this is statistically difficult to defend.
###try with different frequency cutoff
#1.load contrast + both to outgroup - podminky z ###per meiotic genes
#2. calc for each gene - na,0, number Pn, Dn, Ps, Ds, alfa, NI
#3. Play with meiotic vs. the rest of genome - permutation? take all?
library(data.table)
setwd("~/JICAutumn2016/finalAnalysis29Apr/")
con<-c("TETDIP","PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA","CROPAN","CRODIN","CROSEC","CROBAL","CROWCA")

con<-c("TETDIP")
al <- fread("lists/ALcodes.txt",header = F)
#to filter for missing data: 0.9 = 90% nomissing
mffg<-0.49
#to filter for minimum derived allele frequency (0.15 - recommended)
mdafP<-0.15
##pridat az budu vedet GW treshold!!
mdafD<-0.5
quant<-0.95
for (contr in  con){
  system(command = paste("grep -E 'missense_variant|synonymous_variant' data/",contr,"_WS1_MS1_BPM.txt > data/",contr,".synNon_WS1_MS1_BPM.txt",sep=""))
  #  contr<-"DINTET"
  mat <- matrix(nrow = nrow(al), ncol = 9,dimnames = list(c(),c("contr","gene","pN","pS","dN","dS","alfa","NI","dos")))
  pop0<-substr(contr,1,3)
  pop1<-substr(contr,4,6)
  p<-fread(paste("head -n -1 data/",pop0,pop1,"_WS1_MS1_BPM.txt",sep = ""),header = T)
  # Histogram Colored (blue and red)
  png(filename = paste("results/AFS_", contr,"_",mdafP,"_",mdafD,".png"),width = 1300,height = 850,pointsize = 24)
  contrN<-subset(x = p, p$AAS %in% "missense_variant")
  contrS<-subset(x = p, p$AAS %in% "synonymous_variant")
  hist(contrN$AC0/contrN$AN0, col=rgb(1,0.7,0,0.6),xlim=c(0,1), xlab="AF",breaks = 100,main = paste("AFS in ",contr," pair",sep=""))
  hist(contrS$AC0/contrS$AN0, col=rgb(0,0,1,0.6), add=T,breaks = 100)
  legend(legend = c("synonymous", "nonsynonymous"),fill = c("blue", "orange"),x = "topright")
  abline(v = mdafP,col="red")
  dev.off()
  d0<-fread(paste("head -n -1 data/CRO",pop0,"_WS1_MS1_BPM.txt",sep = ""),header = T)
  d1<-fread(paste("head -n -1 data/CRO",pop1,"_WS1_MS1_BPM.txt",sep = ""),header = T)
  dd<-merge(x=d0,y=d1, by=c("scaff","end"))
  ###TADY NUTNO VYBRAT VARIANTU, ZALEZI NA TYPU ANALYZY!!!
  #filter missing data and minDAF (possibly abs((p$AC0+p$AC1)/(p$AN0+p$AN1)) > mdafP)
  #p1<-subset(x = p,subset = p$AN0 > max(p$AN0)*mffg & p$AN1 > max(p$AN1)*mffg & abs(p$AFD) >= quantile(abs(p$AFD),probs = quant))
  p1<-subset(x = p,subset = p$AN0 > max(p$AN0)*mffg & p$AN1 > max(p$AN1)*mffg & abs(p$AFD) >= mdafP)
  #filter missing data and almost fixed variants (possibly: abs((dd$AC0.x/dd$AN0.x) - ((dd$AC1.x+dd$AC1.y)/(dd$AN1.x+dd$AN1.y))))
  #dd1<-subset(x = dd,subset = (dd$AN0.x+dd$AN0.y) > max((dd$AN0.x+dd$AN0.y),na.rm = T)*mffg & (dd$AN1.x+dd$AN1.y) > max(dd$AN1.x+dd$AN1.y,na.rm = T)*mffg & abs((dd$AC0.x/dd$AN0.x)-((dd$AC1.x+dd$AC1.y)/(dd$AN1.x+dd$AN1.y))) >= quantile(abs((dd$AC0.x/dd$AN0.x)-((dd$AC1.x+dd$AC1.y)/(dd$AN1.x+dd$AN1.y))),probs = quant))
  dd1<-subset(x = dd,subset = (dd$AN0.x+dd$AN0.y) > max((dd$AN0.x+dd$AN0.y),na.rm = T)*mffg & (dd$AN1.x+dd$AN1.y) > max(dd$AN1.x+dd$AN1.y,na.rm = T)*mffg & abs((dd$AC0.x/dd$AN0.x)-((dd$AC1.x+dd$AC1.y)/(dd$AN1.x+dd$AN1.y))) >= mdafD)
  
  for (g in readLines("lists/ALcodes.txt"))
  {# g='AL4G46460'
    index<-which(al$V1 %in% paste(g))
    pN<-nrow(subset(x = p1,subset = p1$ALcode %in% g & p1$AAS %in% "missense_variant"))
    pS<-nrow(subset(x = p1,subset = p1$ALcode %in% g & p1$AAS %in% "synonymous_variant"))
    dN<-nrow(subset(x = dd1,subset = dd1$ALcode.x %in% g & dd1$AAS.x %in% "missense_variant"))
    dS<-nrow(subset(x = dd1,subset = dd1$ALcode.x %in% g & dd1$AAS.x %in% "synonymous_variant"))
    if (pS>0 & dN>0) 
    {alfa=1-(dS*pN)/(dN*pS)
    NI=(pN/pS)/(dN/dS)} else 
    {alfa<-"NA"
    NI<-"NA"}
    
    if (pS+pN>0 & dS+dN>0) 
    {dos=dN/(dN+dS)-pN/(pN+pS)} else 
    {dos<-"NA"}
    mat[as.numeric(index),1]<-contr
    mat[as.numeric(index),2]<-g
    mat[as.numeric(index),3]<-pN
    mat[as.numeric(index),4]<-pS
    mat[as.numeric(index),5]<-dN
    mat[as.numeric(index),6]<-dS
    mat[as.numeric(index),7]<-alfa #to estimate the fraction (α) of substitutions driven to fixation by positive selection at the functional sites
    mat[as.numeric(index),8]<-NI
    mat[as.numeric(index),9]<-dos}
  write.table(mat,append = F,file = paste("results/",contr,"_mafdP_",mdafP,"_mafdD_",mdafD,"_MKT.txt",sep=""),quote = F, sep = "\t",col.names = T,row.names = F)}


##################For HFD
library(data.table)
setwd("/home/majda/Desktop/protFeatures/data/")
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
  
  png(filename=paste(contr,"_AFDperNandSHist.png", sep=""), 
      units = "px",
      width=1800, 
      height=1000, 
      pointsize=16)
  hist(abs(s$AFD), freq = F, breaks = 200, xlab =paste("abs(AFD) - ",contr, sep=""), col=rgb(1,0,0,0.5), main =NA,xlim=c(0,1))
  hist(abs(n$AFD), freq=F, breaks = 200, col=rgb(0,0,1,0.5), add=T,xlim=c(0,1))
  abline(v=(quantile(abs(s$AFD),probs=c(.99))),col="red")
  legend("topleft", c("synonymous","nonsynonymous"), fill=c( "red", "blue"))
  box()
  dev.off()}
write.table(msum,append = F,file ="AFSStatsPerContrast.txt",quote = F, sep = "\t",col.names = T,row.names = T)

#########################################################################not needed anymore
###per whole-genome - pro ScanTools pridat do predchozi smycky
for (lin in  c('CRO','PAN','DIN','BAL','SEC','WCA','TET')){
  # lin="CRO"
  inp<-fread(paste(lin,".table.repol.txt",sep=""),header=F)
  inp[inp=="-9"]<-NA
  mat1 <- matrix(nrow = 1, ncol = 6)
  n<-sum(subset(inp,inp$V8 %in% "missense_variant")[,10:ncol(inp)],na.rm = T)
  s<-sum(subset(inp,inp$V8 %in% "synonymous_variant")[,10:ncol(inp)],na.rm = T)
  ns<-as.numeric(n)/as.numeric(s)
  mat1[1,1]<-as.character(inp[1,1])
  mat1[1,2]<-"genome"
  mat1[1,3]<-n
  mat1[1,4]<-s
  mat1[1,5]<-ns
  mat1[1,6]<-"1"
  write.table(mat1,append = T,file = "NSgenomeWide.AllLin.txt",quote = F, sep = "\t",col.names = F,row.names = F)
}

#For ScanTools
########N/S analysis per gene

lineages<-c('+lineages+')
for (lin in  lineages){
  system(command = paste("grep -E 'missense_variant|synonymous_variant' ",lin,".table.repol.txt > ",lin,".synNon.txt",sep=""))
  inp<-read.table(paste(lin,".synNon.txt",sep=""),header=F)
  inp[inp=="-9"]<-NA
  for (g in readLines("'+alcode_path.split("/")[-1]+'"))
  { mat <- matrix(nrow = 1, ncol = 6)
    gene<-subset(x = inp,subset = inp$V7 %in% g)
    if (nrow(subset(gene,gene$V8 %in% "missense_variant"))==0) {n<-0} else 
    {n<-sum(subset(gene,gene$V8 %in% "missense_variant")[,10:ncol(gene)],na.rm = T)}
    if (nrow(subset(gene,gene$V8 %in% "synonymous_variant"))==0) {s<-0} else 
    {s<-sum(subset(gene,gene$V8 %in% "synonymous_variant")[,10:ncol(gene)],na.rm = T)}
    if ((n+s)<5) {check<-NA} else 
    {check<-1}
    ns<-as.numeric(n)/as.numeric(s)
    mat[1,1]<-as.character(gene[1,1])
    mat[1,2]<-g
    mat[1,3]<-n
    mat[1,4]<-s
    mat[1,5]<-ns
    mat[1,6]<-check
    write.table(mat,append = T,file = "NSperGene.AllLin.Allgenes.txt",quote = F, sep = "\t",col.names = F,row.names = F)}
  mat1 <- matrix(nrow = 1, ncol = 6)
  n<-sum(subset(inp,inp$V8 %in% "missense_variant")[,10:ncol(inp)],na.rm = T)
  s<-sum(subset(inp,inp$V8 %in% "synonymous_variant")[,10:ncol(inp)],na.rm = T)
  ns<-as.numeric(n)/as.numeric(s)
  mat1[1,1]<-as.character(inp[1,1])
  mat1[1,2]<-"genome"
  mat1[1,3]<-n
  mat1[1,4]<-s
  mat1[1,5]<-ns
  mat1[1,6]<-"1"
  write.table(mat1,append = T,file = "NSgenomeWide.AllLin.txt",quote = F, sep = "\t",col.names = F,row.names = F)}


lineages<-c("CRO")
for (lin in  lineages){
  inp<-read.table(paste(lin,".synNon.txt",sep=""),header=F)
  inp[inp=="-9"]<-NA
  for (g in readLines("ALcodesAll.txt"))
  { mat <- matrix(nrow = 1, ncol = 6)
  gene<-subset(x = inp,subset = inp$V7 %in% g)
  if (nrow(subset(gene,gene$V8 %in% "missense_variant"))==0) {n<-0} else 
  {n<-sum(subset(gene,gene$V8 %in% "missense_variant")[,10:ncol(gene)],na.rm = T)}
  if (nrow(subset(gene,gene$V8 %in% "synonymous_variant"))==0) {s<-0} else 
  {s<-sum(subset(gene,gene$V8 %in% "synonymous_variant")[,10:ncol(gene)],na.rm = T)}
  if ((n+s)<5) {check<-NA} else 
  {check<-1}
  ns<-as.numeric(n)/as.numeric(s)
  mat[1,1]<-as.character(gene[1,1])
  mat[1,2]<-g
  mat[1,3]<-n
  mat[1,4]<-s
  mat[1,5]<-ns
  mat[1,6]<-check
  write.table(mat,append = T,file = "NSperGene.AllLin.Allgenes.txt",quote = F, sep = "        ",col.names = F,row.names = F)}
  mat1 <- matrix(nrow = 1, ncol = 6)
  n<-sum(subset(inp,inp$V8 %in% "missense_variant")[,10:ncol(inp)],na.rm = T)
  s<-sum(subset(inp,inp$V8 %in% "synonymous_variant")[,10:ncol(inp)],na.rm = T)
  ns<-as.numeric(n)/as.numeric(s)
  mat1[1,1]<-as.character(inp[1,1])
  mat1[1,2]<-"genome"
  mat1[1,3]<-n
  mat1[1,4]<-s
  mat1[1,5]<-ns
  mat1[1,6]<-"1"
  write.table(mat1,append = T,file = "NSgenomeWide.AllLin.txt",quote = F, sep = "       ",col.names = F,row.names = F)}

#For ScanTools
########Fst

con<-c("TETDIP","PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA","CROPAN","CRODIN","CROSEC","CROBAL","CROWCA")
msum<-matrix(nrow = sum(table(con)), ncol = 16,dimnames =list(c(con),c("Smin","S1Q","Smedian","Smean","S3Q","Smax","Ssd", "Nmin","N1Q","Nmedian","Nmean","N3Q","Nmax","Nsd","Squantile95","Squantile99")))
al <- read.table("lists/ALcodesAll.txt",header = F)
mat <- matrix(nrow = nrow(al), ncol = sum(table(con)),dimnames = list(c(al$V1),c(con))) 
for (contr in  con){
  # contr="CROWCA"
  indlin<-which(con %in% paste(contr))
  inp<-head(read.table(paste("data/",contr,"_WS1_MS1_BPM.txt",sep=""),header=T,fill = T),-1)
  s<-subset(inp,inp$ann %in% "synonymous_variant") ###ann
  n<-subset(inp,inp$ann %in% "missense_variant")
  sq <- quantile(s$FstH,probs=c(0.999))
  msum[as.numeric(indlin),1:6]<-summary (s$FstH)
  msum[as.numeric(indlin),7]<-sd(s$FstH)
  msum[as.numeric(indlin),8:13]<-summary (n$FstH)
  msum[as.numeric(indlin),14]<-sd(n$FstH)
  #msum[as.numeric(indlin),15] <-quantile(s$FstH,probs=c(.95)) ###
  msum[as.numeric(indlin),15]<-sq ###
  for (g in readLines("lists/ALcodesAll.txt"))
  { index<-which(al$V1 %in% paste(g))
  gene<-subset(x = inp,subset = inp$ALcode %in% g & inp$ann %in% "missense_variant") ###
  
  if (nrow(gene)==0) {high<-0} else 
  {high<-nrow(subset(gene,gene$FstH >= as.numeric(sq)))}
  mat[index,indlin]<-high
  }
}
write.table(mat,append = F,file="FstHighperGene.AllLin.Allgenes.0.95.txt",quote = F, sep = "	", col.names = T,row.names = T)
write.table(msum,append = F,file ="FstStatsPerContrast.0.95.txt",quote = F, sep = "	", col.names = T,row.names = T)
########
con<-"TETDIP"
msum<-matrix(nrow = 1, ncol = 15,dimnames =list(c(con),c("Smin","S1Q","Smedian","Smean","S3Q","Smax","Ssd", "Nmin","N1Q","Nmedian","Nmean","N3Q","Nmax","Nsd","Squantile_0.995")))
al <- read.table("ALcodesAll.txt",header = F)
mat <- matrix(nrow = nrow(al), ncol = 1,dimnames = list(as.character((al[,1])),c(con)))
inp<-head(read.table(paste("data/",con,".synNon_WS1_MS1_BPM.txt",sep=""),header=T,fill = T),-1)
s<-subset(inp,inp$ann %in% "synonymous_variant" & inp$AFD >= 0.1)
n<-subset(inp,inp$ann %in% "missense_variant")
sq <- quantile(s$FstH,probs=c(0.99))
msum[1,1:6]<-summary (s$FstH)
msum[1,7]<-sd(s$FstH)
msum[1,8:13]<-summary (n$FstH)
msum[1,14]<-sd(n$FstH)
msum[1,15]<-sq
for (g in readLines("ALcodesAll.txt"))
{ index<-which(al$V1 %in% paste(g))
gene<-subset(x = inp,subset = inp$ALcode %in% g & inp$ann %in% "missense_variant")

if (nrow(gene)==0) {high<-0} else 
{high<-nrow(subset(gene,gene$FstH >= as.numeric(sq)))}
mat[index,1]<-high
}
write.table(mat,append = F,file="FstHighperGene.Allgenes.0.995.TETDIP.txt",quote = F, sep = "   ", col.names = T,row.names = T)
write.table(msum,append = F,file ="FstStatsPerContrast.0.995.TETDIP.txt",quote = F, sep = "     ", col.names = T,row.names = T)

#list of identities old
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/")
library(data.table)
library(stringr)
d<-read.table("lists/MeioticCSVIdentity.txt",header = F)
for (g in d$V1){
  dd<-read.table(paste("lists/",g,sep=""),h=T,sep = ",")
  gg<-substr(g,1,9)
  dd$V3<-gg
  write.table(dd[1:nrow(dd),],"Alignment_Identities_allSites.txt",append = T,quote = F,row.names = F,col.names = F,sep = "\t")}
