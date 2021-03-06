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
con<-c("TETDIP")
# con<-c("CROPAN")
al <- fread("lists/ALcodes.txt",header = F)
#to filter for missing data: 0.9 = 90% nomissing
mffg<-0.49
#to filter for minimum derived allele frequency (0.15 - recommended)
mdafP<-0.15
##pridat az budu vedet GW treshold!!
mdafD<-0.5260489 
for (contr in  con){
  system(command = paste("head -n 1 data/",contr,"_WS1_MS1_BPM.txt > data/",contr,".synNon_WS1_MS1_BPM.txt",sep=""))
  system(command = paste("grep -E 'missense_variant|synonymous_variant' data/",contr,"_WS1_MS1_BPM.txt >> data/",contr,".synNon_WS1_MS1_BPM.txt",sep=""))
  #  contr<-"DINTET"
  mat <- matrix(nrow = nrow(al), ncol = 11,dimnames = list(c(),c("contr","gene","pN","pS","dN","dS","alfa","NI","DoS","p-value","name")))
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
  summary(abs(p1$AFD))
  quantile(abs(p1$AFD),0.99)
  
  for (gene in readLines("lists/ALmeiotic.ALcode.tab.Name.txt"))
  {# gene='AL1G10680 PRD3'
    g<-substr(gene,1,9)
    name<-substr(gene,11,30)
    index<-which(al$V1 %in% paste(g))
    pN<-nrow(subset(x = p1,subset = p1$ALcode %in% g & p1$ann %in% "missense_variant"))
    pS<-nrow(subset(x = p1,subset = p1$ALcode %in% g & p1$ann %in% "synonymous_variant"))
    dN<-nrow(subset(x = d1,subset = d1$ALcode %in% g & d1$ann %in% "missense_variant"))
    dS<-nrow(subset(x = d1,subset = d1$ALcode %in% g & d1$ann %in% "synonymous_variant"))
    pval<-fisher.test(matrix(c(pN, pS, dN, dS), nrow=2),alternative = "less")$p.value
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
    mat[as.numeric(index),10]<-pval
    mat[as.numeric(index),11]<-name}
  write.table(mat,append = F,file = paste("results/",contr,"_mafdP_",mdafP,"_mafdD_",mdafD,"_MKT.txt",sep=""),quote = F, sep = "\t",col.names = T,row.names = F)}
# Genome
mat <- matrix(nrow = 1, ncol = 11,dimnames = list(c(),c()))
pN<-nrow(subset(x = p1,p1$ann %in% "missense_variant"))
pS<-nrow(subset(x = p1,p1$ann %in% "synonymous_variant"))
dN<-nrow(subset(x = d1,d1$ann %in% "missense_variant"))
dS<-nrow(subset(x = d1,d1$ann %in% "synonymous_variant"))
pval<-fisher.test(matrix(c(pN, pS, dN, dS), nrow=2),alternative = "less")$p.value
if (pS>0 & dN>0) 
{alfa=1-(dS*pN)/(dN*pS)
NI=(pN/pS)/(dN/dS)} else 
{alfa<-"NA"
NI<-"NA"}
if (pS+pN>0 & dS+dN>0) 
{dos=dN/(dN+dS)-pN/(pN+pS)} else 
{dos<-"NA"}
mat[1,1]<-contr
mat[1,2]<-"Genome"
mat[1,3]<-pN
mat[1,4]<-pS
mat[1,5]<-dN
mat[1,6]<-dS
mat[1,7]<-alfa #to estimate the fraction (α) of substitutions driven to fixation by positive selection at the functional sites
mat[1,8]<-NI
mat[1,9]<-dos #Stoletzki, 2011
mat[1,10]<-pval
mat[1,11]<-"NA"
write.table(mat,append = T,file = paste("results/",contr,"_mafdP_",mdafP,"_mafdD_",mdafD,"_MKT.txt",sep=""),quote = F, sep = "\t",col.names = T,row.names = F)




###########HighFst AASs###########
  ###  -- THIS IS INCORPORATED INTO SCANTOOLS -- ###
  setwd("~/JICAutumn2016/finalAnalysis29Apr/")
  library(data.table)
  con<-c("TETDIP","PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA") #,"CROPAN","CRODIN","CROSEC"  ,"CROBAL","CROWCA"
  #con<-c("TETDIP")
  msum<-matrix(nrow = sum(table(con)), ncol = 18,dimnames =list(c(con),c("Smin","S1Q","Smedian","Smean","S3Q","Smax","Ssd", "Nmin","N1Q","Nmedian","Nmean","N3Q","Nmax","Nsd","Squantile95","Squantile99","Nquantile95","Nquantile99")))
  al <- fread("lists/ALmeiotic.ALcode.tab.Name.txt",header = F) ##for all genes: ALcodesAll.txt
  mat <- matrix(nrow = nrow(al), ncol = sum(table(con)),dimnames = list(c(al$V1),c(con))) 
  for (contr in  con){ #   contr="TETDIP"
    indlin<-which(con %in% paste(contr))
    inp<-fread(paste("data/",contr,"_WS1_MS1_BPM.txt",sep=""),header=T)
    inp<-as.data.frame(na.omit(inp))
    s<-subset(inp,inp$ann %in% "synonymous_variant")
    n<-subset(inp,inp$ann %in% "missense_variant")
    sq <- quantile(s$FstH,probs=c(.99))
    msum[as.numeric(indlin),1:6]<-summary (s$FstH)
    msum[as.numeric(indlin),7]<-sd(s$FstH)
    msum[as.numeric(indlin),8:13]<-summary (n$FstH)
    msum[as.numeric(indlin),14]<-sd(n$FstH)
    msum[as.numeric(indlin),15] <-quantile(s$FstH,probs=c(.95))
    msum[as.numeric(indlin),16]<-sq
    msum[as.numeric(indlin),17] <-quantile(n$FstH,probs=c(.95))
    msum[as.numeric(indlin),18]<-quantile(n$FstH,probs=c(.99))
    sq <- quantile(n$FstH,probs=c(.99))
    
#  } ###comment this
 # write.table(msum,append = F,file ="results/FstStatsPerContrast_10Aug2020.txt",quote = F, sep = "\t",col.names = T,row.names = T) ###comment this!!
  
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
    for (g in readLines("lists/ALmeiotic.ALcode.tab.Name.txt")) ##for all genes: ALcodesAll.txt
    {# g='AL8G38150'
      g<-substr(g,1,9)
      index<-which(al$V1 %in% paste(g))
      gene<-subset(x = inp,subset = inp$ALcode %in% g & inp$ann %in% "missense_variant")
      
      if (nrow(gene)==0) {high<-0} else 
      {high<-nrow(subset(gene,gene$FstH >= as.numeric(sq)))}
      mat[index,indlin]<-high
    }}
  write.table(mat,append = F,file = "results/FstHighperGene.AllLin.AllMeiosis.txt",quote = F, sep = "\t",col.names = T,row.names   = T)
  write.table(msum,append = F,file ="results/FstStatsPerContrast.txt",quote = F, sep = "\t",col.names = T,row.names = T)
  
  msum<-as.data.frame(msum)
  wilcox.test(msum$Squantile99,msum$Nquantile99)
  ###  -- THE END OF SCANTOOLS PART  -- ###
###COPY METACENTRUM OUTPUT TO /DATA!!

#####Histogram of number of high-Fst sites per gene genome-wide - for any contrast
library(data.table)
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/")
system("bash gatherFstperGene.sh")
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/results/posSelection/")
quant = 0.99 #quant = 0.999
hfall<-read.table(paste("../FstHighperGene.Allgenes.",quant,".txt", sep=""),header=T,row.names = 1)
colnames(hfall)[1]<-"TETDIP"
fst<-read.table(paste("../FstStatsPerContrast.",quant,".txt", sep=""),header=T)
colnames(fst)[15]<-"Squantile"
genes<-read.table("../../lists/ALmeiotic.ALcode.tab.Name.txt",header=F)
con<-c("TETDIP")
for (contr in  con){ # contr="TETDIP"
  indCon<-which(con %in% paste(contr))
  png(filename=paste("histGW.meiotic.highFst.",contr,quant,".png",sep=""), 
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
    {}}}
  dev.off()
  
########### Fst, FAAD_G, FAAD_S ###########
  setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/")
  system("bash gatherFstperGene.sh")
  con<-c("TETDIP","PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA")
  for (contr in  con){
    system(command = paste("head -n 1 allBPM/",contr,"_WS1_MS1_BPM.txt > allBPM/",contr,".synNon_WS1_MS1_BPM.txt",sep=""))
    system(command = paste("grep -E 'missense_variant|synonymous_variant' allBPM/",contr,"_WS1_MS1_BPM.txt >> allBPM/",contr,".synNon_WS1_MS1_BPM.txt",sep=""))}
######TETRAPLOIDS#######
#1. GetFstOutliers Fst
library(data.table)
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/results/posSelection/")
quant = 0.99 #quant = 0.999
fst<-read.table(paste("../FstStatsPerContrast.",quant,".txt", sep=""),header=T)
colnames(fst)[15]<-"Squantile"
con<-c("TETDIP")
if (file.exists(paste("tetraploid/fst/FstoutliersGW_",quant,".txt",sep=""))) file.remove(paste("tetraploid/fst/FstoutliersGW_",quant,".txt",sep=""))
for (contr in  con){ # contr="TETDIP"
  indCon<-which(con %in% paste(contr))
  inp<-fread(paste("head -n -1 ../../allBPM/",contr,".synNon_WS1_MS1_BPM.txt",sep=""),header=T)
  q<-fst$Squantile[indCon]
  n<-subset(inp,inp$ann %in% "missense_variant" & inp$FstH >= q)
  write.table(n, file = paste("tetraploid/fst/FstoutliersGW_",quant,".txt",sep=""), quote = F,append = T,row.names = F,col.names = T)}
#2. Overlap Fst, FAAD_S, FAAD_G
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/results/posSelection/")
library(VennDiagram)
library(dplyr)
quant = 0.99 #quant = 0.999
f1<-read.table(paste("tetraploid/fst/FstoutliersGW_",quant,".txt",sep=""),h=T)
#   d<-subset(f1,f1$ALcode %in% genes$V1)
for (lin in c("all_dip_tet")){ # lin = "all_dip_tet"
if (file.exists(paste("tetraploid/",lin,"_",quant,"_overlapsSGF.intervals",sep=""))) file.remove(paste("tetraploid/",lin,"_",quant,"_overlapsSGF.intervals",sep=""))
if (file.exists(paste("tetraploid/",lin,"_",quant,"_overlapsSGF.intervals",sep=""))) file.remove(paste("tetraploid/",lin,"_",quant,"_overlapsSGF.txt",sep=""))
  g1<-read.table(paste("tetraploid/faadSG/",lin,"_2.scores.0.01.grantham",sep=""),h=F)#..., dip AC, tet AC, dip AF, tet AF, V8:grantham, diploid faad, tetraploid faad, delta faad ###Chrom loci AA AF/2x AF/4X DAP_DAF1 DAP_DAF2 grant_matrix_score grant1 grant2 difference in scores
  s1<-read.table(paste("tetraploid/faadSG/",lin,"_2.scores.0.01.sift",sep=""),h=F)
  g<-as.character(interaction(g1$V1, g1$V2))
  s<-as.character(interaction(s1$V1, s1$V2))
  f<-as.character(interaction(f1$scaff, f1$end))
  length(Reduce(intersect, list(g,s,f)))
  length(Reduce(intersect, list(g,s)))
  length(Reduce(intersect, list(g,f)))
  length(Reduce(intersect, list(f,s)))
  aa<-Reduce(intersect, list(g,s,f)) #,s
  aa1<-gsub('.', ':', aa,fixed = T)
  v<-venn.diagram(x=list("Fst_H"=f,"FAAD_S"=s,"FAAD_G"=g),paste("tetraploid/",lin,"_venn99.tiff",sep=""), lty = "blank", fill = c("blue1","red2","green"),alpha=0.3 , main = lin,cex=1, cat.cex=1,sub = "1% outliers",height = 2000,width = 2500)
  write.table(aa1,paste("tetraploid/",lin,"_",quant,"_overlapsSGF.intervals",sep=""),quote = F,col.names = F,row.names = F,append = T)
  d<-read.table(paste("tetraploid/",lin,"_",quant,"_overlapsSGF.intervals",sep=""),h=F,sep=":")
  fa<-subset(f1, as.character(interaction(f1$scaff, f1$end)) %in% as.character(interaction(d$V1, d$V2)))
  ga1<-subset(g1, as.character(interaction(g1$V1, g1$V2)) %in% as.character(interaction(d$V1, d$V2)))
  sa1<-subset(s1, as.character(interaction(s1$V1, s1$V2)) %in% as.character(interaction(d$V1, d$V2)))
  sa<-sa1 %>%  distinct(V1, V2, .keep_all = TRUE) 
  ga<-ga1 %>%  distinct(V1, V2, .keep_all = TRUE)
  all<-cbind(fa[,],sa$V8,ga[,c(8,11)]) #"NA"
  colnames(all)[21:23]<-c("1_SIFT","grantham","deltaFAAD")
  write.table(all,paste("tetraploid/",lin,"_",quant,"_overlapsSGF.txt",sep=""),quote = F,col.names = T,row.names = F,append = T)}
# 3. extract meiosis
library(stringr)
library(readr)
genes<-read.table("../../lists/ALmeiotic.ALcode.tab.Name.txt",header=F) ### .noSIFT.txt
a<-read.table("../../lists/pI_hyd.txt",header = T)
n<-read.table("../../lists/ALmeiotic.ALcode.tab.Name.txt",h=F)
i<-read.table("../../pai_Malvidae/Alignment_Identities_allSites.txt",h=F)
all<-read.table(paste("tetraploid/",lin,"_",quant,"_overlapsSGF.txt",sep=""), h=T)
d<-subset(all,all$ALcode %in% genes$V1)
d$anc<-substr(x = d$AAS,start = 1,stop = 3)
d$der<-str_sub(string = d$AAS, start= -3)
indAnc<-which(colnames(d) %in% "anc")
indDer<-which(colnames(d) %in% "der")
d$anc_pI<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,4])
d$der_pI<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,4])
d$delta_pI<-abs(d$anc_pI-d$der_pI)
d$anc_hydCA<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,2])
d$der_hydCA<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,2])
d$delta_hydCA<-abs(d$anc_hydCA-d$der_hydCA)
d$anc_hydCC<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,3])
d$der_hydCC<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,3])
d$delta_hydCC<-abs(d$anc_hydCC-d$der_hydCC)
indAl<-which(colnames(d) %in% "ALcode")
d$prot<-apply(X = d,1, function(x) subset(x = n, subset = n$V1 %in% x[indAl])[1,2])
d$aas_pos<-as.character(parse_number(as.character(d$AAS)))
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
d1<-d[ order(d[,2], d[,11]), ]
write.table(d1, paste("tetraploid/",lin,"_",quant,"_overlapsSGF.meiosis.txt",sep=""),append = F,quote = F,sep = "\t",row.names = F)

######DIPLOIDS#######
#1. Get Fst outliers genome wide
library(data.table)
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/")
system("bash gatherFstperGene.sh")
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/results/posSelection/")
quant = 0.99 #quant = 0.999
fst<-read.table(paste("../FstStatsPerContrast.",quant,".txt", sep=""),header=T)
colnames(fst)[15]<-"Squantile"
con<-c("PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA")
if (file.exists(paste("diploid/fst/FstoutliersGW_",quant,".txt",sep=""))) file.remove(paste("diploid/fst/FstoutliersGW_",quant,".txt",sep=""))
for (contr in  con){ # contr="PANDIN"
  indCon<-which(con %in% paste(contr))
  inp<-fread(paste("head -n -1 ../../allBPM/",contr,".synNon_WS1_MS1_BPM.txt",sep=""),header=T)
  q<-fst$Squantile[indCon+1]
  n<-subset(inp,inp$ann %in% "missense_variant" & inp$FstH >= q)
  n$FstHrel<-n$FstH/q
  write.table(n, file = paste("diploid/fst/FstoutliersGW_",quant,".txt",sep=""), quote = F,append = T,row.names = F,col.names = T)}
#2. Overlap Fst, FAAD_S, FAAD_G
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/results/posSelection/")
library(VennDiagram)
library(dplyr)
quant = 0.99
f1<-read.table(paste("diploid/fst/FstoutliersGW_",quant,".txt",sep=""),h=T)
f1<-subset(f1, f1$outname != "outname")
l<-c("Baltic","Wcarp","Secarp","Dinaric","Pannonian")
linm<-c("BAL","WCA","SEC","DIN","PAN")
for (lin in l){ # lin = "Wcarp"
  iLin<-which(l %in% paste(lin))
  g1<-read.table(paste("diploid/faadSG/",lin,"_0.01_2.outliers.sorted",sep=""),h=F) #DAF_DAP per lineage,sift,FAAD
  s1<-read.table(paste("diploid/faadSG/",lin,"_0.01_2.outliers.SIFT.SCCetc.sorted",sep=""),h=F)
  g<-as.character(interaction(g1$V1, g1$V2))
  s<-as.character(interaction(s1$V1, s1$V2))
  f<-as.character(interaction(f1$scaff, f1$end))
  length(Reduce(intersect, list(g,s,f)))
  length(Reduce(intersect, list(g,s)))
  length(Reduce(intersect, list(g,f)))
  length(Reduce(intersect, list(f,s)))
  aa<-Reduce(intersect, list(g,s,f)) ###  
  aa1<-gsub('.', ':', aa,fixed = T)
  v<-venn.diagram(x=list("Fst_H"=f,"FAAD_S"=s,"FAAD_G"=g),paste("diploid/",lin,"_venn99.tiff",sep=""), lty = "blank", fill = c("blue1","red2","green"),alpha=0.3 , main = lin,cex=1, cat.cex=1,sub = "1% outliers",height = 2000,width = 2500)
  write.table(aa1,paste("diploid/",lin,"_",quant,"_overlapsSGF.intervals",sep=""),quote = F,col.names = F,row.names = F,append = F)
  d<-read.table(paste("diploid/",lin,"_",quant,"_overlapsSGF.intervals",sep=""),h=F,sep=":")
  d<-d[ order(d[,1], d[,2]), ]
  fa<-subset(f1, as.character(interaction(f1$scaff, f1$end)) %in% as.character(interaction(d$V1, d$V2)))
  fa1<-subset(fa, fa$outname %like% linm[iLin])
  fa1<-fa1[ order(fa1[,2], fa1[,11],as.numeric(fa1[,21]),decreasing = T), ]
  fa2<-fa1[!duplicated(fa1[,c('scaff','end')]),] 
  fa<-fa2[ order(fa2[,2], as.integer(as.character(fa2[,11]))), ]
  d<-fa[,c(2,11)]
  ga1<-subset(g1, as.character(interaction(g1$V1, g1$V2)) %in% as.character(interaction(d$scaff, d$end)))
  sa1<-subset(s1, as.character(interaction(s1$V1, s1$V2)) %in% as.character(interaction(d$scaff, d$end)))
  sa<-sa1 %>%  distinct(V1, V2, .keep_all = TRUE) #might not be needed
  ga<-ga1 %>%  distinct(V1, V2, .keep_all = TRUE)
  all<-cbind(fa[,],sa[,c(9,9+iLin)],ga[,c(9)]) #   "NA","NA"
  colnames(all)[22:24]<-c("1_SIFT","FAAD","grantham")
  write.table(all,paste("diploid/",lin,"_",quant,"_overlapsSGF.txt",sep=""),quote = F,col.names = T,row.names = F,append = F)} ##
# 3. extract meiosis
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/results/posSelection/")
library(stringr)
library(readr)
quant = 0.99
genes<-read.table("../../lists/ALmeiotic.ALcode.tab.Name.txt",header=F) ## noSIFT.
a<-read.table("../../lists/pI_hyd.txt",header = T)
n<-read.table("../../lists/ALmeiotic.ALcode.tab.Name.txt",h=F) ## .noSIFT
i<-read.table("../../pai_Malvidae/Alignment_Identities_allSites.txt",h=F)
l<-c("Baltic","Secarp","Dinaric","Pannonian") ##"Wcarp",
if (file.exists(paste("diploid/all_diploid_",quant,"_overlapsSGF.meiosis.txt",sep=""))) file.remove(paste("diploid/all_diploid_",quant,"_overlapsSGF.meiosis.txt",sep=""))
for (lin in l){ # lin = "Wcarp"
  iLin<-which(l %in% paste(lin))
all<-read.table(paste("diploid/",lin,"_",quant,"_overlapsSGF.txt",sep=""), h=T) ##
d<-subset(all,all$ALcode %in% genes$V1)
d$anc<-substr(x = d$AAS,start = 1,stop = 3)
d$der<-str_sub(string = d$AAS, start= -3)
indAnc<-which(colnames(d) %in% "anc")
indDer<-which(colnames(d) %in% "der")
d$anc_pI<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,4])
d$der_pI<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,4])
d$delta_pI<-abs(d$anc_pI-d$der_pI)
d$anc_hydCA<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,2])
d$der_hydCA<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,2])
d$delta_hydCA<-abs(d$anc_hydCA-d$der_hydCA)
d$anc_hydCC<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indAnc])[1,3])
d$der_hydCC<-apply(X = d,1, function(x) subset(x = a, subset = a$AA %in% x[indDer])[1,3])
d$delta_hydCC<-abs(d$anc_hydCC-d$der_hydCC)
indAl<-which(colnames(d) %in% "ALcode")
d$prot<-apply(X = d,1, function(x) subset(x = n, subset = n$V1 %in% x[indAl])[1,2])
d$aas_pos<-as.character(parse_number(as.character(d$AAS)))
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
d1<-d[ order(d[,2], d[,11]), ]
d1$lin<-lin
write.table(d1, paste("diploid/",lin,"_",quant,"_overlapsSGF.meiosis.txt",sep=""),append = F,quote = F,sep = "\t",row.names = F)
write.table(d1, paste("diploid/all_diploid_",quant,"_overlapsSGF.meiosis.txt",sep=""),append = T,quote = F,sep = "\t",row.names = F)} ##
dd<-read.table(paste("diploid/all_diploid_",quant,"_overlapsSGF.meiosis.txt",sep=""), h=T)
dd1<-subset(dd, dd$outname != "outname")
write.table(dd1, paste("diploid/all_diploid_",quant,"_overlapsSGF.meiosis.txt",sep=""),append = F,quote = F,sep = "\t",row.names = F)
write.table(cbind(as.character(dd1$scaff),":",as.character(dd1$end)), paste("diploid/all_diploid_",quant,"_overlapsSGF.meiosis.intervals",sep=""),append = F,quote = F,sep = "",row.names = F,col.names = F)

#############High effect changes
library(data.table)
setwd(dir = "~/JICAutumn2016/finalAnalysis29Apr/results/posSelection/")
genes<-read.table("../../lists/ALmeiotic.ALcode.tab.Name.txt",header=F)
quant = 0.99 #quant = 0.999
fst<-read.table(paste("../FstStatsPerContrast.",quant,".txt", sep=""),header=T)
colnames(fst)[15]<-"Squantile"
con<-c("TETDIP")
if (file.exists(paste("tetraploid/fst/FstoutliersGW_",quant,".txt",sep=""))) file.remove(paste("tetraploid/fst/FstoutliersGW_",quant,".txt",sep=""))
for (contr in  con){ # contr="TETDIP"
  indCon<-which(con %in% paste(contr))
  inp<-fread(paste("head -n -1 ../../allBPM/",contr,"_WS1_MS1_BPM.txt",sep=""),header=T)
  q<-fst$Squantile[indCon]
  n<-subset(inp,inp$ann %in% "splice_acceptor_variant" | inp$ann %in% "splice_donor_variant" | inp$ann %in% "splice_region_variant" | inp$ann %in% "stop_lost" | inp$ann %in% "start_lost" | inp$ann %in% "stop_gained")
  k<- subset(n,n$FstH >= q)
  d<-subset(k,k$ALcode %in% genes$V1)
  d$prot<-apply(X = d,1, function(x) subset(x = genes, subset = genes$V1 %in% x[7])[1,2])
  write.table(d, file = paste("tetraploid/fst/highEffect.txt",sep=""), quote = F,append = F,row.names = F,col.names = T)}

####Functionality
#alIdent
library(ggpubr)

setwd("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/")
d<-read.table("posSelection/FinalListDiploids7Aug.txt",h=T)
pd<-read.table("posSelection/FinalListDiploidsPan7Aug.txt",h=T)
t<-read.table("posSelection/FinalListTetraploids7Aug.txt", h=T)
a<-read.table("file:///home/aa/JICAutumn2016/finalAnalysis29Apr/pai_Malvidae/Alignment_Identities_allSites.txt",h=F)
a<-subset(a,a$V2<10)
novA <- cbind("A",a$V3)
novT <- cbind("T",t$al_ident)
novPD <- cbind(as.character(pd$contract),pd$al_ident)
novD <- cbind("D",d$al_ident)
all<-as.data.frame(rbind(novT,novD))
allPD<-as.data.frame(rbind(novT,novPD)) #,novA

#Only diploid-tetraploid
allPD$V2<-as.numeric(as.character(allPD$V2))
my_comparisons <- list( c("D", "T"), c("D","P"),c("P","T") )
pdf("alIdent.pdf",width = 4,height = 4,pointsize = 24)
ggviolin(allPD, x = "V1", y = "V2", fill = "V1",
         palette = c("red","orange", "blue"),
         add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = 1.5)
dev.off()


#Only diploid-tetraploid
all$V2<-as.numeric(as.character(all$V2))
my_comparisons <- list( c("D", "T") )
    ggviolin(all, x = "V1", y = "V2", fill = "V1",
             palette = c("red", "blue"),
             add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = 1.5)
    
    
    
    

png(filename=paste("functionality/al_identPloidy.png", sep=""), 
    units = "px",
    width=875, 
    height=875, 
    pointsize=16*2)
boxplot(as.numeric(as.character(all$V2))~all$V1, col=c("red","blue"),ylab="conservatism (PAI)", names = c("In diploids", "In tetraploids"))
title ("Alignment identity")
dev.off()
png(filename = paste("functionality/al_identPloidy.hist.png"),width = 1300,height = 850,pointsize = 24)
hist(t$al_ident, col=rgb(0,0,1,0.6),xlim=c(0,1), xlab="Identity",breaks = 100,main = "Conservatism (PAI)",freq = F)
hist(d$al_ident, col=rgb(1,0,0,0.6), add=T,breaks = 100,freq = F)
legend(legend = c("tetraploids", "diploids"),fill = c("blue", "red"),x = "topleft")
dev.off()
summary(all)
t.test(t$al_ident, d$al_ident)
t.test(t$al_hydr, d$al_hydr)
t.test(t$al_pI, d$al_pI)
wilcox.test(t$al_ident, d$al_ident)
wilcox.test(t$al_hydr, d$al_hydr)
wilcox.test(t$al_pI, d$al_pI)
#HydCC
novT <- cbind("T",t$delta_hydCC)
novD <- cbind("D",d$delta_hydCC)
all<-as.data.frame(rbind(novT,novD))
png(filename=paste("functionality/delta_hydCCPloidy.png", sep=""), 
    units = "px",
    width=875, 
    height=875, 
    pointsize=16*2)
boxplot(as.numeric(as.character(all$V2))~all$V1, col=c("red","blue"),ylab="Hydrophobicity CC", names = c("In diploids", "In tetraploids"))
title ("Hydrophobicity")
dev.off()
png(filename = paste("functionality/delta_hydCCPloidy.hist.png"),width = 1300,height = 850,pointsize = 24)
hist(t$delta_hydCC, col=rgb(0,0,1,0.6), xlab="Hydrophobicity",breaks = 100,main = "Hydrophobicity CC",freq = F)
hist(d$delta_hydCC, col=rgb(1,0,0,0.6), add=T,breaks = 100,freq = F)
legend(legend = c("tetraploids", "diploids"),fill = c("blue", "red"),x = "topright")
dev.off()
wilcox.test(t$delta_hydCC, d$delta_hydCC)
#HydCA
novT <- cbind("T",t$delta_hydCA)
novD <- cbind("D",d$delta_hydCA)
all<-as.data.frame(rbind(novT,novD))
png(filename=paste("functionality/delta_hydCAPloidy.png", sep=""), 
    units = "px",
    width=875, 
    height=875, 
    pointsize=16*2)
boxplot(as.numeric(as.character(all$V2))~all$V1, col=c("red","blue"),ylab="Hydrophobicity CA", names = c("In diploids", "In tetraploids"))
title ("Hydrophobicity")
dev.off()
png(filename = paste("functionality/delta_hydCAPloidy.hist.png"),width = 1300,height = 850,pointsize = 24)
hist(t$delta_hydCA, col=rgb(0,0,1,0.6), xlab="Hydrophobicity",breaks = 100,main = "Hydrophobicity CA",freq = F)
hist(d$delta_hydCA, col=rgb(1,0,0,0.6), add=T,breaks = 100,freq = F)
legend(legend = c("tetraploids", "diploids"),fill = c("blue", "red"),x = "topright")
dev.off()
t.test(t$delta_hydCA, d$delta_hydCA)
wilcox.test(t$delta_hydCA, d$delta_hydCA)
#pI
novT <- cbind("T",t$delta_pI)
novD <- cbind("D",d$delta_pI)
all<-as.data.frame(rbind(novT,novD))
png(filename=paste("functionality/delta_pIPloidy.png", sep=""), 
    units = "px",
    width=875, 
    height=875, 
    pointsize=16*2)
boxplot(as.numeric(as.character(all$V2))~all$V1, col=c("red","blue"),ylab="pI", names = c("In diploids", "In tetraploids"))
title ("pI")
dev.off()
png(filename = paste("functionality/delta_pIPloidy.hist.png"),width = 1300,height = 850,pointsize = 24)
hist(t$delta_pI, col=rgb(0,0,1,0.6), xlab="pI",breaks = 100,main = "pI",freq = F)
hist(d$delta_pI, col=rgb(1,0,0,0.6), add=T,breaks = 100,freq = F)
legend(legend = c("tetraploids", "diploids"),fill = c("blue", "red"),x = "topright")
dev.off()
t.test(t$delta_pI, d$delta_pI)
wilcox.test(t$delta_pI, d$delta_pI)

###Make barplots for tree
setwd("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/")
d<-read.table("posSelection/FinalListDiploids7Aug.txt",h=T)
t<-read.table("posSelection/FinalListTetraploids7Aug.txt", h=T)
a<-data.frame(t(table(t$prot)))
aa<-subset(a, a$Freq != 1)
df <- aa[order(aa$Freq,decreasing = TRUE),]
xx<-nrow(df)
yy<-max(df$Freq)+3
png(filename = paste("tetraploids.png"),width = 850,height = 850,pointsize = 32,bg = NA)
dd<-barplot(df$Freq,names.arg = df$Var2,las=2,col="blue3",ylim = c(0,yy),xlim = c(0.7,13.3),axes = F,mgp = c(3, 0.3, 0))
text(x = dd,labels = df$Freq,y = df$Freq+1.65,cex = 1.15)
dev.off()

pan<-subset(d,d$FAAD_SFAAD_Goutlier %in% "Pannonian")
a<-data.frame(t(table(pan$prot)))
aa<-subset(a, a$Freq != 1 & a$Freq != 0)
df <- aa[order(aa$Freq,decreasing = TRUE),]
png(filename = paste("pannonian.png"),width = 850,height = 850,pointsize = 32,bg = NA)
dd<-barplot(df$Freq,names.arg = df$Var2,las=2,col="orange",ylim = c(0,yy),xlim = c(0.7,13.3),axes = F,mgp = c(3, 0.3, 0))
text(x = dd,labels = df$Freq,y = df$Freq+1.65,cex = 1.15)
dev.off()

bal<-subset(d,d$FAAD_SFAAD_Goutlier %in% "Baltic")
a<-data.frame(t(table(bal$prot)))
aa<-subset(a, a$Freq != 1 & a$Freq != 0)
df <- aa[order(aa$Freq,decreasing = TRUE),]
png(filename = paste("baltic.png"),width = 850,height = 850,pointsize = 32,bg = NA)
dd<-barplot(df$Freq,names.arg = df$Var2,las=2,col="magenta4",ylim = c(0,yy),xlim = c(0.7,13.3),axes = F,mgp = c(3, 0.3, 0))
text(x = dd,labels = df$Freq,y = df$Freq+1.65,cex = 1.15)
dev.off()

sec<-subset(d,d$FAAD_SFAAD_Goutlier %in% "Secarp")
a<-data.frame(t(table(sec$prot)))
aa<-subset(a, a$Freq != 1 & a$Freq != 0)
df <- aa[order(aa$Freq,decreasing = TRUE),]
png(filename = paste("Secarp.png"),width = 850,height = 850,pointsize = 32,bg = NA)
dd<-barplot(df$Freq,names.arg = df$Var2,las=2,col="steelblue1",ylim = c(0,yy),xlim = c(0.7,13.3),axes = F,mgp = c(3, 0.3, 0))
text(x = dd,labels = df$Freq,y = df$Freq+1.65,cex = 1.15)
dev.off()


####GO enrichment analysis
#install.packages("hash")
#source("http://bioconductor.org/biocLite.R")
#biocLite("GenomeInfoDb")
#biocLite("goProfiles")
#biocLite("GenomicRanges")
#install.packages('devtools')
#library(devtools)
#install("/home/aa/Desktop/programs/SNP2GO/")
#install_url("https://www.rdocumentation.org/packages/GenomeInfoDb/")
setwd("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/SNP2GO/")
system("cut -f2,11 ../../data/TETDIP.synNon_WS1_MS1_BPM.txt | awk 'NR == 1 || NR % 1000 == 0' > nonsel.txt") #for final better NR % 500 == 0
library(SNP2GO)
###DIPLOIDS
l<-c("Baltic","Wcarp","Secarp","Dinaric","Pannonian")
for (lin in l){ # lin = "Wcarp"
  s<-read.table(paste("../posSelection/diploid/",lin,"_0.99_overlapsSGF.txt",sep=""),h=T) 
  snps<-as.data.frame(s$scaff)
  snps$end <- as.numeric(s$end)
  cand <- GRanges(seqnames=snps[,1],ranges=IRanges(snps[,2],snps[,2]))
  noncand <- read.table("nonsel.txt",h=F)
  noncand[,2] <- as.numeric(noncand[,2])
  noncand <- GRanges(seqnames=noncand[,1],ranges=IRanges(noncand[,2],noncand[,2]))
  y <- snp2go(gtf="/home/aa/Desktop/references/lyrataV2/LyV2.gtf", goFile="/home/aa/Desktop/references/lyrataV2/mart_Alyrata_V2", candidateSNPs=cand, noncandidateSNPs=noncand, FDR=0.05, runs=10000, extension=0, min.regions=1)
  write.table(file=paste(lin,"_snp2go.tsv",sep=""),y$enriched,sep="\t",row.names=F)}

###TETRAPLOIDS
l<-c("all_dip_tet")
for (lin in l){ # lin = "all_dip_tet"
  s<-read.table(paste("../posSelection/tetraploid/all_dip_tet_0.99_overlapsSGF.txt",sep=""),h=T) 
  snps<-as.data.frame(s$scaff)
  snps$end <- as.numeric(as.character(s$end))
  cand <- GRanges(seqnames=snps[,1],ranges=IRanges(snps[,2],snps[,2]))
  noncand <- read.table("nonsel.txt",h=F)
  noncand[,2] <- as.numeric(noncand[,2])
  noncand <- GRanges(seqnames=noncand[,1],ranges=IRanges(noncand[,2],noncand[,2]))
  y <- snp2go(gtf="/home/aa/Desktop/references/lyrataV2/LyV2.gtf", goFile="/home/aa/Desktop/references/lyrataV2/mart_Alyrata_V2", candidateSNPs=cand, noncandidateSNPs=noncand, FDR=0.05, runs=10000, extension=0, min.regions=1)
  write.table(file=paste(lin,"_snp2go.tsv",sep=""),y$enriched,sep="\t",row.names=F)}
##To explore more:
# Get all enriched GO terms of GFF analysis:
gff.significant.terms <- y$enriched$GO
# Get the first of the enriched GO terms:
first.term <- gff.significant.terms[1] 
# Print all regions associated with at least one GO term:
print(y$regions)
# Print the regions associated with the first of the enriched GO terms:
# There are two possibilities to do so:
# version 2:
#print(y$regions[unlist(as.list(y$go2ranges[["regions"]][first.term]))])
#genes associated with enriched GO terms:
z <- nrow(y$enriched)
sgo<-y$regions[unlist(as.list(y$go2ranges[["regions"]][gff.significant.terms[1:z]]))]
## Although version 2 seems more complicated, it allows to get the regions
## associated with more than one term. In the following example, all regions associated
## with the first ten enriched GO terms (gff.significant.terms[1:10]) are printed:
print(y$regions[unlist(as.list(y$go2ranges[["regions"]][gff.significant.terms[1:10]]))])
# Print the candidate SNPs associated with the first of the enriched GO terms:
# Like for the GO regions, there are also two possibilities to do that:
# version 2:
print(cand[unlist(as.list(y$go2ranges[["candidates"]][first.term]))])
## Print the noncandidate SNPs associated with the first of the enriched GO terms:
## version 2:
print(noncand[unlist(as.list(x$go2ranges[["noncandidates"]][first.term]))])
# Get number of informative candidates of the GTF analysis:
z <- y$informative.candidateSNPs
#y$enriched
# Print the list of all GO terms associated to at least one gene in the
write.table(file="snp2go_gtf.tsv",y$enriched,sep="\t",row.names=F)
#write.table(file = "SNP2GOassocSNPs.txt",x = cand[unlist(as.list(y$go2ranges[["candidates"]][1:271]))])

################### Rank
setwd("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/SNP2GO/")
#1. get list of GO terms - grep them from mart
a<-read.table("GOlistTetsLarge.txt",h=T,sep = "\t")
b<-read.table("../../../../Desktop/references/lyrataV2/mart_Alyrata_V2.txt")
c<-read.table("../posSelection/tetraploid/all_dip_tet_0.99_overlapsSGF.txt",h=T)
bb<-subset(b,b$V1 %in% c$ALcode)
b<-subset(bb,bb$V2 %in% a$term_ID)
bb<-b[order(b$V2,decreasing = F),]
bb$V3<-apply(X = bb,1, function(x) subset(x = a, subset = a$term_ID %in% x[2])[1,2])
#FstH
cc<-c[order(c$FstH,decreasing = T),]
cc$rank<-seq(1,nrow(cc),1)
ccc<-cc[!duplicated(cc$ALcode), ]
bb$FstH<-apply(X = bb,1, function(x) subset(x = ccc, subset = ccc$ALcode %in% x[1])[1,24])
#AFD
cc<-c[order(abs(c$AFD),decreasing = T),]
cc$rank<-seq(1,nrow(cc),1)
ccc<-cc[!duplicated(cc$ALcode), ]
bb$AFD<-apply(X = bb,1, function(x) subset(x = ccc, subset = ccc$ALcode %in% x[1])[1,24])
#SIFT
cc<-c[order(c$X1_SIFT,decreasing = T),]
cc$rank<-seq(1,nrow(cc),1)
ccc<-cc[!duplicated(cc$ALcode), ]
bb$X1_SIFT<-apply(X = bb,1, function(x) subset(x = ccc, subset = ccc$ALcode %in% x[1])[1,24])
#Grantham
cc<-c[order(c$grantham,decreasing = T),]
cc$rank<-seq(1,nrow(cc),1)
ccc<-cc[!duplicated(cc$ALcode), ]
bb$grantham<-apply(X = bb,1, function(x) subset(x = ccc, subset = ccc$ALcode %in% x[1])[1,24])
bb$sumAFD_Fst<-bb$FstH+bb$AFD
bb$sumAll<-bb$FstH+bb$AFD+bb$X1_SIFT+bb$grantham
bbb<-bb[order(bb$sumAll,decreasing = F),]
bbb<-bb[order(bb$sumAFD_Fst,decreasing = F),]
summary(bbb)
write.table(bbb,"rankPathwaysTetraploids.txt",quote=F,row.names = F,sep="\t")

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
### 2 haplotypes
#get AF per lineage, combine into one matrix, also number anc, der
library (dplyr)
library(stringr)
library(readr)
setwd(dir = "/home/aa/JICAutumn2016/finalAnalysis29Apr/results/haplotypes/")
pops<-c('HNE', 'KZL', 'SZI', 'BEL', 'BIH', 'FOJ', 'GOR', 'RZA', 'VID', 'MIE', 'PRE', 'SNO', 'TRD', 'VEL', 'DRA', 'LAC', 'TZI', 'BRD', 'WEK', 'CHO', 'RFT', 'SWA', 'BGS', 'GUL', 'HOC', 'KAS', 'KOS', 'SCH', 'SPI', 'TKO', 'TRE', 'TRT', 'ZAP', 'KOW', 'STE', 'TBG')
#diploid
mm1<- read.table("../posSelection/diploid/all_diploid_0.99_overlapsSGF.meiosis.txt",h=T)
gf<- read.table("../posSelection/diploid/all_diploid_0.99_overlapsGF.meiosis.txt",h=T)
mm<-rbind(mm1,gf)
aa<-mm[!duplicated(mm[,c('ALcode')]),]
ff<-read.table("../posSelection/diploid/fst/FstoutliersGW_0.99.txt",h=T)
ff<-ff[!ff$ALcode %in% names(which(table(ff$ALcode) == 1)), ]
rr<-subset(ff, ff$ALcode %in% aa$ALcode & ff$outname %in% aa$outname)
n<-read.table("ALmeiotic.ALcode.tab.Name.InclZyp.txt")
indAl<-which(colnames(rr) %in% "ALcode")
rr$prot<-apply(X = rr,1, function(x) subset(x = n, subset = n$V1 %in% x[indAl])[1,2])
#write.table(mm,"xxx",row.names = F)
#OR
rr<-read.table(file = "diploidHIF_FINAL.txt",h=T)
#tetraploid
mm1<- read.table("../posSelection/tetraploid/all_dip_tet_0.99_overlapsSGF.meiosis.txt",h=T)
gf<- read.table("../posSelection/tetraploid/all_dip_tet_0.99_overlapsGF.meiosis.txt",h=T)
mm<-rbind(mm1,gf)
aa<-mm[!duplicated(mm[,c('ALcode')]),][,7]
ff<-read.table("../posSelection/tetraploid/fst/FstoutliersGW_0.99.txt",h=T)
ff<-ff[!ff$ALcode %in% names(which(table(ff$ALcode) == 1)), ]
rr<-subset(ff, ff$ALcode %in% aa)
write.table(file = 'quantile0.99.intervals',x = cbind(as.character(rr$scaff),rr$end),sep=":",col.names = F,row.names = F,quote = F)
system("java -Xmx10g -jar ~/Desktop/programs/GATK/GenomeAnalysisTK.jar -T SelectVariants -R ~/Desktop/references/lyrataV2/alygenomes.fasta -V ~/JICAutumn2016/heatmapPipelineAll/heatmapPipeline9March/data/All.missense.vcf -o /home/aa/alpine/scanTools/ScanTools/haplotypesDiploids/analys.vcf -L quantile0.99.intervals")
system("java -Xmx10g -jar ~/Desktop/programs/GATK/GenomeAnalysisTK.jar -T SelectVariants -R /home/aa/Desktop/references/arenosa/SCC1-4.MSH4.SMC6a.SHOC1.renamed.fasta -V /home/aa/alpine/scanTools/ScanTools/arenosa_meiosis/SHOC1.ann.vcf.gz -o /home/aa/alpine/scanTools/ScanTools/haplotypes300SHOC/SHOC.analys.vcf -L SHOC.quantile0.99.intervals")
system("java -Xmx10g -jar ~/Desktop/programs/GATK/GenomeAnalysisTK.jar -T SelectVariants -R ~/Desktop/references/lyrataV2/alygenomes.fasta -V /home/aa/alpine/scanTools/ScanTools/zyp1a/zyp1a.ann.variant.vcf -o /home/aa/alpine/scanTools/ScanTools/haplotypes300ZYP/ZYP.analys.vcf -L quantile0.99.intervals")

#WHY IT DOES NOT IMPORT??
#setwd("~/alpine/scanTools/ScanTools/")
#library(reticulate)
#use_python("/usr/bin/python3.4")
#repl_python()
#import ScanTools
#test = ScanTools.scantools("/home/aa/alpine/scanTools/ScanTools")
#test.splitVCFsAnn("haplotypes300",pops=['HNE', 'KZL', 'SZI', 'BEL', 'BIH', 'FOJ', 'GOR', 'RZA', 'VID', 'MIE', 'PRE', 'SNO', 'TRD', 'VEL', 'DRA', 'LAC', 'TZI', 'BRD', 'WEK', 'CHO', 'RFT', 'SWA', 'BGS', 'GUL', 'HOC', 'KAS', 'KOS', 'SCH', 'SPI', 'TKO', 'TRE', 'TRT', 'ZAP', 'KOW', 'STE', 'TBG'], repolarization_key="/home/aa/alpine/scanTools/ScanTools/repolarized.lookupKey.perSpeciesThreshold.txt", min_dp=0,mffg=0.9,partition='medium', mem=2000,overwrite=True, keep_intermediates=True, print1=False)
##test.splitVCFsAnn("haplotypes300ZYP",pops=['HNE', 'KZL', 'SZI', 'BEL', 'BIH', 'FOJ', 'GOR', 'RZA', 'VID', 'MIE', 'PRE', 'SNO', 'TRD', 'VEL', 'DRA', 'LAC', 'TZI', 'BRD', 'WEK', 'CHO', 'RFT', 'SWA', 'BGS', 'GUL', 'HOC', 'KAS', 'KOS', 'SCH', 'SPI', 'TKO', 'TRE', 'TRT', 'ZAP', 'KOW', 'STE', 'TBG'], repolarization_key="/home/aa/alpine/scanTools/ScanTools/repolarized.lookupKey.perSpeciesThreshold.txt", min_dp=0,mffg=0.9,partition='medium', mem=2000,overwrite=True, keep_intermediates=True, print1=False)
## Change PopKey
##test.splitVCFsAnn(vcf_dir="haplotypes300SHOC",repolarization_key="-99",min_dp=0,mffg=0.9,ref_path="/home/aa/Desktop/references/arenosa/SCC1-4.MSH4.SMC6a.SHOC1.renamed.fasta", gatk_path="~/Desktop/programs/GATK/GenomeAnalysisTK.jar", mem=10000, time='0-04:00', numcores=1, print1=False, overwrite=True, partition="nbi-long", keep_intermediates=False, use_scratch=False, scratch_path="/vcf", partition2="nbi-medium", time2="0-12:00",ref_name="SCC1-4.MSH4.SMC6a.SHOC1.renamed.fasta", pops=['HNE', 'KZL', 'SZI', 'BEL', 'BIH', 'FOJ', 'GOR', 'RZA', 'VID', 'MIE', 'PRE', 'SNO', 'TRD', 'VEL', 'DRA', 'LAC', 'TZI', 'BRD', 'WEK', 'CHO', 'RFT', 'SWA', 'BGS', 'GUL', 'HOC', 'KAS', 'KOS', 'SCH', 'SPI', 'TKO', 'TRE', 'TRT', 'ZAP', 'KOW', 'STE', 'TBG'])
#exit
#------------------------------------
# python3
# import ScanTools
# test = ScanTools.scantools("/home/aa/alpine/scanTools/ScanTools")
# test.splitVCFsAnn("haplotypes300",pops=['HNE', 'KZL', 'SZI', 'BEL', 'BIH', 'FOJ', 'GOR', 'RZA', 'VID', 'MIE', 'PRE', 'SNO', 'TRD', 'VEL', 'DRA', 'LAC', 'TZI', 'BRD', 'WEK', 'CHO', 'RFT', 'SWA', 'BGS', 'GUL', 'HOC', 'KAS', 'KOS', 'SCH', 'SPI', 'TKO', 'TRE', 'TRT', 'ZAP', 'KOW', 'STE', 'TBG'], repolarization_key="/home/aa/alpine/scanTools/ScanTools/repolarized.lookupKey.perSpeciesThreshold.txt", min_dp=0,mffg=0.9,partition='medium', mem=2000,overwrite=True, keep_intermediates=False, print1=False)
# mv VCF_haplotypes300_DP0.M0.9/*.table.repol.txt /home/aa/JICAutumn2016/finalAnalysis29Apr/results/haplotypes/
##Choose one:
contr<-c("diploid","tetraploid","other")
colors<-c("red3", "blue2", "grey")
#1.neutral - diploid
c1<-c('HNE', 'KZL', 'SZI', 'BEL', 'BIH', 'FOJ', 'GOR', 'RZA', 'VID', 'MIE', 'PRE', 'SNO', 'TRD', 'VEL')
#2.selected
c2<-c('DRA', 'LAC', 'TZI', 'BRD', 'WEK', 'CHO', 'RFT', 'SWA', 'BGS', 'GUL', 'HOC', 'KAS', 'KOS', 'SCH', 'SPI', 'TKO', 'TRE', 'TRT', 'ZAP', 'KOW', 'STE', 'TBG')

contr<-c("non-SCarp","SCarp","other")
colors<-c("darkslateblue", "steelblue1", "grey")
c1<-c('HNE', 'KZL', 'SZI', 'BEL', 'BIH', 'FOJ', 'MIE', 'PRE', 'SNO', 'TRD', 'VEL','DRA', 'LAC', 'TZI', 'BRD', 'WEK', 'CHO', 'RFT', 'SWA', 'BGS', 'GUL', 'HOC', 'KAS', 'KOS', 'SCH', 'SPI', 'TKO', 'TRE', 'TRT', 'ZAP', 'KOW', 'STE', 'TBG')
#2.selected
c2<-c('GOR', 'RZA', 'VID')

contr<-c("non-baltic","baltic","other")
colors<-c("darkslateblue", "magenta4", "grey")
c1<-c('HNE', 'KZL', 'SZI', 'BEL', 'BIH', 'FOJ', 'GOR', 'RZA', 'VID', 'SNO', 'TRD', 'VEL','DRA', 'LAC', 'TZI', 'BRD', 'WEK', 'CHO', 'RFT', 'SWA', 'BGS', 'GUL', 'HOC', 'KAS', 'KOS', 'SCH', 'SPI', 'TKO', 'TRE', 'TRT', 'ZAP', 'KOW', 'STE', 'TBG')
#2.selected
c2<-c('MIE', 'PRE')

a<-read.table(paste("/home/aa/alpine/scanTools/ScanTools/VCF_haplotypes300_DP0.M0.9/",pops[1],".table.repol.txt",sep=""),h=F)
#z<-read.table(paste("/home/aa/alpine/scanTools/ScanTools/VCF_haplotypes300ZYP_DP0.M0.9/",pops[1],".table.repol.txt",sep=""),h=F)
#s<-read.table(paste("/home/aa/alpine/scanTools/ScanTools/VCF_haplotypes300SHOC_DP0.M0.9/",pops[1],".table.recode.txt",sep=""),h=F)
#nnn<-rbind(z,a,s)
nnn<-a
#baltic:
nnn<-subset(a,a$V7 %in% "AL8G40620")
#scarp:
nnn<-subset(a,a$V7 %in% c("AL6G50600","AL8G38150"))

mtot<-matrix(nrow = nrow(nnn), ncol = sum(table(pops))+4,dimnames =list(c(),c("scaff","pos","alcode","AAS",pops)))
mtot[,1]<-as.character(nnn[,3])
mtot[,2]<-as.character(nnn[,4])
mtot[,3]<-as.character(nnn[,7])
mtot[,4]<-as.character(nnn[,9])
for (pop in  pops){ # pop="KZL"
  ip<-which(pops %in% paste(pop))
  a1<-read.table(paste("/home/aa/alpine/scanTools/ScanTools/VCF_haplotypes300_DP0.M0.9/",pop,".table.repol.txt",sep=""),h=F)
  #z<-read.table(paste("/home/aa/alpine/scanTools/ScanTools/VCF_haplotypes300ZYP_DP0.M0.9/",pop,".table.repol.txt",sep=""),h=F)
  #s<-read.table(paste("/home/aa/alpine/scanTools/ScanTools/VCF_haplotypes300SHOC_DP0.M0.9/",pop,".table.recode.txt",sep=""),h=F)
  #a<-rbind(z,a1,s)
  #a<-subset(a1,a1$V7 %in% "AL8G40620") #rmi1
  a<-subset(a1,a1$V7 %in% c("AL6G50600","AL8G38150")) #smc6b, pms1
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
for (i in  1:nrow(k)){ # i=1
  if (k$c1[i]>k$c2[i])
  {k[i,5:(ncol(k)-10)] <- 1 - k[i,5:(ncol(k)-10)]
  a<-k$anc[i]
  k$anc[i]<-k$der[i]
  k$der[i]<-a
  } else {}}
write.table(k[1:(ncol(k)-2)],"AFsperLinRepol.txt",append = F,quote = F,row.names = F,col.names = T)
#calculate haplotype allele frequencies
k<-read.table("AFsperLinRepol.txt",h=T)
g<-read.table("PopGeogr.txt",h=T)
for (al in unique(k$alcode)){ # al<-"AL1G10680"
 a1<-subset(k,k$alcode %in% al)
 a2<-as.data.frame(t(a1[,5:(ncol(a1)-3)]))
 #neutral
 g$c1<-1-(apply(a2, 1, max))
 #selected
 g$c2<-apply(a2, 1, min)
 #recombined
 g$oth<-abs(1-g$c1-g$c2)
 g$nAA<-nrow(a1)
 g$length<-abs(max(a1$pos)-min(a1$pos))
 n<-read.table("ALmeiotic.ALcode.tab.Name.InclZyp.txt")
 name<-as.matrix(subset(n,n$V1 %in% al))[1,2]
 write.table(g,paste(al,"_",name,".hapl.txt",sep=""),quote = F,sep = "\t",row.names = F)}
#MAP
library(RgoogleMaps)
library(argosfilter)
library(mapplots)
list<-list.files(path = ".", pattern = "*.hapl.txt")
terrmap<-GetMap.bbox(lon,lat,maptype= "terrain", destfile = "terrain.png",path = paste0("&style=feature:road|element:all|visibility:off","&style=feature:administrative.country|element:geometry.stroke|visibility:off","&style=feature:all|element:labels|visibility:off"), MINIMUMSIZE = T)
size<-terrmap$size
#pdf(file=paste("figures/","all_scarp",".hapl.pdf", sep=""),width = size[1]/72, height = size[2]/72)
for (l in list){ # l="AL8G38150_SMC6B.hapl.txt"
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
  terrmap<-legend("bottomleft", inset=.005, legend = c(paste("Length = ",df$length[1]," bp",sep=""),paste("Markers = ",df$nAA[1]," AAs",sep="")), horiz=F, cex=1,bg = 'white')
  dev.off()
  }
#dev.off()

### 3 haplotypes
library (dplyr)
library(stringr)
library(readr)
library(data.table)
setwd(dir = "/home/aa/JICAutumn2016/finalAnalysis29Apr/results/haplotypes/3hapl/")
pops<-c('HNE', 'KZL', 'SZI', 'BEL', 'BIH', 'FOJ', 'GOR', 'RZA', 'VID', 'MIE', 'PRE', 'SNO', 'TRD', 'VEL', 'DRA', 'LAC', 'TZI', 'BRD', 'WEK', 'CHO', 'RFT', 'SWA', 'BGS', 'GUL', 'HOC', 'KAS', 'KOS', 'SCH', 'SPI', 'TKO', 'TRE', 'TRT', 'ZAP', 'KOW', 'STE', 'TBG')
contr<-c("baltic","diploid","tetraploid","other")
colors<-c("magenta4","red3", "blue2", "grey")
#1.neutral
c0<-c('HNE', 'KZL', 'SZI','BEL', 'BIH', 'FOJ', 'GOR', 'RZA', 'VID', 'MIE', 'PRE', 'SNO', 'TRD', 'VEL')
#2.selected
c1<-c('MIE', 'PRE')
c2<-c('DRA', 'LAC', 'TZI', 'BRD', 'WEK', 'CHO', 'RFT', 'SWA', 'BGS', 'GUL', 'HOC', 'KAS', 'KOS', 'SCH', 'SPI', 'TKO', 'TRE', 'TRT', 'ZAP', 'KOW', 'STE', 'TBG')
#diploid
mm<- read.table("../../posSelection/diploid/all_diploid_0.99_overlapsSGF.meiosis.txt",h=T)
mp<-subset(mm, mm$lin %in% "Pannonian")
aa<-mp[!duplicated(mp[,c('ALcode')]),][,7]
ff<-read.table("../../posSelection/diploid/fst/FstoutliersGW_0.99.txt",h=T)
ff<-ff[!ff$ALcode %in% names(which(table(ff$ALcode) == 1)), ]
rd<-subset(ff, ff$ALcode %in% aa & ff$outname %like% "PAN")
rd$end<-as.character(as.integer(as.character(rd$end)))
rd<-rd[ order(rd[,2], rd[,11]), ]
##OR:
a1<-read.table("../diploidHIF_FINAL.txt",h=T)
rd<-subset(a1,as.character(a1$ALcode) %in% c("AL4G46460","AL6G15380","AL6G30890","AL8G25590"))
#tetraploid
mm<- read.table("../../posSelection/tetraploid/all_dip_tet_0.99_overlapsSGF.meiosis.txt",h=T)
aa<-mm[!duplicated(mm[,c('ALcode')]),][,7]
ff<-read.table("../../posSelection/tetraploid/fst/FstoutliersGW_0.99.txt",h=T)
ff<-ff[!ff$ALcode %in% names(which(table(ff$ALcode) == 1)), ]
rt<-subset(ff, ff$ALcode %in% aa)
rt$end<-as.character(as.integer(as.character(rt$end)))
rt<-rt[ order(rt[,2], rt[,11]), ]
#to get interval list of HIA
rr<-rbind(rd[,c(2,11)],rt[,c(2,11)])
rr$lin<-"PAN"
rr$lin[(nrow(rt)+1):nrow(rr)]<-"TET"
rr<-rr[ order(rr[,1], rr[,2]), ]
rr<-rr[!duplicated(interaction(rr$scaff, rr$end,rr$lin)), ]
write.table(file = 'quantile0.99.intervals',x = cbind(as.character(rr$scaff),rr$end),sep=":",col.names = F,row.names = F,quote = F)
write.table(file = 'quantile0.99.lin.txt',x = cbind(as.character(rr$scaff),rr$end,rr$lin),sep="\t",col.names = F,row.names = F,quote = F)
system("java -Xmx10g -jar ~/Desktop/programs/GATK/GenomeAnalysisTK.jar -T SelectVariants -R ~/Desktop/references/lyrataV2/alygenomes.fasta -V ~/JICAutumn2016/heatmapPipelineAll/heatmapPipeline9March/data/All.missense.vcf -o ~/alpine/scanTools/ScanTools/haplotypes300/analys.vcf -L quantile0.99.intervals")
# python3
# import ScanTools
# test = ScanTools.scantools("/home/aa/alpine/scanTools/ScanTools")
# test.splitVCFsAnn("haplotypes300",pops=['HNE', 'KZL', 'SZI', 'BEL', 'BIH', 'FOJ', 'GOR', 'RZA', 'VID', 'MIE', 'PRE', 'SNO', 'TRD', 'VEL', 'DRA', 'LAC', 'TZI', 'BRD', 'WEK', 'CHO', 'RFT', 'SWA', 'BGS', 'GUL', 'HOC', 'KAS', 'KOS', 'SCH', 'SPI', 'TKO', 'TRE', 'TRT', 'ZAP', 'KOW', 'STE', 'TBG'], repolarization_key="/home/aa/alpine/scanTools/ScanTools/repolarized.lookupKey.perSpeciesThreshold.txt", min_dp=0,mffg=0.9,partition='medium', mem=10000,overwrite=True, keep_intermediates=False, print1=False)
# mv VCF_haplotypes300_DP0.M0.9/*.table.repol.txt /home/aa/JICAutumn2016/finalAnalysis29Apr/results/haplotypes/3hapl
nnn<-read.table(paste("/home/aa/alpine/scanTools/ScanTools/VCF_haplotypes300_DP0.M0.9/",pops[1],".table.repol.txt",sep=""),h=F)
#to get h1P,h2T and h3PT haplotype
rr<-read.table('quantile0.99.lin.txt')
#divide to unique and duplicated
hh<-rr[!duplicated(rr[,c('V1','V2')]),]
h1<-subset(hh,hh$V3 %in% "PAN")
h2<-subset(hh,hh$V3 %in% "TET")
h3<-rr[duplicated(rr[,c('V1','V2')]),]
write.table(h1[,1:2],"h1.txt",col.names = F,row.names = F,quote = F)
write.table(h2[,1:2],"h2.txt",col.names = F,row.names = F,quote = F)
write.table(h3[,1:2],"h3.txt",col.names = F,row.names = F,quote = F)
#get per pop allele freq for each haplotype separately
for (h in 1:3){ # h=1
  hx<-read.table(paste("h",h,".txt",sep=""),h=F)
  nnh<-subset(nnn, as.character(interaction(nnn$V3, nnn$V4)) %in% as.character(interaction(hx$V1, hx$V2))) 
  mtot<-matrix(nrow = nrow(nnh), ncol = sum(table(pops))+4,dimnames =list(c(),c("scaff","pos","alcode","AAS",pops)))
  mtot[,1]<-as.character(nnh[,3])
  mtot[,2]<-as.character(nnh[,4])
  mtot[,3]<-as.character(nnh[,7])
  mtot[,4]<-as.character(nnh[,9])
  for (pop in  pops){ # pop="KZL"
    ip<-which(pops %in% paste(pop))
    a1<-read.table(paste("/home/aa/alpine/scanTools/ScanTools/VCF_haplotypes300_DP0.M0.9/",pop,".table.repol.txt",sep=""), header =F)
    a<-subset(a1,as.character(interaction(a1$V3, a1$V4)) %in% as.character(interaction(hx$V1, hx$V2)))
    a[ a == "-9" ] <- NA
    a$na<-apply(is.na(a), 1, sum)
    a$an<-((ncol(a)-10)-a$na)*as.numeric(a[1,2])
    a$af<-rowSums(a[,10:(ncol(a)-2)],na.rm = T) / a$an
    mtot[,ip+4]<-a$af}
  write.table(mtot,paste("AFsperLin_h",h,".txt",sep=""),append = F,quote = F,row.names = F,col.names = T)}
#repolarize h1 - PAN ###########Odsud dal s jinymi populacemi v kontrastech!!!
k<-read.table(paste("AFsperLin_h1.txt",sep=""),h=T)
k$anc<-substr(x = k$AAS,start = 1,stop = 3)
k$der<-str_sub(string = k$AAS, start= -3)
k$aas_pos<-as.character(parse_number(as.character(k$AAS)))
k1<-select(.data = k,c1) #c1 - P
k2<-select(.data = k,c(c0,c2)) # T + D
k$c1<-rowSums(k1,na.rm = T)/ncol(k1)
k$c2<-rowSums(k2,na.rm = T)/ncol(k2)
for (i in  1:nrow(k)){ # i=1
  if (k$c1[i]<k$c2[i])
  {k[i,5:(ncol(k)-5)] <- 1 - k[i,5:(ncol(k)-5)]
  print(i)
  a<-k$anc[i]
  k$anc[i]<-k$der[i]
  k$der[i]<-a
  } else {}}
write.table(k[1:(ncol(k)-2)],"AFsperLinRepol_h1.txt",append = F,quote = F,row.names = F,col.names = T)
#repolarize h2 - TET
k<-read.table(paste("AFsperLin_h2.txt",sep=""),h=T)
k$anc<-substr(x = k$AAS,start = 1,stop = 3)
k$der<-str_sub(string = k$AAS, start= -3)
k$aas_pos<-as.character(parse_number(as.character(k$AAS)))
k1<-select(.data = k,c2) #c1 - T
k2<-select(.data = k,c(c0,c1)) # P + D
k$c1<-rowSums(k1,na.rm = T)/ncol(k1)
k$c2<-rowSums(k2,na.rm = T)/ncol(k2)
for (i in  1:nrow(k)){ # i=1
  if (k$c1[i]<k$c2[i])
  {k[i,5:(ncol(k)-5)] <- 1 - k[i,5:(ncol(k)-5)]
  print(i)
  a<-k$anc[i]
  k$anc[i]<-k$der[i]
  k$der[i]<-a
  } else {}}
write.table(k[1:(ncol(k)-2)],"AFsperLinRepol_h2.txt",append = F,quote = F,row.names = F,col.names = T)
#repolarize h3 - PAN + TET
k<-read.table(paste("AFsperLin_h3.txt",sep=""),h=T)
k$anc<-substr(x = k$AAS,start = 1,stop = 3)
k$der<-str_sub(string = k$AAS, start= -3)
k$aas_pos<-as.character(parse_number(as.character(k$AAS)))
k1<-select(.data = k,c(c1,c2)) #T+D
k2<-select(.data = k,c(c0)) # P + D
k$c1<-rowSums(k1,na.rm = T)/ncol(k1)
k$c2<-rowSums(k2,na.rm = T)/ncol(k2)
for (i in  1:nrow(k)){ # i=1
  if (k$c1[i]<k$c2[i])
  {k[i,5:(ncol(k)-5)] <- 1 - k[i,5:(ncol(k)-5)]
  print(i)
  a<-k$anc[i]
  k$anc[i]<-k$der[i]
  k$der[i]<-a
  } else {}}
write.table(k[1:(ncol(k)-2)],"AFsperLinRepol_h3.txt",append = F,quote = F,row.names = F,col.names = T)
#overall HAF calculation
g<-read.table("../PopGeogr.txt",h=T)
h1<-read.table("AFsperLinRepol_h1.txt",h=T)
h2<-read.table("AFsperLinRepol_h2.txt",h=T) 
h3<-read.table("AFsperLinRepol_h3.txt",h=T) 
for (al in unique(rbind(as.character(h1$alcode),as.character(h2$alcode),as.character(h3$alcode)))){ # al<-"AL6G15380"
  a1<-subset(h1,h1$alcode %in% al)
  h1p<-as.data.frame(t(a1[,5:(ncol(a1)-3)]))
  a2<-subset(h2,h2$alcode %in% al)
  h2t<-as.data.frame(t(a2[,5:(ncol(a2)-3)]))
  a3<-subset(h3,h3$alcode %in% al)
  h3pt<-as.data.frame(t(a3[,5:(ncol(a3)-3)]))
  g$h1p <-apply(h1p, 1, min, na.rm=T)
  g$h1dt<-1-(apply(h1p, 1, max, na.rm=T))
  g$h2t<-apply(h2t, 1, min, na.rm=T)
  g$h2dp<-1-(apply(h2t, 1, max, na.rm=T))
  g$h3pt<-apply(h3pt, 1, min, na.rm=T)
  g$h3d<-1-(apply(h3pt, 1, max, na.rm=T))
  #For equally long haplotypes:
  #g$haf1p<-apply((cbind(g$h1p, g$h2dp,g$h3pt)),1,min)
  #g$haf2t<-apply((cbind(g$h2t, g$h1dt,g$h3pt)),1,min)
  #g$haf0d<-apply((cbind(g$h3d, g$h1dt,g$h2dp)),1,min)
  g$haf1p<-apply((cbind(g$h1p, g$h3pt)),1,min)
  g$haf2t<-apply((cbind(g$h2t, g$h3pt)),1,min)
  g$haf0d<-apply((cbind(g$h3d, g$h1dt,g$h2dp)),1,min)
  g$oth<-1-g$haf1p-g$haf2t-g$haf0d
  g$nAAt<-nrow(a2)+nrow(a3)
  g$nAAp<-nrow(a1)+nrow(a3)
  g$nAAd<-nrow(a1)+nrow(a2)+nrow(a3)
  g$lengthT<-abs(max(rbind(a2$pos,a3$pos))-min(rbind(a2$pos,a3$pos)))
  g$lengthP<-abs(max(rbind(a1$pos,a3$pos))-min(rbind(a1$pos,a3$pos)))
  g$lengthD<-abs(max(rbind(a1$pos,a2$pos,a3$pos))-min(rbind(a1$pos,a2$pos,a3$pos)))
  n<-read.table("../ALmeiotic.ALcode.tab.Name.InclZyp.txt")
  name<-as.matrix(subset(n,n$V1 %in% al))[1,2]
  write.table(g,paste(al,"_",name,".hapl.txt",sep=""),quote = F,sep = "\t",row.names = F)}
#MAP
library(RgoogleMaps)
library(argosfilter)
library(mapplots)
list<-list.files(path = ".", pattern = "*.hapl.txt")
for (l in list){ # l="AL8G25590_DYAD.hapl.txt"
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
    terrmap<-add.pie(z=c(df$haf1p[i],df$haf0d[i],df$haf2t[i],df$oth[i]), x=a3, y=b3, radius=10, col = c("orange","red3", "blue2", "grey"), labels=NA)
    terrmap<-text( x=a3, y=b3, labels = df$pop[i],pos = df$pos[i],col = as.character(df$col[i]), offset = 0.83,cex = 0.7)}
  terrmap<-legend("topleft", inset=.005, title=paste(prot," haplotypes",sep=""),c("pannonian","diploid","tetraploid","other"),fill = c("orange","red3", "blue2", "grey"), horiz=F, cex=1,bg = 'white')
  terrmap<-legend("bottomleft", inset=.005, title = "Haplotype parameters",legend = c(paste("length pannonian = ",df$lengthP[1]," bp",sep=""),paste("length diploid = ",df$lengthD[1]," bp",sep=""),paste("length tetraploid = ",df$lengthT[1]," bp",sep=""),paste("markers pannonian = ",df$nAAp[1]," AAs",sep=""),paste("markers diploid = ",df$nAAd[1]," AAs",sep=""),paste("markers tetraploid = ",df$nAAt[1]," AAs",sep="")), horiz=F, cex=0.8,bg = 'white')
  dev.off()}

###BARPLOT
setwd("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/haplotypes/forBarplot/")
#ALL
pdf(file="AllHaplotypes.pdf",width = 11.69, height = 11.69,title = "All haplotypes",pointsize = 18) #8.27*2
par(mai=c(0,0.15,0.1,0)) #par(mai=c(0,0.15,0,0))
par(mfrow=c(14,1)) #plot.new()
prot<-as.data.frame(c('PRD3','ZYP1A','ZYP1B','ASY1','PDS5B','DYAD','SCC4','SHOC1'))
for(i in 1:nrow(prot)){ # i=1
  name=as.character(paste(as.character(prot[i,1]),'.hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[,8:10]
  rownames(d1) <- d0[,1]
  barplot(t(as.matrix(d1)),col = c("red3", "blue2", "grey"),border="NA", las=2,xlab = NA,ylab=NA, axes = F, axisnames = T,space=0)
  title(ylab=title, line=-1.8, cex.lab=1.2)}
prot<-as.data.frame(c('ASY3','SYN1','SMG7'))
for(i in 1:nrow(prot)){ # i=1
  name=as.character(paste(as.character(prot[i,1]),'.hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[,14:17]
  rownames(d1) <- d0[,1]
  barplot(t(as.matrix(d1)),col = c("orange", "blue2","red3", "grey"),border="NA", las=2,xlab = NA,ylab=NA, axes = F, axisnames = T,space=0)
  title(ylab=title, line=-1.8, cex.lab=1.2)}
prot<-as.data.frame(c('SMC6B'))
for(i in 1:nrow(prot)){ # i=1
  name=as.character(paste(as.character(prot[i,1]),'.hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[,8:10]
  rownames(d1) <- d0[,1]
  barplot(t(as.matrix(d1)),col = c("darkslateblue", "steelblue1", "grey"),border="NA", las=2,xlab = NA,ylab=NA, axes = F, axisnames = T,space=0)
  title(ylab=title, line=-1.8, cex.lab=1.2)}
prot<-as.data.frame(c('RMI1'))
for(i in 1:nrow(prot)){ # i=1
  name=as.character(paste(as.character(prot[i,1]),'.hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[,8:10]
  rownames(d1) <- d0[,1]
  barplot(t(as.matrix(d1)),col = c("darkslateblue", "magenta4", "grey"),border="NA", las=2,xlab = NA,ylab=NA, axes = F, axisnames = T,space=0)
  title(ylab=title, line=-1.8, cex.lab=1.2)}
dev.off()


###########COEVOLUTION 2020
library(psych)
corr.test(mtcars)$ci

###BARPLOT
setwd("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/haplotypes/forBarplot/")
#ALL
aa<-as.data.frame(matrix(nrow=22,ncol=11)) # 11
prot<-as.data.frame(c('PRD3','ZYP1A','ZYP1B','ASY1','PDS5B','DYAD','SCC4','SHOC1')) #'PRD3','ZYP1A','ZYP1B','ASY1','PDS5B','DYAD','SCC4','SHOC1'
for(i in 1:nrow(prot)){ # i=1
  name=as.character(paste(as.character(prot[i,1]),'.hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[c(15:36),8:10]  ##15:36 15:17,29:36
  rownames(aa) <- d0[c(15:36),1]
  aa[,i]<-d1[,2]
  colnames(aa)[i]<-title}
prot<-as.data.frame(c('ASY3','SYN1','SMG7')) #'ASY3','SYN1','SMG7'

for(i in 1:nrow(prot)){ # i=1
  name=as.character(paste(as.character(prot[i,1]),'.hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[c(15:36),14:17] ##15:36 15:17,29:36
  rownames(aa) <- d0[c(15:36),1]
  i=i+8 ### 8
  aa[,i]<-d1[,2]
  colnames(aa)[i]<-title}

#ALL
library(corrplot)

pdf("HAFcorrelations_corrplot_6prot.pdf")
M<-cor(aa)
corrplot(M, method="circle", order="hclust",title = "All 11 admixed populations")
#RUDERAL
aaa<-aa[8:11,-2]
M<-cor(aaa)
corrplot(M, method="circle", order="hclust",title = "4 ruderal populations")
#WCARP
aaa<-aa[4:7,]
M<-cor(aaa)
corrplot(M, method="circle", order="hclust",title = "4 Wcarp. populations")
#SCARP
aaa<-aa[1:3,]
M<-cor(aaa)
corrplot(M, method="circle", order="hclust",title = "3 Scarp. populations")
dev.off()

#ALLLLLLLL
aaa<-aa[3:14,]
m<-apply(aa,2,mean)
min(m)
max(m)
typeof(aaa)
a1<-as.numeric(as.matrix(aaa))
min(aaa)
max(aaa)
median(a1)
mean(a1)
summary(a1)
mean(as.numeric(droplevels(aaa)))
lapply(aaa, mean)
a2<-subset(a1,a1==1)

library(psych)
aaa<-corr.test(aa,adjust = "none",use = "pairwise")

pdf("HAFcorrelations_corrplot_6prot.pdf")

M<-cor(aa)
corrplot(aaa$r, method="circle", order="hclust",title = "All 22 tetraploid populations",p.mat = aaa$p,sig.level = 0.05,type = "lower",insig = "label_sig")
dev.off()

M<-cor(aa,method = "pearson")




library(RColorBrewer)
my_palette <- colorRampPalette(c("white","gold","red"))(n = 10)

a<-corr.test(aa)$r
write.table(a,"HAFcorrelations.txt",quote = F)
pdf("HAFcorrelations.pdf")
heatmap(a,symm = T,keep.dendro = F,col=my_palette)
heatmap(abs(a),symm = T,keep.dendro = F,col=my_palette)

dev.off()


###

########TETRAPLOID
setwd("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/haplotypes/forBarplot/")
aa<-as.data.frame(matrix(nrow=22,ncol=11))
prot<-as.data.frame(c('PRD3','ZYP1A','ZYP1B','ASY1','PDS5B','DYAD','SCC4','SHOC1'))
for(i in 1:nrow(prot)){ # i=1
  name=as.character(paste(as.character(prot[i,1]),'.hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[c(15:36),8:10]
  rownames(aa) <- d0[c(15:36),1]
  aa[,i]<-d1[,2]
  colnames(aa)[i]<-title}
prot<-as.data.frame(c('ASY3','SYN1','SMG7'))

for(i in 1:nrow(prot)){ # i=1
  name=as.character(paste(as.character(prot[i,1]),'.hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[c(15:36),14:17]
  rownames(aa) <- d0[c(15:36),1]
  i=i+8
  aa[,i]<-d1[,2]
  colnames(aa)[i]<-title}
#ALLLLLLLL
aaa<-aa[3:14,]
m<-apply(aaa,2,mean)
min(m)

#####DIPLOID
setwd("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/haplotypes/forBarplot/")
aa<-as.data.frame(matrix(nrow=14,ncol=11))
prot<-as.data.frame(c('PRD3','ZYP1A','ZYP1B','ASY1','PDS5B','DYAD','SCC4','SHOC1'))
for(i in 1:nrow(prot)){ # i=1
  name=as.character(paste(as.character(prot[i,1]),'.hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[c(1:14),8:10]
  rownames(aa) <- d0[c(1:14),1]
  aa[,i]<-d1[,1]
  colnames(aa)[i]<-title}
prot<-as.data.frame(c('ASY3','SYN1','SMG7'))

for(i in 1:nrow(prot)){ # i=1
  name=as.character(paste(as.character(prot[i,1]),'.hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[c(1:14),14:17]
#  rownames(aa) <- d0[c(15:36),3]
  i=i+8
  aa[,i]<-d1[,3]
  colnames(aa)[i]<-title}

aaa<-aa
m<-apply(aaa,2,mean)
min(m)
max(m)

#####PANNONIAN
setwd("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/haplotypes/forBarplot/")
aa<-as.data.frame(matrix(nrow=14,ncol=3))
prot<-as.data.frame(c('ASY3','SYN1','SMG7'))
for(i in 1:nrow(prot)){ # i=1
  name=as.character(paste(as.character(prot[i,1]),'.hapl.txt',sep=''))
  title=as.character(paste(as.character(prot[i,1])))
  d0<-read.table(name, header=T)
  d1<-d0[c(1:14),14:17]
  #  rownames(aa) <- d0[c(15:36),3]
 # i=i+8
  aa[,i]<-d1[,1]
  colnames(aa)[i]<-title}

aaa<-aa[1:3,]
m<-apply(aaa,2,mean)
min(m)
max(m)






  ############# Ts/Tv ratio
setwd("~/JICAutumn2016/finalAnalysis29Apr/results/TsTv ratios/")

s1<-read.table("scaf1.txt",h=F)
s2<-read.table("scaf2.txt",h=F)
s3<-read.table("scaf3.txt",h=F)
s4<-read.table("scaf4.txt",h=F)
s5<-read.table("scaf5.txt",h=F)
s6<-read.table("scaf6.txt",h=F)
s7<-read.table("scaf7.txt",h=F)
s8<-read.table("scaf8.txt",h=F)
s1$prop<-s1$V5/(s1$V4+s1$V5)
s2$prop<-s2$V5/(s2$V4+s2$V5)
s3$prop<-s3$V5/(s3$V4+s3$V5)
s4$prop<-s4$V5/(s4$V4+s4$V5)
s5$prop<-s5$V5/(s5$V4+s5$V5)
s6$prop<-s6$V5/(s6$V4+s6$V5)
s7$prop<-s7$V5/(s7$V4+s7$V5)
s8$prop<-s8$V5/(s8$V4+s8$V5)

aa<-cbind(s1$prop,s2$prop,s3$prop,s4$prop,s5$prop,s6$prop,s7$prop,s8$prop)
tot<-apply(aa,1,mean)

prot<-c("ASY1","ASY3","DYAD","PDS5B","PRD3","SCC4","SDS","SHOC1","SMG7","SYN1","ZYP1A","ZYP1B")
for (p in prot) { #    p = "SMG7"
a<-read.table(paste(p,".tstv_counts.txt",sep=""),h=T)
a1<-a[,2:3]
pr<-prop.test(x = a1$X.7.nTransitions,p = tot,n = a1$X.7.nTransitions+a1$X.8.nTransversions)

tot1<-mean(tot)
pr<-prop.test(x = 21,p = tot1,n = 42)


print(pr$p.value)

}


###CANDIDATE PROTEINS

A<-read.table("file:///home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/tetraploid/all_dip_tet_0.99_overlapsSGF.txt",h=T)
aa<-as.data.frame(table(A$ALcode))
aaa<-subset(aa,aa$Freq>2)
write.table(aaa,"file:///home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/tetraploid/candidates_SGF_3moreAAS.txt",col.names = F,row.names = F,quote = F)


aa<-read.table("file:///home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/tetraploid/candidates_SGF_3moreAAS.txt")
bb<-read.table("/home/aa/Desktop/references/lyrataV2/genesIDsLyV2.gff")
bb$V10<-substr(bb$V9,4,12)
cc<-subset(bb,V10 %in% aa$V1)
cc$V11<-paste(cc$V1,":",cc$V4,"-",cc$V5,sep="")
write.table(cc$V11,"file:///home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/tetraploid/candidates_SGF_3moreAAS.intervals",col.names = F,row.names = F,quote = F)

# ANNOTATE
library(data.table)
aa<-read.table("file:///home/aa/JICAutumn2016/finalAnalysis29Apr/results/candidates_SGF_3moreAAS.txt")
ann<-fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20181231.txt",h=T,quote="")
mm3<-subset(ann,ann$AL %in% aa$V1)
mm<-as.data.frame(mm3[!duplicated(mm3[,c(1)]),])

ann<-fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR10orth_des_20150927.txt",h=T)
mm3<-subset(ann,ann$`Version-2` %in% aa$V1)
mm<-as.data.frame(mm3[!duplicated(mm3[,c(1)]),])
write.table(mm3,"file:///home/aa/JICAutumn2016/finalAnalysis29Apr/results/candidates_SGF_3moreAAS_withAnnotations.txt",quote = F,col.names = F,row.names = F)

# HAMMING SUMMARY
aa<-read.table("file:///home/aa/JICAutumn2016/finalAnalysis29Apr/results/Hamming_all candidates/special_genes.txt")
aa1<-subset(aa,aa$V2<aa$V4)
fisher.test(matrix(c(4, 6, 52, 56), nrow=2))
fisher.test(matrix(c(4, 52, 6, 56), nrow=2))

56/108
62/118

  #### DOTPLOT
library(data.table)
aa<-fread("file:///home/aa/JICAutumn2016/finalAnalysis29Apr/data/subsampled_lyrataMapped/TETDIP_WS1_MS1_BPM.txt",h=T)

ge<-"AL4G46460"  #AL6G15380
a1<-subset(aa,aa$ALcode %in% ge)
a<-subset(aa,aa$scaff %in% a1$scaff[1] & aa$end > (min(a1$end)-6000) & aa$end < (max(a1$end)+6000))

bb<-read.table("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/FinalListTetraploids7Aug.txt",h=T)
b<-subset(bb,bb$ALcode %in% ge)



    cc<-ggplot() + 
      geom_point(data = a, aes(x=as.numeric(as.character(end)), y=abs(as.numeric(as.character(AFD)))),size = 3, alpha = 0.5) +
      geom_point(data = b, aes(x=as.numeric(as.character(end)), y=abs(as.numeric(as.character(AFD)))),size = 3, alpha = 0.6,col="red3") + 
      geom_rect(aes(xmin = min(a1$end), xmax = max(a1$end), ymin = -Inf, ymax = Inf),alpha=0.1) + 
      labs(x=paste("Position on ",a1$scaff[1],sep=""), y ='AFD (diploids - tetraploids)',title = a1$ALcode[1]) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size=12),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))
    pdf(paste("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/",ge,".pdf",sep=""),height = 4,width = 8,pointsize = 12)
    cc
    dev.off()
    
    
    ### Number of candidate meiosis 
    library(data.table)
    a<-fread("/home/aa/JICAutumn2016/finalAnalysis29Apr/data/SEC.table.repol.txt",h=F)
    a<-subset(a,a$V8 %in% "missense_variant")
    ag<-as.data.frame(table(a$V7))
    ag1<-subset(ag, ag$Freq>1)
    
s=4
p<-read.table("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/diploid/Pannonian_0.99_overlapsSGF.txt",h=T)
pg<-as.data.frame(table(substr(p$ALcode,1,9)))
pg1<-subset(pg, pg$Freq > 1)
fisher.test(matrix(c(20000-78-nrow(pg1)-s,78-s , nrow(pg1)-s, s), nrow=2))

s=0
p<-read.table("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/diploid/Dinaric_0.99_overlapsSGF.txt",h=T)
pg<-as.data.frame(table(substr(p$ALcode,1,9)))
pg1<-subset(pg, pg$Freq>1)
fisher.test(matrix(c(20000-78-nrow(pg1)-s,78-s , nrow(pg1)-s, s), nrow=2))

s=2
p<-read.table("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/diploid/Secarp_0.99_overlapsSGF.txt",h=T)
pg<-as.data.frame(table(substr(p$ALcode,1,9)))
pg1<-subset(pg, pg$Freq>1)
fisher.test(matrix(c(20000-78-nrow(pg1)-s,78-s , nrow(pg1)-s, s), nrow=2))

s=3
p<-read.table("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/diploid/Baltic_0.99_overlapsSGF.txt",h=T)
pg<-as.data.frame(table(substr(p$ALcode,1,9)))
pg1<-subset(pg, pg$Freq>1)
fisher.test(matrix(c(20000-78-nrow(pg1)-s,78-s , nrow(pg1)-s, s), nrow=2))

s=0
p<-read.table("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/diploid/Wcarp_0.99_overlapsSGF.txt",h=T)
pg<-as.data.frame(table(substr(p$ALcode,1,9)))
pg1<-subset(pg, pg$Freq>1)
fisher.test(matrix(c(20000-78-nrow(pg1)-s,78-s , nrow(pg1)-s, s), nrow=2))

s=11
p<-read.table("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/tetraploid/all_dip_tet_0.99_overlapsSGF.txt",h=T)
pg<-as.data.frame(table(substr(p$ALcode,1,9)))
pg1<-subset(pg, pg$Freq>1)
fisher.test(matrix(c(20000-78-nrow(pg1)-s,78-s , nrow(pg1)-s, s), nrow=2))

###No. AAS per candidate meiosis/no meiosis
#CANDIDATE PROTEINS
setwd("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/tetraploid/")
A<-read.table("file:///home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/tetraploid/all_dip_tet_0.99_overlapsSGF.txt",h=T)
aa<-as.data.frame(table(A$ALcode))
aaa<-subset(aa,aa$Freq>2)
m<-read.table("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/tetraploid/all_dip_tet_0.99_overlapsSGF.meiosis.txt",h=T)
mm<-as.data.frame(table(m$ALcode))
mmm<-subset(mm,mm$Freq>2)
nm<-subset(aaa,!aaa$Var1 %in% mmm$Var1)
zyp <- data.frame("ZYP", 25)      
names(zyp) <- c("Var1", "Freq")  
mmm1 <- rbind(mmm, zyp) 
shoc <- data.frame("SHOC", 41)      
names(shoc) <- c("Var1", "Freq")  
mmm1 <- rbind(mmm1, shoc) 
mmm2<-mmm[-c(2,7,9),]
# 9 SGF
summary(mmm)
# 11 all
summary(mmm1)
# 6 denovo
summary(mmm2)



library(ggpubr)

novA <- cbind("Meiotic",mmm2$Freq)
novT <- cbind("Non-meiotic",nm$Freq)
all<-as.data.frame(rbind(novA,novT))
#Only diploid-tetraploid
all$V2<-as.numeric(as.character(all$V2))
my_comparisons <- list( c("Meiotic", "Non-meiotic") )
pdf("NoAASs_6Denovo.pdf",width = 3.5,height = 5,pointsize = 24)
ggviolin(all, x = "V1", y = "V2", fill = "V1",
         palette = c("gold", "grey"),
         add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = 1.5)
dev.off()



# Rarefaction
library(data.table)
library(Publish)
library(plotrix)
setwd("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/")
m<-read.table("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/tetraploid/all_dip_tet_0.99_overlapsSGF.txt",h=T)
d<-fread("DIP.table.repol.txt",h=F,na.strings = "-9")
#head(d)
d1<-subset(d, paste(d$V3,d$V4) %in% paste(m$scaff,m$end))
nrep<-100
nind<-80
seq<-c(seq(1,20,1),seq(25,80,5))
mat <- matrix(nrow = length(seq), ncol = nrep+4)
mat[,1]<-seq
for (i in 1:nrow(mat)) {  #  i = 32
  for (num in 1:nrep) {
    

  subs<-mat[i,1]
  aaa<-sample.int(n = nind, size = subs, replace = F)  
  aaa<-10+aaa
  d2<-d1[,..aaa]
  d2$sum<-rowSums(d2,na.rm = T)
  s<-subset(d2,d2$sum == 0)
  mat[i,num+1]<-nrow(d2)-nrow(s)

} 

  ci<-ci.mean(mat[i,2:(ncol(mat)-3)])
  mat[i,ncol(mat)-2]<-ci$lower
  mat[i,ncol(mat)-1]<-ci$upper
  mat[i,ncol(mat)-0]<-ci$mean
  }

pdf("rarefactionDiploids_80.pdf",width = 7,height = 4,pointsize = 10)
plotCI(mat[,1], y = mat[,ncol(mat)-0], ui=mat[,ncol(mat)-1], li=mat[,ncol(mat)-2],xlab="Number of individuals sampled",ylab="Number of detected segregating sites")
dev.off()



###PAI - number of AASs
setwd("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/")
t<-read.table("posSelection/FinalListTetraploids7Aug.txt", h=T)
a<-read.table("file:///home/aa/JICAutumn2016/finalAnalysis29Apr/pai_Malvidae/Alignment_Identities_allSites.txt",h=F)
l<-read.table("../lists/ALmeiotic.ALcode.tab.Name.txt")
l$pai<-""
l$aas<-""
for (i in 1:nrow(l)) { #   i=1
  l1<-as.character(droplevels(l[i,1]))
  a1<-subset(a,a$V6 %in% l1)
  l$pai[i]<-mean(a1$V3)
  t1<-subset(t,t$ALcode %in% l1)
  l$aas[i]<-nrow(t1)
  
}

#plot(l$aas~l$pai)

plot(l$pai~l$aas)
#abline(lm(l$pai~l$aas))
pdf("PAI_AAS.pdf",5,5)
plot(x = l$aas, y = l$pai, xlab= "Number of tetraploid-differentiated AASs", ylab= "Alignment identity", col= "blue", pch = 19, cex = 1, lty = "solid", lwd = 2)
text(x = l$aas, y = as.numeric(l$pai), labels=l$V2,col = "black",pos = 4)
abline(h=mean(as.numeric(l$pai)))
dev.off()


#Subsampling 16-16
library(data.table)
aa<-read.table("file:///home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/FinalListTetraploids7Aug.txt",h=T)
bb<-fread("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/TETDIP_WS1_MS1_BPM.txt",h=T)
b1<-subset(bb,bb$ann %in% "missense_variant")
b2<-b1[ order(b1[,20],decreasing = T), ]
b3<-b2[1:10861,]
b4<-subset(b3,paste(b3$scaff,b3$AAS) %in% paste(aa$scaff,aa$AAS))
b5<-b4[ order(b4[,7],decreasing = F), ]
table(b5$ALcode)

a1<-subset(aa,!paste(aa$scaff,aa$AAS) %in% paste(b4$scaff,b4$AAS))
a2<-a1[1:7,]
write.table(a2,"/home/aa/JICAutumn2016/finalAnalysis29Apr/results/missingSubsampling16_16.txt",quote = F,row.names = F)




library(data.table)
a<-fread("file:///home/aa/JICAutumn2016/finalAnalysis29Apr/data/subsampled_lyrataMapped/TETDIP.synNon_WS1_MS1_BPM.txt",h=T)
s<-subset(a,a$ann %in% "synonymous_variant")
ss<-subset(s,abs(s$AFD) == 1)
n<-subset(a,a$ann %in% "missense_variant")
nn<-subset(n,abs(n$AFD) == 1)
hist(abs(s$AFD), breaks = 100, main = "AFD diploids-tetraploids synonymous")
hist(abs(n$AFD), breaks = 100,main = "AFD diploids-tetraploids nonsynonymous")
summary(abs(s$AFD))
quantile(abs(s$AFD),0.99)
summary(abs(n$AFD))
quantile(abs(n$AFD),0.99)


###MNM

setwd("/home/aa/JICAutumn2016/finalAnalysis29Apr/data/")
library(data.table)

CHO_001_c - 16
BRD_001_a - 10
HOC_004_1 - 25
SWA_002_1 - 47


a<-fread("TET.table.repol.txt",h=F)
a1<-subset(a,a$V8 %in% "missense_variant")
a1<-a
a2<-subset(a1,!a1$V47 %in% "-9")
a3<-subset(a2,!a2$V47 %in% 0)
a3<-subset(a2,a2$V47 %in% 4)

setkey(a3, V7, V4)
a3[, diff := V4 - shift(V4, fill = first(V4)), by = V7]
a4<-subset(a3,!a3$diff %in% 0)
a5<-subset(a4,!a4$V7 %like% "-")
a6<-subset(a5,a5$diff < 9255)

#   write.table(a5[,c(1,2,3,4,5,6,7,8,52)],"SNPdistances.txt",quote = F,row.names = F,col.names = F)

pdf("diff_allSNPs_SWA_002_1.fixed.pdf")
hist(a5$diff,breaks = 20000,xlab = "distance (bp)",main = "Distance between adjescent missense SNPs within a gene",xlim = c(0,10000))
hist(a5$diff,breaks = 100000,xlab = "distance (bp)",main = "Distance between adjescent missense SNPs within a gene",xlim = c(0,2000))
hist(a5$diff,breaks = 400000,xlab = "distance (bp)",main = "Distance between adjescent missense SNPs within a gene",xlim = c(0,500))
dev.off()



for(p in c("AL1G10680","AL1G35730","AL4G46460","AL6G30890","AL8G25590","AL8G25600","AL2G25920","AL2G37810")){ # p="AL4G46460"  p="AL6G30890"
  
  b2<-subset(a6,as.character(a6$V7) %in% p)
  pdf(paste("boxplot.",p,".fixed.pdf"))
  boxplot(a6$diff, b2$diff)
  print(summary(a6$diff))
  print(summary(b2$diff))
  hist(b2$diff,nclass = 100)
  
  dev.off()
  print(p)
  print(wilcox.test(a6$diff, b2$diff))
}

b3<-subset(b1,!b1$diff %in% 0)
pdf(paste("hist.all.fixed.pdf"))
boxplot(a4$diff, b3$diff,names = c('genome-wide',"meiosis FstH candidates"),ylab="distance (bp)")
dev.off()
print(wilcox.test(a4$diff, b3$diff))
















b<-fread("High_Fst_TET.txt")
b<-fread("../results/haplotypes/2HaplDipTet/tetraploidHIF.intervals",sep = ":",h=F)

b1<-subset(a1,paste(a1$V7,a1$V9) %in% paste(b$ALcode,b$AASS))
b1<-subset(a1,paste(a1$V3,a1$V4) %in% paste(b$V1,b$V2))

setkey(b1, V7, V4)
b1[, diff := V4 - shift(V4, fill = first(V4)), by = V7]
b1$prot<-b$Protein

for(p in unique(b1$V7)){
  
  b2<-subset(b1,as.character(b1$V7) %in% p)
  b3<-subset(b2,!b2$diff %in% 0)
  pdf(paste("hist.",p,".fixed.pdf"))
  boxplot(a4$diff, b3$diff)
  dev.off()
  print(p)
  print(wilcox.test(a4$diff, b3$diff))
}

b3<-subset(b1,!b1$diff %in% 0)
pdf(paste("hist.all.fixed.pdf"))
boxplot(a4$diff, b3$diff,names = c('genome-wide',"meiosis FstH candidates"),ylab="distance (bp)")
dev.off()
print(wilcox.test(a4$diff, b3$diff))


