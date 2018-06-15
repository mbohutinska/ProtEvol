# ProtEvol

###synNon analysis
#Input file:
/holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/bordel/All.inclAdded.synonNonsyn.NOfiltr.vcf.gz
#incl. missing data, rare variants, only nonsyn and synon sites, all the 300 individuals
#devided into chromosomes: (/holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/variants.sh)
```
#!/bin/bash -e
#PBS -N SelectVar
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=2:mem=8gb:scratch_local=25gb
#PBS -m abe
#PBS -j oe

module add gatk-3.7 
trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 
DATADIR="/storage/plzen1/home/holcovam/ScanTools/synNon"
cp /storage/plzen1/home/holcovam/references/lyrataV2/alygenomes* $SCRATCHDIR || exit 1
cp $DATADIR/*.intervals $SCRATCHDIR || exit 1
cp $DATADIR/All.inclAdded.synonNonsyn.NOfiltr.vcf.gz* $SCRATCHDIR || exit 1
cd $SCRATCHDIR || exit 2
echo data loaded at `date`

##FUNGUJE:
#java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar \
#  -T CombineVariants \
#  -R alygenomes.fasta \
#  --variant:non added7March.missense.vcf.gz \
#  --variant:syn All.missense.vcf.gz \
#  -o All.inclAdded.missense.vcf.gz \
#  -genotypeMergeOptions PRIORITIZE \
#  -priority non,syn
#rm All.missense.vcf.gz*
#
#java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar \
#  -T CombineVariants \
#  -R alygenomes.fasta \
#  --variant:non added7March.synon.vcf.gz \
#  --variant:syn All.synon.vcf.gz \
#  -o All.inclAdded.synon.vcf.gz \
#  -genotypeMergeOptions PRIORITIZE \
#  -priority non,syn
#rm added7March.*
#rm All.synon.vcf.gz*
#
#java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar \
#  -T CombineVariants \
#  -R alygenomes.fasta \
#  --variant:non All.inclAdded.missense.vcf.gz \
#  --variant:syn All.inclAdded.synon.vcf.gz \
#  -o All.inclAdded.synonNonsyn.NOfiltr.vcf.gz \
#  -genotypeMergeOptions PRIORITIZE \
#  -priority non,syn

java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.inclAdded.synonNonsyn.NOfiltr.vcf.gz -o All.inclAdded.synonNonsyn.NOfiltr.scaff1.vcf.gz -L s1.intervals
java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.inclAdded.synonNonsyn.NOfiltr.vcf.gz -o All.inclAdded.synonNonsyn.NOfiltr.scaff2.vcf.gz -L s2.intervals
java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.inclAdded.synonNonsyn.NOfiltr.vcf.gz -o All.inclAdded.synonNonsyn.NOfiltr.scaff3.vcf.gz -L s3.intervals
java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.inclAdded.synonNonsyn.NOfiltr.vcf.gz -o All.inclAdded.synonNonsyn.NOfiltr.scaff4.vcf.gz -L s4.intervals
java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.inclAdded.synonNonsyn.NOfiltr.vcf.gz -o All.inclAdded.synonNonsyn.NOfiltr.scaff5.vcf.gz -L s5.intervals
java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.inclAdded.synonNonsyn.NOfiltr.vcf.gz -o All.inclAdded.synonNonsyn.NOfiltr.scaff6.vcf.gz -L s6.intervals
java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.inclAdded.synonNonsyn.NOfiltr.vcf.gz -o All.inclAdded.synonNonsyn.NOfiltr.scaff7.vcf.gz -L s7.intervals
java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.inclAdded.synonNonsyn.NOfiltr.vcf.gz -o All.inclAdded.synonNonsyn.NOfiltr.scaff8.vcf.gz -L s8.intervals

rm *.intervals
rm All.inclAdded.synonNonsyn.NOfiltr.vcf.gz*
rm alygenomes*
cp $SCRATCHDIR/* $DATADIR/ || export CLEAN_SCRATCH=false
printf "\nFinished\n\n"
```

#output:
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/All.inclAdded.synonNonsyn.NOfiltr.scaff8.vcf.gz.tbi
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/All.inclAdded.synonNonsyn.NOfiltr.scaff8.vcf.gz
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/All.inclAdded.synonNonsyn.NOfiltr.scaff7.vcf.gz.tbi
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/All.inclAdded.synonNonsyn.NOfiltr.scaff7.vcf.gz
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/All.inclAdded.synonNonsyn.NOfiltr.scaff6.vcf.gz.tbi
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/All.inclAdded.synonNonsyn.NOfiltr.scaff6.vcf.gz
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/All.inclAdded.synonNonsyn.NOfiltr.scaff5.vcf.gz.tbi
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/All.inclAdded.synonNonsyn.NOfiltr.scaff5.vcf.gz
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/All.inclAdded.synonNonsyn.NOfiltr.scaff4.vcf.gz.tbi
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/All.inclAdded.synonNonsyn.NOfiltr.scaff4.vcf.gz
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/All.inclAdded.synonNonsyn.NOfiltr.scaff3.vcf.gz.tbi
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/All.inclAdded.synonNonsyn.NOfiltr.scaff3.vcf.gz
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/All.inclAdded.synonNonsyn.NOfiltr.scaff2.vcf.gz.tbi
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/All.inclAdded.synonNonsyn.NOfiltr.scaff2.vcf.gz
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/All.inclAdded.synonNonsyn.NOfiltr.scaff1.vcf.gz.tbi
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/synNon/All.inclAdded.synonNonsyn.NOfiltr.scaff1.vcf.gz

#ScanTools
###Protein evolution
module add python34-modules-gcc
python3
import ScanTools
test = ScanTools.scantools("/storage/plzen1/home/holcovam/ScanTools") 
	#necessary to paste in all values as strings "" instead of intagers!! double-check if you have ALL individuals specified in the PopKey!!!!
#1.			For whole genome
test.splitVCFsAnn(vcf_dir="synNon",repolarization_key="repolarized.lookupKey.perSpeciesThreshold.txt", min_dp="8",mffg="0.5", mem="8", time_scratch='03:00:00', ncpu="2",overwrite=True, scratch_gb="1",keep_intermediates=False, use_scratch=True,scratch_path="$SCRATCHDIR", pops=['CRO','PAN','DIN','BAL','SEC','WCA','TET'], print1=False)

scp holcovam@nympha.metacentrum.cz:/storage/plzen1/home/holcovam/ScanTools/VCF_synNon_DP8.M1/*.table.repol.txt .
#finish killed:
test.splitVCFsAnn(vcf_dir="synNon",repolarization_key="repolarized.lookupKey.perSpeciesThreshold.txt", min_dp="8",mffg="0.5", mem="8", time_scratch='03:00:00', ncpu="2",overwrite=True, scratch_gb="1",keep_intermediates=False, use_scratch=True,scratch_path="$SCRATCHDIR", pops=['BAL','TET'], print1=False)

#to delete all completely missing sites: (otherwise devision by 0 error) - does not help much.. - better to filter for missingnes < 0.5?? - To discuss
grep -v  -e "-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9." SEC
grep -v  -e "-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9." BAL
grep -v  -e "-9.-9.-9.-9." CRO
grep -v  -e "-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9." DIN
grep -v  -e "-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9." WCA
grep -v  -e "-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9.-9." TET

#for a single pair
#test.calcbpmAnn(recode_dir= "VCF_bordel_DP8.M0.1", pops=['PAN','TET'], output_name="PAN_TET", window_size=1, min_snps=1, mem=4, ncpu=1, use_repol=True, keep_intermediates=False, time_scratch="0:30:00",scratch_gb=1, print1=False)
#2.			pairwise - timinf fine-tuned
test.calcPairwisebpmAnn(recode_dir= "VCF_synNon_DP8.M0.5", pops=['CRO','PAN','DIN','BAL','SEC','WCA','TET'], window_size=1, min_snps=1, mem=1, ncpu=1, use_repol=True, keep_intermediates=False, time_scratch="4:00:00",scratch_gb=1, print1=False)
#3.			download and run locally
#/home/aa/JICAutumn2016/finalAnalysis29Apr/data/

#4.			make missing annotation files - once, then use is for all analyses with current data
#get list af all genes in lyrata
cat LyV2.gff | grep -v "#" | grep -v "gene" | cut -f3 -d '=' | sort | uniq | grep -v "U" > ALcodesAll.txt
#holcovam/ScanTools/synNon/bordel/snpSIFT.sh - to get AAS - genomic position dictionary

```
#!/bin/bash -e
#PBS -N snpSIFT
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=2:mem=8gb:scratch_local=1gb
#PBS -m abe
#PBS -j oe

module add jdk-8
trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 
DATADIR="/storage/plzen1/home/holcovam/ScanTools/synNon/bordel"
cp /storage/plzen1/home/holcovam/programs/snpEff/SnpSift.jar $SCRATCHDIR || exit 1
cp $DATADIR/All.inclAdded.synonNonsyn.NOfiltr.vcf.gz* $SCRATCHDIR || exit 1
cd $SCRATCHDIR || exit 2
echo data loaded at `date`

java -Xmx8g -jar SnpSift.jar extractFields -s "," -e "." All.inclAdded.synonNonsyn.NOfiltr.vcf.gz CHROM POS ANN[*].GENE ANN[*].HGVS_P > Pos_AAS_300AllSynNon.tab.txt

rm All.inclAdded.synonNonsyn.NOfiltr.vcf.gz*
rm SnpSift.jar
cp $SCRATCHDIR/* $DATADIR/ || export CLEAN_SCRATCH=false
printf "\nFinished\n\n"
```

##run R script /home/aa/JICAutumn2016/finalAnalysis29Apr/protScanTools.r
#outputs:
#MK test
#highFst AASs 

####NEW - to get number of outlier high-Fst SNPs genome-wide
test.outlierAASs(contrasts=str('"WCATET", "SECWCA", "SECTET", "PANWCA", "PANTET", "PANSEC", "PANDIN", "PANBAL", "DINWCA", "DINTET", "DINSEC", "DINBAL", "BALWCA", "BALTET", "BALSEC"'), recode_dir = "VCF_synNon_DP8.M0.5", outlier_quantile = "0.99", alcode_path = "../references/lyrataV2/ALcodesAll.txt", time_scratch="48:00:00", ncpu=1, mem=4, scratch_gb = 4, print1 = True)
##more cpu and ram - i.e. 8 ncpu, 24 RAM??, 1 cpu, 4 ram - 39 hours for .99 outliers
#from interval list (from meiotic.highFstFinalSummaryAnn.txt - not earlier!) run heatmap pipeline
#ScanTools/synNon/bordel/variants.sh:

```
#!/bin/bash -e
#PBS -N SynNon
#PBS -l walltime=05:00:00
#PBS -l select=1:ncpus=4:mem=16gb:scratch_local=1gb
#PBS -m abe
#PBS -j oe

module add gatk-3.7 
trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 
DATADIR="/storage/plzen1/home/holcovam/ScanTools/synNon/bordel"
cp /storage/plzen1/home/holcovam/references/lyrataV2/alygenomes* $SCRATCHDIR || exit 1
cp $DATADIR/quantile99.intervals $SCRATCHDIR || exit 1
cp $DATADIR/All.inclAdded.missense.vcf.gz* $SCRATCHDIR || exit 1
cd $SCRATCHDIR || exit 2
echo data loaded at `date`

java -Xmx16g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.inclAdded.missense.vcf.gz -o MeioticHighFst.inclAdded.missense.vcf.gz -L quantile99.intervals

rm *.intervals
rm All.inclAdded.missense.vcf.gz*
rm alygenomes*
cp $SCRATCHDIR/* $DATADIR/ || export CLEAN_SCRATCH=false
printf "\nFinished\n\n"

```

###problem with annotation
aa@majda:~/JICAutumn2016/finalAnalysis29Apr/data> cat tet.txt dic.txt | sort | uniq -d | wc -l
5494407
aa@majda:~/JICAutumn2016/finalAnalysis29Apr/data> wc -l tet.txt 
5499857 tet.txt
aa@majda:~/JICAutumn2016/finalAnalysis29Apr/data> wc -l dic.txt 
5996693 dic.txt

##Get SNPeff annotation
version 4_2 - Aly1 database - for genes not having AL1GXXXX code
verion 4_3 - can build LyV2 database (bug in 4_2)
1. BUILD LYV2 DATABASE: version 4_3
java -jar snpEff.jar build -v LyV2 
I put this:
# Arabidopsis lyrata V2
LyV2.genome : Lyrata
under these lines in config file:
# Arabidopsis lyrata
alyrata107.genome : Arabidopsis_lyrata
alyrata107.reference : http://www.phytozome.net/search.php?method=Org_Alyrata
alyrata1.genome : Arabidopsis_lyrata
alyrata1.reference : http://genome.jgi-psf.org/Araly1/Araly1.download.ftp.html

Than you need to put LyV2.fa into snpEff/data/genomes folder and genes.gff into snpEff/data/LyV2 folder.
Then you just use the command to create database

#####################21st of May - after Pirita
holcovam@nympha.metacentrum.cz/auto/plzen1/home/holcovam/ScanTools/all300/ - non-annotated files, all SNPs, no filtration
snpEff.sh: to annotate

###piS
ScanTools/all300/bordel/filterAN_4d.sh:
###ty blbecku musi tam byt i nonvariable!!!
```
#!/bin/bash -e
#PBS -N SelectVar4d
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=2:mem=8gb:scratch_local=19gb
#PBS -m abe
#PBS -j oe
module add gatk-3.7 
trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 
DATADIR="/storage/plzen1/home/holcovam/ScanTools/all300/bordel/"
cp /storage/plzen1/home/holcovam/references/lyrataV2/alygenomes* $SCRATCHDIR || exit 1
cp $DATADIR/*GATK.intervals $SCRATCHDIR || exit 1
cp $DATADIR/All.NOfiltr.scaff* $SCRATCHDIR || exit 1
cd $SCRATCHDIR || exit 2
echo data loaded at `date`
java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.NOfiltr.scaff1.vcf.gz -o All.NOfiltr.4dg.scaff1.vcf.gz -L four.fold.degenerate.sites.GATK.intervals
java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.NOfiltr.scaff2.vcf.gz -o All.NOfiltr.4dg.scaff2.vcf.gz -L four.fold.degenerate.sites.GATK.intervals
java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.NOfiltr.scaff3.vcf.gz -o All.NOfiltr.4dg.scaff3.vcf.gz -L four.fold.degenerate.sites.GATK.intervals
java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.NOfiltr.scaff4.vcf.gz -o All.NOfiltr.4dg.scaff4.vcf.gz -L four.fold.degenerate.sites.GATK.intervals
java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.NOfiltr.scaff5.vcf.gz -o All.NOfiltr.4dg.scaff5.vcf.gz -L four.fold.degenerate.sites.GATK.intervals
java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.NOfiltr.scaff6.vcf.gz -o All.NOfiltr.4dg.scaff6.vcf.gz -L four.fold.degenerate.sites.GATK.intervals
java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.NOfiltr.scaff7.vcf.gz -o All.NOfiltr.4dg.scaff7.vcf.gz -L four.fold.degenerate.sites.GATK.intervals
java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.NOfiltr.scaff8.vcf.gz -o All.NOfiltr.4dg.scaff8.vcf.gz -L four.fold.degenerate.sites.GATK.intervals
rm *.intervals
rm All.NOfiltr.scaff*
rm alygenomes*
cp $SCRATCHDIR/* $DATADIR/ || export CLEAN_SCRATCH=false
printf "\nFinished\n\n"
```
#instead: scp /run/media/aa/Transcend/4dgData/4d_DM_HM_BP_BI_all_PASS.vcf.gz* holcovam@nympha.metacentrum.cz:~/ScanTools/300_fourfold


#TODO
cp to all3004d
test.splitVCFs(vcf_dir="all3004d",repolarization_key="repolarized.lookupKey.perSpeciesThreshold.txt", min_dp="8",mffg="0.5", mem="16", time_scratch='3:00:00', ncpu="2",overwrite=True, scratch_gb="2",keep_intermediates=False, use_scratch=True,scratch_path="$SCRATCHDIR", pops=['PAN','DIN','BAL','SEC','WCA','TET','DIP'], print1=False)

test.calcwpm(recode_dir = "VCF_all3004d_DP8.M0.5", window_size = 50000, min_snps = 1, pops=['PAN','DIN','BAL','SEC','WCA'], mem=1, ncpu=1, scratch_gb=1, use_repol=True, time_scratch="1:00:00", overwrite=False, sampind="15", print1=False)
test.calcwpm(recode_dir = "VCF_all3004d_DP8.M0.5", window_size = 50000, min_snps = 1, pops=['TET','DIP'], mem=1, ncpu=1, scratch_gb=1, use_repol=True, time_scratch="2:00:00", overwrite=False, sampind="-99", print1=False)

#only genome>
cat BAL.WS50.0k_MS1_15ind_WPM.txt | head -n1 > genome.300subsampled.WPM.txt 
cat PAN.WS50.0k_MS1_15ind_WPM.txt | grep "Genome" >> genome.300subsampled.WPM.txt
cat DIN.WS50.0k_MS1_15ind_WPM.txt | grep "Genome" >> genome.300subsampled.WPM.txt
cat SEC.WS50.0k_MS1_15ind_WPM.txt | grep "Genome" >> genome.300subsampled.WPM.txt
cat BAL.WS50.0k_MS1_15ind_WPM.txt | grep "Genome" >> genome.300subsampled.WPM.txt
cat WCA.WS50.0k_MS1_15ind_WPM.txt | grep "Genome" >> genome.300subsampled.WPM.txt
cat DIP.WS50.0k_MS1_39ind_WPM.txt | grep "Genome" >> genome.300subsampled.WPM.txt
cat TET.WS50.0k_MS1_39ind_WPM.txt | grep "Genome" >> genome.300subsampled.WPM.txt

cat BAL.WS50.0k_MS100_15ind_WPM.txt | head -n1 >>      genome.300subsampled.WPM.txt 
cat PAN.WS50.0k_MS100_15ind_WPM.txt | grep "Genome" >> genome.300subsampled.WPM.txt
cat DIN.WS50.0k_MS100_15ind_WPM.txt | grep "Genome" >> genome.300subsampled.WPM.txt
cat SEC.WS50.0k_MS100_15ind_WPM.txt | grep "Genome" >> genome.300subsampled.WPM.txt
cat BAL.WS50.0k_MS100_15ind_WPM.txt | grep "Genome" >> genome.300subsampled.WPM.txt
cat WCA.WS50.0k_MS100_15ind_WPM.txt | grep "Genome" >> genome.300subsampled.WPM.txt
cat DIP.WS50.0k_MS100_39ind_WPM.txt | grep "Genome" >> genome.300subsampled.WPM.txt
cat TET.WS50.0k_MS100_39ind_WPM.txt | grep "Genome" >> genome.300subsampled.WPM.txt

###1.Subsampled, all variants
#screen
test.splitVCFsAnn(vcf_dir="all300",repolarization_key="repolarized.lookupKey.perSpeciesThreshold.txt", min_dp="8",mffg="0.5", mem="24", time_scratch='13:00:00', ncpu="2",overwrite=True, scratch_gb="4",keep_intermediates=False, use_scratch=True,scratch_path="$SCRATCHDIR", pops=['CRO','PAN','DIN','BAL','SEC','WCA','TET'], print1=False)
#crashed
test.splitVCFsAnn(vcf_dir="all300",repolarization_key="repolarized.lookupKey.perSpeciesThreshold.txt", min_dp="8",mffg="0.5", mem="12", time_scratch='13:00:00', ncpu="2",overwrite=True, scratch_gb="20",keep_intermediates=False, use_scratch=True,scratch_path="$SCRATCHDIR", pops=['BAL','SEC','WCA','TET','DIP'], print1=False)

#TODO
###add diploid to PopKey to filtered and full
#download and run locally
#/home/aa/JICAutumn2016/finalAnalysis29Apr/data/
##run R script /home/aa/JICAutumn2016/finalAnalysis29Apr/protScanTools.r
#outputs:
#N/S
scp holcovam@nympha.metacentrum.cz:/storage/plzen1/home/holcovam/ScanTools/VCF_all300_DP8.M0.5/*.table.repol.txt .

#requires>1gb memory!!
#to test
test.N_SperGene(lineages=str('"CRO"'), recode_dir = "VCF_all300_DP8.M0.5", alcode_path = "../references/lyrataV2/ALcodesAll.txt", time_scratch="00:10:00", ncpu=1, mem=1, scratch_gb = 1, print1 = True)
#final comand
test.N_SperGene(lineages=str('"CRO","PAN","DIN","BAL","SEC","WCA","DIP","TET"'), recode_dir = "VCF_all300_DP8.M0.5", alcode_path = "../references/lyrataV2/ALcodesAll.txt", time_scratch="48:00:00", ncpu=1, mem=16, scratch_gb = 1, print1 = False)
#to test if meiosis genes are significantly higher/lower N/S - Pearson’s Chi-squared test 

#Fst_H
test.calcPairwisebpmAnn(recode_dir= "VCF_all300_DP8.M0.5", pops=['TET','DIP'], window_size=1, min_snps=1, mem=1, ncpu=1, use_repol=True, keep_intermediates=False, time_scratch="4:00:00",scratch_gb=1, print1=False)
#max time
test.calcPairwisebpmAnn(recode_dir= "VCF_all300_DP8.M0.5", pops=['CRO','PAN','DIN','BAL','SEC','WCA'], window_size=1, min_snps=1, mem=1, ncpu=1, use_repol=True, keep_intermediates=False, time_scratch="10:00:00",scratch_gb=1, print1=False)
#max time 1:27, 15 mb RAM
scp holcovam@nympha.metacentrum.cz:/storage/plzen1/home/holcovam/ScanTools/VCF_all300_DP8.M0.5/*_BPM.txt .

####NEW - to get number of outlier high-Fst SNPs genome-wide
test.outlierAASs(pops=["PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA","CROPAN","CRODIN","CROSEC","CROBAL","CROWCA"], recode_dir = "VCF_all300_DP8.M0.5", outlier_quantile = "0.99", alcode_path = "../references/lyrataV2/ALcodesAll.txt", time_scratch="3:00:00", ncpu=1, mem=4, scratch_gb = 1, print1 = False)
test.outlierAASs(pops=["PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA","CROPAN","CRODIN","CROSEC","CROBAL","CROWCA"], recode_dir = "VCF_all300_DP8.M0.5", outlier_quantile = "0.95", alcode_path = "../references/lyrataV2/ALcodesAll.txt", time_scratch="2:00:00", ncpu=1, mem=4, scratch_gb = 1, print1 = False)
test.outlierAASs(pops=["PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA","CROPAN","CRODIN","CROSEC","CROBAL","CROWCA"], recode_dir = "VCF_all300_DP8.M0.5", outlier_quantile = "0.999", alcode_path = "../references/lyrataV2/ALcodesAll.txt", time_scratch="2:00:00", ncpu=1, mem=4, scratch_gb = 1, print1 = False)
#max 02:00:10 - 1 cpu, 4 ram - 2 ram enough!!
test.outlierAASs(pops=["PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA","CROPAN","CRODIN","CROSEC","CROBAL","CROWCA"], recode_dir = "VCF_all300_DP8.M0.5", outlier_quantile = "0.99", alcode_path = "../references/lyrataV2/ALcodesAll.txt", time_scratch="2:10:00", ncpu=1, mem=2, scratch_gb = 1, min_afd = 0.1, print1 = False)
test.outlierAASs(pops=["PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA","CROPAN","CRODIN","CROSEC","CROBAL","CROWCA"], recode_dir = "VCF_all300_DP8.M0.5", outlier_quantile = "0.999", alcode_path = "../references/lyrataV2/ALcodesAll.txt", time_scratch="2:10:00", ncpu=1, mem=2, scratch_gb = 1, min_afd = 0.1, print1 = False)
#killed:
test.outlierAASs(pops=["PANDIN","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA"], recode_dir = "VCF_all300_DP8.M0.5", outlier_quantile = "0.99", alcode_path = "../references/lyrataV2/ALcodesAll.txt", time_scratch="5:10:00", ncpu=1, mem=2, scratch_gb = 1, min_afd = 0.1, print1 = False)

test.outlierAASs(pops=["DINSEC","DINWCA","SECWCA","BALWCA"], recode_dir = "VCF_all300_DP8.M0.5", outlier_quantile = "0.999", alcode_path = "../references/lyrataV2/ALcodesAll.txt", time_scratch="5:10:00", ncpu=1, mem=2, scratch_gb = 1, min_afd = 0.1, print1 = False)



test.outlierAASs(pops=["TETDIP"], recode_dir = "VCF_all300_DP8.M0.5", outlier_quantile = "0.99", alcode_path = "../references/lyrataV2/ALcodesAll.txt", time_scratch="6:00:00", ncpu=1, mem=4, scratch_gb = 1, min_afd = 0, print1 = False)
#test.outlierAASs(pops=["TETDIP"], recode_dir = "VCF_all300_DP8.M0.5", outlier_quantile = "0.95", alcode_path = "../references/lyrataV2/ALcodesAll.txt", time_scratch="5:00:00", ncpu=1, mem=4, scratch_gb = 1, min_afd = 0.1, print1 = False)
test.outlierAASs(pops=["TETDIP"], recode_dir = "VCF_all300_DP8.M0.5", outlier_quantile = "0.999", alcode_path = "../references/lyrataV2/ALcodesAll.txt", time_scratch="6:00:00", ncpu=1, mem=4, scratch_gb = 1, min_afd = 0, print1 = False)
#test.outlierAASs(pops=["TETDIP"], recode_dir = "VCF_all300_DP8.M0.5", outlier_quantile = "0.995", alcode_path = "../references/lyrataV2/ALcodesAll.txt", time_scratch="5:00:00", ncpu=1, mem=4, scratch_gb = 1, print1 = False)
#max 04:47:01 - 1 cpu, 4 ram - 4 ram ok

#old version:
test.outlierAASs(contrasts=str('"TETDIP","PANDIN","PANSEC","PANBAL","PANWCA","DINSEC","DINBAL","DINWCA","SECWCA","BALSEC","BALWCA","CROPAN","CRODIN","CROSEC","CROBAL","CROWCA"'), recode_dir = "VCF_all300_DP8.M0.5", outlier_quantile = "0.99", alcode_path = "../references/lyrataV2/ALcodesAll.txt", time_scratch="72:00:00", ncpu=4, mem=32, scratch_gb = 4, print1 = False)
##1 cpu, 4 ram - 39 hours for .999 outliers - it's enough, it won't go to more
##4 ncpu, 32 ram - 98:37:49	.99 outliers
scp holcovam@nympha.metacentrum.cz:/storage/plzen1/home/holcovam/ScanTools/VCF_all300_DP8.M0.5/Fst* . #(data)
cd ..
bash /home/aa/JICAutumn2016/finalAnalysis29Apr/gatherFstperGene.sh

#MK test
#TODO: how to summarise it?



#from interval list (from meiotic.highFstFinalSummaryAnn.txt - not earlier!) run heatmap pipeline
#ScanTools/synNon/bordel/variants.sh:


```
#!/bin/bash -e
#PBS -N SynNon
#PBS -l walltime=05:00:00
#PBS -l select=1:ncpus=4:mem=16gb:scratch_local=1gb
#PBS -m abe
#PBS -j oe

module add gatk-3.7 
trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 
DATADIR="/storage/plzen1/home/holcovam/ScanTools/synNon/bordel"
cp /storage/plzen1/home/holcovam/references/lyrataV2/alygenomes* $SCRATCHDIR || exit 1
cp $DATADIR/quantile99.intervals $SCRATCHDIR || exit 1
cp $DATADIR/All.inclAdded.missense.vcf.gz* $SCRATCHDIR || exit 1
cd $SCRATCHDIR || exit 2
echo data loaded at `date`

java -Xmx16g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V All.inclAdded.missense.vcf.gz -o MeioticHighFst.inclAdded.missense.vcf.gz -L quantile99.intervals

rm *.intervals
rm All.inclAdded.missense.vcf.gz*
rm alygenomes*
cp $SCRATCHDIR/* $DATADIR/ || export CLEAN_SCRATCH=false
printf "\nFinished\n\n"

```

#Phytozome12.1
#keyword search with ALcodes
#By default, multiple terms are combined with OR.  Prefix each term with "+" for AND searches.  To match an exact multi-word phrase  (e.g., "heat shock"), surround the phrase in quotation marks. 
#download - peptide sequences - all but ZYP!A
#The Expect value (E) is a parameter that describes the number of hits one can "expect" to see by chance when searching a database of a particular size. It decreases exponentially as the Score (S) of the match increases. Essentially, the E value describes the random background noise. For example, an E value of 1 assigned to a hit can be interpreted as meaning that in a database of the current size one might expect to see 1 match with a similar score simply by chance.
#The lower the E-value, or the closer it is to zero, the more "significant" the match is. However, keep in mind that virtually identical short alignments have relatively high E values. This is because the calculation of the E value takes into account the length of the query sequence. These high E values make sense because shorter sequences have a higher probability of occurring in the database purely by chance. For more details please see the calculations in the BLAST Course.
#The Expect value can also be used as a convenient way to create a significance threshold for reporting results. You can change the Expect value threshold on most BLAST search pages. When the Expect value is increased from the default value of 10, a larger list with more low-scoring hits can be reported.
#In the context of sequence alignments, a score is a numerical value that describes the overall quality of an alignment. Higher numbers correspond to higher similarity. The score scale depends on the scoring system used (substitution matrix, gap penalty)

#home dir: home/aa/JICAutumn2016/finalAnalysis29Apr/pai_Malvidae
#downloaded orthologs from Phytozome: file:///home/aa/JICAutumn2016/finalAnalysis29Apr/pai_Malvidae/orthologsTop50.txt - query protein sequences (multifasta) of all lyrata meiotic proteins: 
file:///home/aa/JICAutumn2016/finalAnalysis29Apr/pai_Malvidae/meiosisProteinsPeptide_part1a.txt
file:///home/aa/JICAutumn2016/finalAnalysis29Apr/pai_Malvidae/meiosisProteinsPeptide_part1b.txt
file:///home/aa/JICAutumn2016/finalAnalysis29Apr/pai_Malvidae/meiosisProteinsPeptide_part2.txt
Blosum62, target proteome, BLASTP, E treshold -1, 50 best hits, target malvidae
output:
Blast output 10: csv format, Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
file:///home/aa/JICAutumn2016/finalAnalysis29Apr/pai_Malvidae/part1a.bin
file:///home/aa/JICAutumn2016/finalAnalysis29Apr/pai_Malvidae/part1b.bin
file:///home/aa/JICAutumn2016/finalAnalysis29Apr/pai_Malvidae/part2.bin
merged together: 
file:///home/aa/JICAutumn2016/finalAnalysis29Apr/pai_Malvidae/orthologsTop50.txt
ctrlH: , \t
awk '{print $(NF-9)","$(NF-8)","$(NF-7)","$(NF-6)","$(NF-5)","$(NF-4)","$(NF-3)","$(NF-2)","$(NF-1)","$NF}'  orthologsTop50.txt > orthologsTop50numbers.txt
ctrlH: (PAC: @

R - makeOrtologList.r - part1
...
output:
file:///home/aa/JICAutumn2016/finalAnalysis29Apr/pai_Malvidae/orthologIDtoBLAST.txt
file:///home/aa/JICAutumn2016/finalAnalysis29Apr/pai_Malvidae/orthologsbestHitsPerSpecies.txt

Use orthologIDtoBLAST.txt for keyoword search in Phytozome: target Malvidae, no trailing wildcard - 1285  hits out of 1334.??
view (by 1000, to cart, download):
file:///home/aa/JICAutumn2016/finalAnalysis29Apr/pai_Malvidae/bestHits.fa
file:///home/aa/JICAutumn2016/finalAnalysis29Apr/pai_Malvidae/bestHitsDetails.txt

bestHits.fa: ctrlH:
\)\n	\)~~~\n
\n	NA
>	\n
 
R - makeOrtologList.r - part2

Geneious - batch alignment (default), graphs, show identity, coverage, mean PI, mean hydrophobicity
export - csvs


