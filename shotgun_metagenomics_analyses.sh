workdir="/home/projects/ku-cbd/people/antalb/Fishing_bats_GM_shotgun"

####
# Transfer data
####

mkdir ${workdir}/0-Rawdata
cd ${workdir}/0-Rawdata
erda
cd rawdata/20160831_modern_guano_metagenomics/
get Sample_TOG_MXDA_ES166/*R1_00[1-9].fastq
get Sample_TOG_MXDA_ES172/*R1_00[1-9].fastq
get Sample_TOG_MXDA_ES174/*R1_00[1-9].fastq
get Sample_TOG_MXDA_ES244/*R1_00[1-9].fastq
get Sample_TOG_MXDA_ES354/*R1_00[1-9].fastq
get Sample_TOG_MXDA_ES365/*R1_00[1-9].fastq
get Sample_TOG_MXDA_ES437/*R1_00[1-9].fastq
get Sample_TOG_MXDA_ES441/*R1_00[1-9].fastq
get Sample_TOG_MXDA_ES466/*R1_00[1-9].fastq
get Sample_TOG_MXDA_ES468/*R1_00[1-9].fastq
get Sample_TOG_MXDA_ES470/*R1_00[1-9].fastq
get Sample_TOG_MXDA_ES472/*R1_00[1-9].fastq
get Sample_TOG_MXDA_NO101/*R1_00[1-9].fastq
get Sample_TOG_MXDA_NO102/*R1_00[1-9].fastq
get Sample_TOG_MXDA_ES264/*R1_00[1-9].fastq

cd ../20170329_modern_guano_metagenomics/
get Sample_TOG_MMDC_NO104/*R1_00[1-9].fastq
get Sample_TOG_MMDC_NO106/*R1_00[1-9].fastq
get Sample_TOG_MMDC_NO107/*R1_00[1-9].fastq
get Sample_TOG_MMDC_ES190/*R1_00[1-9].fastq
get Sample_TOG_MMDC_ES610/*R1_00[1-9].fastq
get Sample_TOG_MMDC_ES612/*R1_00[1-9].fastq
get Sample_TOG_MMDC_ES619/*R1_00[1-9].fastq
get Sample_TOG_MMDC_ES621/*R1_00[1-9].fastq

####
# Concatenate samples
####

cd ${workdir}/0-Rawdata
sample=NO102
cat *${sample}*_R1* > ${sample}.fastq
rm TOG_*

####
# Quality-filtering
####

mkdir ${workdir}/1-QualityFilter
qsub -V -A ku-cbd -W group_list=ku-cbd -v "workdir=${workdir}"  -d `pwd` -e ${workdir}/QualityFiltering.err -o ${workdir}/QualityFiltering.out -l nodes=1:ppn=40,mem=50gb,walltime=1:00:00:00 -N QualityFiltering ${workdir}/qualityfiltering.sh

#qualityfiltering.sh
module load tools gcc adapterremoval/2.2.4
while read SAMPLE;do
SAMPLE2=$(echo $SAMPLE | sed 's/\.fastq//')
AdapterRemoval --file1 ${workdir}/0-Rawdata/${SAMPLE} --basename ${workdir}/1-QualityFilter/${SAMPLE2} --minquality 30 --minlength 80 --trimqualities --trimns --maxns 5 --threads 40
done < <(ls ${workdir}/0-Rawdata)

#Modify names
cd ${workdir}/1-QualityFilter/
while read SAMPLE;do
SAMPLE2=$(echo $SAMPLE | sed 's/truncated/fastq/')
mv ${SAMPLE} ${SAMPLE2}
done < <(ls *truncated)
rm *discarded
rm *settings
cd ${workdir}

####
# Map to host
####

mkdir ${workdir}/2-Maptohost
sample=NO107
qsub -V -A ku-cbd -W group_list=ku-cbd -v "workdir=${workdir},sample=${sample}"  -d `pwd` -e ${workdir}/MapToHost.err -o ${workdir}/MapToHost.out -l nodes=1:ppn=40,mem=50gb,walltime=1:00:00:00 -N MapToHost ${workdir}/maptohost.sh

#maptohost.sh
module load tools samtools/1.10 bwa/0.7.15
bwa mem -t 40 -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" /home/projects/ku-cbd/people/antalb/reference_genomes/allChiroptera.fasta ${workdir}/1-QualityFilter/${sample}.fastq | samtools view -f4 -T ${workdir}/2-Maptohost/${sample} -b - | samtools sort -T ${workdir}/2-Maptohost/${sample} - | samtools fastq - > ${workdir}/2-Maptohost/${sample}.fastq
#

####
# Map to fish
####

mkdir ${workdir}/3-Maptofish
sample=NO107
qsub -V -A ku-cbd -W group_list=ku-cbd -v "workdir=${workdir},sample=${sample}"  -d `pwd` -e ${workdir}/MapToFish.err -o ${workdir}/MapToFish.out -l nodes=1:ppn=40,mem=50gb,walltime=1:00:00:00 -N MapToFish ${workdir}/maptofish.sh


#maptofish.sh
module load tools samtools/1.10 bwa/0.7.15
bwa mem -t 40 -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" /home/projects/ku-cbd/people/antalb/reference_genomes/allPoeciliidae.fasta ${workdir}/2-Maptohost/${sample}.fastq | samtools view -F4 -T ${workdir}/3-Maptofish/${sample} -b - | samtools sort -T ${workdir}/3-Maptofish/${sample} - | samtools fastq - > ${workdir}/3-Maptofish/${sample}.fastq
#

####
# Get fish statistics
####

sample=NO107
host=$(grep -c "^@" ${workdir}/2-Maptohost/${sample}.fastq)
fish=$(grep -c "^@" ${workdir}/3-Maptofish/${sample}.fastq)
((fishperc = $fish * 1000000 / $host))
echo $fishperc

####
# Map to fish (metagenomic)
####

mkdir ${workdir}/4-Metagenomic
sample=NO107
qsub -V -A ku-cbd -W group_list=ku-cbd -v "workdir=${workdir},sample=${sample}"  -d `pwd` -e ${workdir}/Metagenomic.err -o ${workdir}/Metagenomic.out -l nodes=1:ppn=40,mem=50gb,walltime=1:00:00:00 -N Metagenomic ${workdir}/metagenomic.sh

#metagenomic.sh
module load tools samtools/1.10 bwa/0.7.15
bwa mem -t 40 -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" /home/projects/ku-cbd/people/antalb/reference_genomes/allPoeciliidae.fasta ${workdir}/2-Maptohost/${sample}.fastq | samtools view -f4 -T ${workdir}/4-Metagenomic/${sample} -b - | samtools sort -T ${workdir}/4-Metagenomic/${sample} - | samtools fastq - > ${workdir}/4-Metagenomic/${sample}.fastq
#

####
# Metagenomic assembly
####

mkdir ${workdir}/5-Assembly
cat ${workdir}/4-Metagenomic/*fastq > ${workdir}/5-Assembly/allsamples.fastq

qsub -V -A ku-cbd -W group_list=ku-cbd -v "workdir=${workdir}"  -d `pwd` -e ${workdir}/Assembly.err -o ${workdir}/Assembly.out -l nodes=1:ppn=40,mem=50gb,walltime=1:00:00:00 -N Assembly ${workdir}/assembly.sh

#assembly.sh
module load tools megahit/1.1.1
megahit -r ${workdir}/5-Assembly/allsamples.fastq -o ${workdir}/5-Assembly/megahit
#

#Move and rename assembly
cp ${workdir}/5-Assembly/megahit/final.contigs.fa ${workdir}/5-Assembly/assembly.fa

####
# Map to assembly
####

mkdir ${workdir}/6-MapToAssembly
qsub -V -A ku-cbd -W group_list=ku-cbd -v "workdir=${workdir}"  -d `pwd` -e ${workdir}/MapToAssembly.err -o ${workdir}/MapToAssembly.out -l nodes=1:ppn=40,mem=50gb,walltime=2:00:00:00 -N MapToAssembly ${workdir}/maptoassembly.sh

#maptoassembly.sh
module load tools samtools/1.10 bwa/0.7.15

#Index assembly
samtools faidx ${workdir}/5-Assembly/assembly.fa
bwa index ${workdir}/5-Assembly/assembly.fa

#Run mappings
cd ${workdir}/4-Metagenomic
while read SAMPLE;do
SAMPLE2=$(echo $SAMPLE | sed 's/\.fastq//')
bwa mem -t 40 -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" ${workdir}/5-Assembly/assembly.fa ${workdir}/4-Metagenomic/${SAMPLE} | samtools view -F4 -T ${workdir}/6-MapToAssembly/${SAMPLE2} -b - | samtools sort -T ${workdir}/6-MapToAssembly/${SAMPLE2} - > ${workdir}/6-MapToAssembly/${SAMPLE2}.bam
done < <(ls *.fastq)

#Create depth table
module load tools perl/5.20.2 metabat/2.12.1 maxbin/2.2.7 fraggenescan/1.31
jgi_summarize_bam_contig_depths --outputDepth ${workdir}/contig.depth.txt ${workdir}/6-MapToAssembly/*.bam
#

####
# Prodigal
####

mkdir ${workdir}/7-Prodigal
#Split assembly
lines=$(cat ${workdir}/5-Assembly/assembly.fa | wc -l)
chunksize=$(($lines / 40))
split -l ${chunksize} ${workdir}/5-Assembly/assembly.fa ${workdir}/7-Prodigal/assembly.fa_


qsub -V -A ku-cbd -W group_list=ku-cbd -v "workdir=${workdir}"  -d `pwd` -e ${workdir}/Prodigal.err -o ${workdir}/Prodigal.out -l nodes=1:ppn=40,mem=50gb,walltime=2:00:00:00 -N Prodigal ${workdir}/prodigal.sh

#prodigal.sh

#Run prodigal in parallel
module load tools prodigal/2.6.3 parallel/20170822

cd ${workdir}/7-Prodigal
function prodigal_parallel() {
chunck=${1}
workdir=${2}
prodigal -p meta -i ${chunck} -o ${chunck}.gbk -a ${chunck}.faa -d ${chunck}.fna
}

export -f prodigal_parallel
parallel -j 40 -k prodigal_parallel {} ${workdir} < <(ls ${workdir}/7-Prodigal/assembly.fa_*)
#
cat ${workdir}/7-Prodigal/assembly.fa_*.faa > ${workdir}/7-Prodigal/assembly.orf.faa
cat ${workdir}/7-Prodigal/assembly.fa_*.fna > ${workdir}/7-Prodigal/assembly.orf.fna
cat ${workdir}/7-Prodigal/assembly.fa_*.gbk > ${workdir}/7-Prodigal/assembly.orf.gbk
rm ${workdir}/7-Prodigal/assembly.fa_*

####
# Diamond database
####

mkdir 8-Enterotoxins
qsub -V -A ku-cbd -W group_list=ku-cbd -v "workdir=${workdir}"  -d `pwd` -e ${workdir}/Diamond.err -o ${workdir}/Diamond.out -l nodes=1:ppn=40,mem=50gb,walltime=2:00:00:00 -N Diamond ${workdir}/diamond.sh

# diamond.sh
module load tools diamond/0.9.29
diamond makedb --in ${workdir}/7-Prodigal/assembly.orf.faa -d ${workdir}/7-Prodigal/assembly.orf
diamond blastp -d ${workdir}/7-Prodigal/assembly.orf -q ${workdir}/8-Enterotoxins/aeromonas_enterotoxins.faa -o ${workdir}/8-Enterotoxins/matches.tsv --outfmt 6
#

sort -k3,3n ${workdir}/8-Enterotoxins/matches.tsv

#Top Matched refs
#tr|B5AGV9|B5AGV9_AERJA - Aerolysin (https://www.ebi.ac.uk/interpro/entry/InterPro/IPR005830/)
grep -A5 "tr|B5AGV9|B5AGV9_AERJA" ${workdir}/8-Enterotoxins/aeromonas_enterotoxins.faa
grep -A5 "tr|A0A125R5P4|A0A125R5P4_9GAMM" ${workdir}/8-Enterotoxins/aeromonas_enterotoxins.faa
grep -A5 "tr|J9Q028|J9Q028_AERHY" ${workdir}/8-Enterotoxins/aeromonas_enterotoxins.faa

#Matched genes: k79_89359_1, k79_106877_1
grep -A5 "k79_89359_1" ${workdir}/7-Prodigal/assembly.orf.faa
grep -A5 "k79_106877_1" ${workdir}/7-Prodigal/assembly.orf.faa

grep "k79_89359" ${workdir}/contig.depth.txt
grep "k79_106877" ${workdir}/contig.depth.txt

mkdir 8-Enterotoxins
qsub -V -A ku-cbd -W group_list=ku-cbd -v "workdir=${workdir}"  -d `pwd` -e ${workdir}/Diamond2.err -o ${workdir}/Diamond2.out -l nodes=1:ppn=40,mem=50gb,walltime=2:00:00:00 -N Diamond2 ${workdir}/diamond2.sh

# diamond2.sh
module load tools diamond/0.9.29
diamond blastp -d ${workdir}/7-Prodigal/assembly.orf -q ${workdir}/8-Enterotoxins/aeromonas.faa -o ${workdir}/8-Enterotoxins/matches_aeromonas.tsv --outfmt 6
#

#Identify 100% Aeromonas contigs
awk '$3 == "100.0" {print $2}' 8-Enterotoxins/matches_aeromonas.tsv | sed 's/_[1-9]$//' | sed 's/_[1-9][1-9]$//' | uniq > ${workdir}/8-Enterotoxins/aeromonas_contigs.txt

grep -w -f ${workdir}/8-Enterotoxins/aeromonas_contigs.txt ${workdir}/contig.depth.txt > ${workdir}/8-Enterotoxins/aeromonas.depth.txt

#T.test
#mean aeromonas contig coverage
aeromonas <- c(0.00137075,0.001193034,0.006156701,0.004882668,0.045648883,0.09518639,0.84134901,1.066127024,0.014080819,0.012956483,0.007286074,0.006719403,0.009394763,0.011356108,0.043379981,0.08597787,0.741331114,1.290935477,12.08785807,17.32587333,53.25491051,82.10130173,55.50530348,87.30866046,2.696369453,8.33054786,82.55743924,140.929274,17.25476361,18.92231086,12.70370662,13.36130202,19.62067046,24.26040817,0.023618183,0.016171559,0.043674228,0.184545785,0.401036216,0.040991619,0.126643788,0.745301308,0.001832268,0.000745508,0.007199385,0.011062036)
#mean ash3 coverage
ash3 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.714286,1.19888,16.3896,17.4289,58.3961,113.169,53.1364,79.1382,2.54545,0.837789,100.786,78.4701,19.7987,58.2403,12.0455,6.7757,21.5844,5.1595,0,0,0,0,0,0,0,0,0,0,0,0)
#t.test
t.test(aeromonas,ash3,paired=TRUE)
pdf("/Users/anttonalberdi/Google\ Drive/Projects/2020_Fishing_bats_mb/figures/aeromona_ash.pdf",width=7,height=6)
reg1 <- lm(aeromonas~ash3)
plot(aeromonas,ash3)
abline(reg1)
dev.off()



####
# MetaErg
####

metaergdir="/home/projects/ku-cbd/people/antalb/software/metaerg"

cd ${workdir}
mkdir ${workdir}/7-Metaerg
qsub -V -A ku-cbd -W group_list=ku-cbd -v "workdir=${workdir},metaergdir=${metaergdir}"  -d `pwd` -e ${workdir}/Metaerg.err -o ${workdir}/Metaerg.out -l nodes=1:ppn=40,mem=150gb,walltime=3:00:00:00 -N Metaerg ${workdir}/metaerg.sh

#metaerg.sh
source /home/people/antalb/.bashrc
module load tools anaconda2/4.4.0 perl/5.30.2 aragorn/1.2.36 ncbi-blast/2.10.0+ diamond/0.9.29 hmmer/3.2.1 java/1.8.0 minced/0.2.0 minpath/1.4 prodigal/2.6.3 signalp/4.1c tmhmm/2.0c
perl ${metaergdir}/bin/metaerg.pl --cpus 40 --force --outdir ${workdir}/7-Metaerg --depth ${workdir}/contig.depth.txt ${workdir}/5-Assembly/assembly.fa
#

# error in generating all.gff

module load tools perl/5.30.2
mkdir ${workdir}/7-Metaerg2
perl ${metaergdir}/bin/output_reports.pl -g ${workdir}/7-Metaerg/tmp/features.annot.gff -f ${workdir}/7-Metaerg/metaerg.pl_05182020.fna -o ${workdir}/7-Metaerg2
####
# MetaErg local docker
####

docker run --shm-size 2g --rm -u $(id -u):$(id -g) -it -v my_local_dir:/data/ xiaolidong/docker-metaerg metaerg.pl -h

####
# MetaErg singularity
####

module load tools anaconda3/4.4.0 singularity/3.5.3
source activate anvio6.1

singularity shell docker://xiaolidong/docker-metaerg
singularity run docker://xiaolidong/docker-metaerg
singularity exec docker://ubuntu:latest echo "Hello Dinosaur!"

singularity pull docker://ubuntu:latest
singularity build ubuntu.img docker://ubuntu:latest

####
# KO annotation using GhostKoala
####

####
# Map to cds
####

mkdir ${workdir}/8-MapToCDS
qsub -V -A ku-cbd -W group_list=ku-cbd -v "workdir=${workdir}"  -d `pwd` -e ${workdir}/MapToCDS.err -o ${workdir}/MapToCDS.out -l nodes=1:ppn=40,mem=50gb,walltime=2:00:00:00 -N MapToCDS ${workdir}/maptocds.sh

#maptocds.sh
module load tools samtools/1.10 bwa/0.7.15

#Index assembly
cp ${workdir}/7-Prodigal/assembly.orf.fna ${workdir}/8-MapToCDS/assembly.orf.fna
samtools faidx ${workdir}/8-MapToCDS/assembly.orf.fna
bwa index ${workdir}/8-MapToCDS/assembly.orf.fna

#Run mappings
cd ${workdir}/4-Metagenomic
while read SAMPLE;do
SAMPLE2=$(echo $SAMPLE | sed 's/\.fastq//')
bwa mem -t 40 -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" ${workdir}/8-MapToCDS/assembly.orf.fna ${workdir}/4-Metagenomic/${SAMPLE} | samtools view -F4 -T ${workdir}/8-MapToCDS/${SAMPLE2} -b - | samtools sort -T ${workdir}/8-MapToCDS/${SAMPLE2} - > ${workdir}/8-MapToCDS/${SAMPLE2}.bam
done < <(ls *.fastq)

#Create depth table
module load tools perl/5.20.2 metabat/2.12.1 maxbin/2.2.7 fraggenescan/1.31
jgi_summarize_bam_contig_depths --outputDepth ${workdir}/cds.depth.txt ${workdir}/8-MapToCDS/*.bam
#


####
# KEGG analyses
####

depth <- c(4918364,4331755,3007944,5892973,4901059,7279808,405643,2864710,7087328,4172882,5970708,6307895,1270124,5110039,3766320,5074952,3241742,440725)
#with noctilio
#depth <- c(4918364,4331755,3007944,5892973,4901059,7279808,405643,2864710,7087328,4172882,5970708,6307895,1270124,5110039,3766320,5074952,3241742,440725,25255636,18676600,25132296,43833236,26637564)
norm_factor <- mean(depth) / depth

setwd("/Users/anttonalberdi/Google\ Drive/Projects/2020_Fishing_bats_mb/shotgun/")

cov_table <- read.table("cds.depth.noNL.txt",sep="\t",header=TRUE)
cov_table <- cov_table[,c(1,2,3,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38)]
#With Noctilio
#cov_table <- cov_table[,c(1,2,3,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48)]

gene_ko_table <- read.csv("user_ko.noNL.txt",sep="\t",header=FALSE,na.strings="NA")
colnames(gene_ko_table) <- c("contigName","ko")
cov_table_ko <- merge(cov_table,gene_ko_table,by="contigName")

library(data.table)
kegg_meta <- data.frame(fread("KEGG.metabolism.csv",header=TRUE,sep=",",colClasses="character"))
ko_kegg <- data.frame(fread("ko_path.txt",header=TRUE,sep="\t",colClasses="character"))
ko_kegg_meta <- merge(kegg_meta,ko_kegg,by.x="Pathway",by.y="path")



######
# KO analysis
######

sampleinfo <- read.csv("samples.csv",sep=",")

cov_table_ko <- aggregate(cov_table_ko[,c(4:21)],by=list(cov_table_ko$ko),FUN=sum)
colnames(cov_table_ko) <- gsub(".bam","",colnames(cov_table_ko))
rownames(cov_table_ko) <- cov_table_ko[,1]
cov_table_ko <- cov_table_ko[,-1]
cov_table_ko_norm <- sweep(cov_table_ko, 2, norm_factor, FUN="*")
cov_table_ko_norm_groups <- merge(t(cov_table_ko_norm),sampleinfo,by.x="row.names",by.y="Samples")

#Chitinase (https://www.kegg.jp/dbget-bin/www_bget?ko:K01183)
ko="K01183"
cov_table_ko_norm_groups_sub <- cov_table_ko_norm_groups[(cov_table_ko_norm_groups$Group == "Insectivorous") | (cov_table_ko_norm_groups$Group == "Fishing_capaccinii"),c(ko,"diet")]
colnames(cov_table_ko_norm_groups_sub)[1] <- "value"
t.test(value ~ diet, data = cov_table_ko_norm_groups_sub)

#B12 (https://www.genome.jp/kegg-bin/show_module?M00122)
#not found:
ko=c("K00798","K19221","K02232","K02227","K00768","K02226","K02233","K02231","K02225","K22316")
cov_table_ko_norm_groups_sub <- cov_table_ko_norm_groups[(cov_table_ko_norm_groups$Group == "Insectivorous") | (cov_table_ko_norm_groups$Group == "Fishing_capaccinii"),]
cov_table_ko_norm_groups_sub <- cov_table_ko_norm_groups[(cov_table_ko_norm_groups$Group == "Insectivorous") | (cov_table_ko_norm_groups$Group == "Fishing_capaccinii") | (cov_table_ko_norm_groups$Group == "Myotis_capaccinii"),]
ko_in=ko[ko %in% colnames(cov_table_ko_norm_groups_sub)]
ko_out=ko[!ko %in% colnames(cov_table_ko_norm_groups_sub)]
cov_table_ko_norm_groups_sub2 <- cov_table_ko_norm_groups_sub[,c(ko_in,"diet")]
cov_table_ko_norm_groups_sub2$value <- rowSums(cov_table_ko_norm_groups_sub2[,ko_in])
t.test(value ~ diet, data = cov_table_ko_norm_groups_sub2)

#B2 (https://www.kegg.jp/dbget-bin/www_bget?md:M00125)
ko=c("K01497","K14652","K01498","K00082","K11752","K22912","K20860","K20861","K20862","K21063","K21064","K02858","K14652","K00794","K00793","K00861","K20884","K00953","K22949","K11753")
cov_table_ko_norm_groups_sub <- cov_table_ko_norm_groups[(cov_table_ko_norm_groups$Group == "Insectivorous") | (cov_table_ko_norm_groups$Group == "Fishing_capaccinii"),]
cov_table_ko_norm_groups_sub <- cov_table_ko_norm_groups[(cov_table_ko_norm_groups$Group == "Insectivorous") | (cov_table_ko_norm_groups$Group == "Fishing_capaccinii") | (cov_table_ko_norm_groups$Group == "Myotis_capaccinii"),]
ko_in=ko[ko %in% colnames(cov_table_ko_norm_groups_sub)]
ko_out=ko[!ko %in% colnames(cov_table_ko_norm_groups_sub)]
cov_table_ko_norm_groups_sub2 <- cov_table_ko_norm_groups_sub[,c(ko_in,"diet")]
cov_table_ko_norm_groups_sub2$value <- rowSums(cov_table_ko_norm_groups_sub2[,ko_in])
t.test(value ~ diet, data = cov_table_ko_norm_groups_sub2)

#B6 (https://www.kegg.jp/kegg-bin/show_module?M00126)
ko=c("K03472","K03473","K00831","K00097","K03474","K00275","K23998")
cov_table_ko_norm_groups_sub <- cov_table_ko_norm_groups[(cov_table_ko_norm_groups$Group == "Insectivorous") | (cov_table_ko_norm_groups$Group == "Fishing_capaccinii"),]
cov_table_ko_norm_groups_sub <- cov_table_ko_norm_groups[(cov_table_ko_norm_groups$Group == "Insectivorous") | (cov_table_ko_norm_groups$Group == "Fishing_capaccinii") | (cov_table_ko_norm_groups$Group == "Myotis_capaccinii"),]
ko_in=ko[ko %in% colnames(cov_table_ko_norm_groups_sub)]
ko_out=ko[!ko %in% colnames(cov_table_ko_norm_groups_sub)]
cov_table_ko_norm_groups_sub2 <- cov_table_ko_norm_groups_sub[,c(ko_in,"diet")]
cov_table_ko_norm_groups_sub2$value <- rowSums(cov_table_ko_norm_groups_sub2[,ko_in])
t.test(value ~ diet, data = cov_table_ko_norm_groups_sub2)


#B9 (https://www.kegg.jp/kegg-bin/show_module?M00126)
ko=c("K01495","K09007","K22391","K01077","K01113","K08310","K19965","K13939","K13940","K01633","K00950","K00796","K01633","K13941","K11754","K20457","K00287","K13998")
cov_table_ko_norm_groups_sub <- cov_table_ko_norm_groups[(cov_table_ko_norm_groups$Group == "Insectivorous") | (cov_table_ko_norm_groups$Group == "Fishing_capaccinii"),]
cov_table_ko_norm_groups_sub <- cov_table_ko_norm_groups[(cov_table_ko_norm_groups$Group == "Insectivorous") | (cov_table_ko_norm_groups$Group == "Fishing_capaccinii") | (cov_table_ko_norm_groups$Group == "Myotis_capaccinii"),]
ko_in=ko[ko %in% colnames(cov_table_ko_norm_groups_sub)]
ko_out=ko[!ko %in% colnames(cov_table_ko_norm_groups_sub)]
cov_table_ko_norm_groups_sub2 <- cov_table_ko_norm_groups_sub[,c(ko_in,"diet")]
cov_table_ko_norm_groups_sub2$value <- rowSums(cov_table_ko_norm_groups_sub2[,ko_in])
t.test(value ~ diet, data = cov_table_ko_norm_groups_sub2)

######
# B12 contigs
######
ko=c("K00798","K19221","K02232","K02227","K00768","K02226","K02233","K02231","K02225","K22316")
b12_contigs <- gene_ko_table[gene_ko_table$ko %in% ko,]

b12_cov_table <- cov_table[cov_table$contigName %in% b12_contigs[,1],]
rownames(b12_cov_table) <- b12_cov_table[,1]
b12_cov_table <- b12_cov_table[,-c(1:3)]
sort(rowSums(b12_cov_table))

perl ../scripts/FastaToOnliner.pl 8-MapToCDS_noNL/cds.fna > 8-MapToCDS_noNL/cds.oneliner.fna

grep -A1 "k79_96446_cds_1" 8-MapToCDS_noNL/cds.oneliner.fna



######
# Pathway analysis
######


cov_table_kegg <- merge(cov_table_ko,ko_kegg_meta,by.x="row.names", by.y="ko")
cov_table_pathway <- aggregate(cov_table_kegg[,c(2:19)],by=list(cov_table_kegg$Pathway),FUN=sum)

rownames(cov_table_pathway) <- cov_table_pathway[,1]
cov_table_pathway <- cov_table_pathway[,-1]
cov_table_pathway_norm <- sweep(cov_table_pathway, 2, norm_factor, FUN="*")

sampleinfo <- read.csv("samples.csv")

cov_table_pathway_norm_groups <- merge(t(cov_table_pathway_norm),sampleinfo,by.x="row.names",by.y="Samples")
cov_table_pathway_norm_groups$diet <- as.character(cov_table_pathway_norm_groups$diet)

pathways <- colnames(cov_table_pathway_norm_groups)[2:134]

table <- c()
for (p in pathways){
pathway=p
cov_table_pathway_norm_groups_sub <- cov_table_pathway_norm_groups[(cov_table_pathway_norm_groups$diet == "Insectivorous") | (cov_table_pathway_norm_groups$diet == "Fishing"),c(pathway,"diet")]
colnames(cov_table_pathway_norm_groups_sub)[1] <- "value"
row <- c(pathway,t.test(value ~ diet, data = cov_table_pathway_norm_groups_sub)$estimate[1],t.test(value ~ diet, data = cov_table_pathway_norm_groups_sub)$estimate[2],t.test(value ~ diet, data = cov_table_pathway_norm_groups_sub)$p.value)
table <- rbind(table,row)
}
colnames(table) <- c("Pathway","Fishing","Insectivorous","pvalue")
table <- as.data.frame(table)
table[,1] <- as.character(table[,1])
table[,2] <- as.numeric(as.character((table[,2])))
table[,3] <- as.numeric(as.character((table[,3])))
table[,4] <- as.numeric(as.character((table[,4])))

significant <- table[table[,4] < 0.05,]
significant <- merge(significant,kegg_meta,by="Pathway")
significant

pathway="005112"
cov_table_pathway_norm_groups_sub <- cov_table_pathway_norm_groups[(cov_table_pathway_norm_groups$diet == "Insectivorous") || (cov_table_pathway_norm_groups$diet == "Fishing"),c(pathway,"diet")]
colnames(cov_table_pathway_norm_groups_sub)[1] <- "value"
t.test(value ~ diet, data = cov_table_pathway_norm_groups_sub)

####
# Heatmap
####

cov_table_pathway_norm_groups
rownames(cov_table_pathway_norm_groups) <- cov_table_pathway_norm_groups[,1]
cov_table_pathway_norm_groups_significant <- cov_table_pathway_norm_groups[,significant[,1]]

cov_table_pathway_norm_groups_significant2 <- sweep(cov_table_pathway_norm_groups_significant, 2, colMeans(cov_table_pathway_norm_groups_significant), FUN="/")
library(gplots)
heatmap.2(as.matrix(cov_table_pathway_norm_groups_significant2))

######
# KO analysis (without non-fishing capaccinii)
######
