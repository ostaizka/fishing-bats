#Fishing bats: 2x250 sequencing, problematic for DADA2

#TERMINAL: prepare the sequences for DADA2
#Demultiplexing
mkdir ${workdir}/1-QualityFiltered_DADA2
module load AdapterRemoval/2.2.2

fastqpath=
fastqfile=

cd Mads
AdapterRemoval --file1 ${fastqpath}/${fastqfile}_1.fq.gz \
--file2 ${fastqpath}/${fastqfile}_2.fq.gz \
--basename $WORK/"$a"/ \
--barcode-list $WORK/"$a".txt --barcode-mm-r1 2 --barcode-mm-r2 2


#GM1
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM1.err -o ${workdir}/0-Log/AdapterRemoval_GM1.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM1 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_VFXF_GM1_S1_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_VFXF_GM1_S1_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM1 --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM1_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM1.err -o ${workdir}/0-Log/AdapterRemoval_GM1.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM1 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_3NJ7_Microbiome1_S7_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_3NJ7_Microbiome1_S7_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM1b --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM1_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

#GM2
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM2.err -o ${workdir}/0-Log/AdapterRemoval_GM2.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM2 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_VFXF_GM2_S2_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_VFXF_GM2_S2_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM2 --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM2_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM2.err -o ${workdir}/0-Log/AdapterRemoval_GM2.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM2 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_TJRJ_Microbiome_2_S2_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_TJRJ_Microbiome_2_S2_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM2b --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM2_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

#GM3
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM3.err -o ${workdir}/0-Log/AdapterRemoval_GM3.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM3 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_VFXF_GM3_S3_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_VFXF_GM3_S3_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM3 --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM3_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM3.err -o ${workdir}/0-Log/AdapterRemoval_GM3.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM3 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_TJRJ_Microbiome_3_S3_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_TJRJ_Microbiome_3_S3_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM3b --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM3_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

#GM4
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM4.err -o ${workdir}/0-Log/AdapterRemoval_GM4.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM4 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_VFXF_GM4_S4_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_VFXF_GM4_S4_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM4 --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM4_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM4.err -o ${workdir}/0-Log/AdapterRemoval_GM4.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM4 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_TJRJ_Microbiome_4_S4_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_TJRJ_Microbiome_4_S4_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM4b --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM4_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

#GM5
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM5.err -o ${workdir}/0-Log/AdapterRemoval_GM5.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM5 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_VFXF_GM5_S5_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_VFXF_GM5_S5_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM5 --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM5_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM5.err -o ${workdir}/0-Log/AdapterRemoval_GM5.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM5 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_TJRJ_Microbiome_5_S5_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_TJRJ_Microbiome_5_S5_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM5b --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM5_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

#GM6
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM6.err -o ${workdir}/0-Log/AdapterRemoval_GM6.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM6 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_VFXF_GM6_S6_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_VFXF_GM6_S6_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM6 --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM6_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM6.err -o ${workdir}/0-Log/AdapterRemoval_GM6.out -l nodes=1:ppn=1,mem=20gb,walltime=0:01:00:00 -N AdapterRemoval_GM6 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_TJRJ_Microbiome_6_S6_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_TJRJ_Microbiome_6_S6_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM6b --minquality 30 --trimns --maxns 5 --trimqualities --threads 1 --gzip --barcode-list ${workdir}/0-Documents/tags_GM6_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

#GM7
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM7.err -o ${workdir}/0-Log/AdapterRemoval_GM7.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM7 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_VFXF_GM7_S7_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_VFXF_GM7_S7_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM7 --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM7_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

#GM8
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM8.err -o ${workdir}/0-Log/AdapterRemoval_GM8.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM8 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_VFXF_GM8_S8_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_VFXF_GM8_S8_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM8 --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM8_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

#GM9
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM9.err -o ${workdir}/0-Log/AdapterRemoval_GM9.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM9 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_VFXF_GM9_S9_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_VFXF_GM9_S9_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM9 --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM9_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

#GM10
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM10.err -o ${workdir}/0-Log/AdapterRemoval_GM10.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM10 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_VFXF_GM10_S10_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_VFXF_GM10_S10_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM10 --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM10_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

#GM11
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM11.err -o ${workdir}/0-Log/AdapterRemoval_GM11.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM11 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_VFXF_GM11_S11_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_VFXF_GM11_S11_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM11 --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM11_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

#GM12
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/0-Log/AdapterRemoval_GM12.err -o ${workdir}/0-Log/AdapterRemoval_GM12.out -l nodes=1:ppn=28,mem=8gb,walltime=1:00:00:00 -N AdapterRemoval_GM12 -de AdapterRemoval --file1 ${workdir}/0-Rawdata/TOG_VFXF_GM12_S12_L001_R1_001.fastq.gz --file2 ${workdir}/0-Rawdata/TOG_VFXF_GM12_S12_L001_R2_001.fastq.gz --basename ${workdir}/1-QualityFiltered_DADA2/GM12 --minquality 30 --trimns --maxns 5 --trimqualities --threads 28 --gzip --barcode-list ${workdir}/0-Documents/tags_GM12_AR.txt --barcode-mm-r1 2 --barcode-mm-r2 2

#Remove unnecesary files
rm 1-QualityFiltered_DADA2/*.discarded.gz
rm 1-QualityFiltered_DADA2/*.settings
rm 1-QualityFiltered_DADA2/*.singleton.truncated.gz
rm 1-QualityFiltered_DADA2/*unidentified*
  
#Rename files
samples=$(ls 1-QualityFiltered_DADA2/*.pair1.truncated.gz | sed 's/.pair1.truncated.gz//')
for sample in ${samples}; do
mv ${sample}.pair1.truncated.gz ${sample}.1.fq.gz
mv ${sample}.pair2.truncated.gz ${sample}.2.fq.gz
done

#Merge duplicates
samples=$(ls 1-QualityFiltered_DADA2/*b.*)

for sample in ${samples}; do
sample_tmp=$(echo ${sample} | sed 's/b\./a\./')
sample_final=$(echo ${sample} | sed 's/b\./\./')
mv $sample_final $sample_tmp
cat ${sample} ${sample_tmp} > ${sample_final}
rm $sample
rm $sample_tmp
done

# make a new directory
cd fishingbats
mkdir sequences
#copy the file to the new folder
cp -a  ~/Desktop/raw/. ~/fishingbats/sequences
#1 Get samples name
ls *.1.fq.gz | sed "s/\.1\.fq\.gz//g" > samples

#REMOVE PRIMERS
# Iberian bats: nextera
# Fishing bats: TagSteady
grep "_S" samples > samples_nextera
grep -v "_S" samples > samples_tagsteady
wc -l samples
wc -l samples_nextera
wc -l samples_tagsteady

#samples_nextera
for sample in $(cat samples_nextera)
do
echo "On sample: $sample"
cutadapt -e 0.15 -g ^CCTANGGGNNGCANCAG -G ^GGACTACNNGGGTATCTAAT \
--discard-untrimmed \
-o ${sample}_1_trimmed.fq.gz -p ${sample}_2_trimmed.fq.gz \
${sample}.1.fq.gz ${sample}.2.fq.gz \
>> cutadapt_primer_trimming_stats.txt 2>&1
done

#samples_tagsteady
for sample in $(cat samples_tagsteady)
do
echo "On sample: $sample"
cutadapt -e 0.15 -g ^CTANGGGNNGCANCAG -G ^GACTACNNGGGTATCTAAT \
--discard-untrimmed \
-o ${sample}_1a_trimmed.fq.gz -p ${sample}_2a_trimmed.fq.gz \
${sample}.1.fq.gz ${sample}.2.fq.gz \
>> cutadapt_primer_trimming_stats.txt 2>&1

echo "On sample: $sample"
cutadapt -e 0.15 -g ^GGACTACNNGGGTATCTAAT -G ^CCTANGGGNNGCANCAG \
--discard-untrimmed \
-o ${sample}_1b_trimmed.fq.gz -p ${sample}_2b_trimmed.fq.gz \
${sample}.1.fq.gz ${sample}.2.fq.gz \
>> cutadapt_primer_trimming_stats.txt 2>&1

#Merge both files
cat ${sample}_1a_trimmed.fq.gz ${sample}_2b_trimmed.fq.gz > ${sample}_1_trimmed.fq.gz
cat ${sample}_2a_trimmed.fq.gz ${sample}_1b_trimmed.fq.gz > ${sample}_2_trimmed.fq.gz

#Remove intermediate files
rm ${sample}_1a_trimmed.fq.gz ${sample}_2b_trimmed.fq.gz ${sample}_2a_trimmed.fq.gz ${sample}_1b_trimmed.fq.gz
done

### R

#Before running it, change files path and decide the level of taxonomic assigment:

library(dada2)
packageVersion("dada2") # 1.11.5 when this was put together

setwd("~/fishingbats/Final_copy/")###Change it according to the study

list.files() # make sure what we think is here is actually here

## first we're setting a few variables we're going to use ##
# one with all sample names, by scanning our "samples" file we made earlier
samples <- scan("samples", what="character")

# one holding the file names of all the forward reads
forward_reads <- paste0(samples, "_1_trimmed.fq.gz")
# and one with the reverse
reverse_reads <- paste0(samples, "_2_trimmed.fq.gz")

# and variables holding file names for the forward and reverse
filtered_forward_reads <- paste0(samples, "_1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples, "_2_filtered.fq.gz")

#when some input samples had no read that pass the filter
ls *_1_filtered.fq.gz | sed "s/_1_filtered\.fq\.gz//g" > samples1
samples1 <- scan("samples1", what="character")
filtered_forward_reads <- paste0(samples1, "_1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples1, "_2_filtered.fq.gz")


#Quality trimming/filtering
#plotQualityProfile(forward_reads) #problems to shhow all together
#plotQualityProfile(reverse_reads)
#more dramatic drop in reverse reads
# and just plotting the last 4 samples of the reverse reads
#plotQualityProfile(reverse_reads[18:19])
#plotQualityProfile(forward_reads[1:4])

#Filtering
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=c(223,223), compress=TRUE, multithread=TRUE)

#forward_reads, reverse_reads: input
#filtered_forward_reads, filtered_reverse_reads: output
#maxEE=Maximum amount of estimated errors that you expect to see and read
#rm.phix: removes any reads that match the PhiX bacteriophage genome, which is typically added to Illumina sequencing runs for quality monitoring
#truncLen: parameter setting the minimum size to trim the forward and reverse reads to in order to keep the quality scores roughly above 30 overall
#minLen: is setting the minimum length reads we want to keep after trimming
#no minLen value in dada2 tutorial
#maxN: additional filtering default parameter that is removing any sequences containing any Ns
#truncq= truncates your reads at the first base pair with a quality score equal to or less than 2
#compress=TRUE: FASTQ files gzipped
#multithread=TRUE: FASTQ files process in parallel

#it also generated an object in R called filtered_out. And that’s a matrix holding how many reads went in and how many reads made it out

filtered_out
saveRDS(filtered_out, "~/fishingbats/many_species.RData")

#look at the filtered reads
#plotQualityProfile(filtered_forward_reads)
#plotQualityProfile(filtered_reverse_reads)
#plotQualityProfile(filtered_reverse_reads[17:20])
#plotQualityProfile(filtered_forward_reads[1:2])

#Generating an error model of our data by learning the specific error-signature of our dataset
err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE)
saveRDS(err_forward_reads, "~/fishingbats/many_species/err_forward_reads.RData")
saveRDS(err_reverse_reads, "~/fishingbats/many_species/err_reverse_reads.RData")

#The developers have incorporated a plotting function to visualize how well the estimated error rates match up with the observed:
#plotErrors(err_forward_reads, nominalQ=TRUE)
#plotErrors(err_reverse_reads, nominalQ=TRUE)
#The red line is what is expected based on the quality score, the black line represents the estimate, and the black dots represent the observed.
#generally speaking, you want the observed (black dots) to track well with the estimated (black line)
#In geneal, error rate decreases as Q scores increases

#Dereplication:keep/process one, and just attach the identical sequences to it
#When DADA2 dereplicates sequences, it also generates a new quality-score profile of each unique sequence
#based on the average quality scores of each base of all of the sequences that were replicates of it.
#the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples 
saveRDS(derep_forward, "~/fishingbats/many_species/derep_forward.RData")
saveRDS(derep_reverse, "~/fishingbats/many_species/derep_reverse.RData")

#Inferring ASVs
#It does this by incorporating the consensus quality profiles and abundances of each unique sequence, and then figuring out if each sequence more likely to be of biological origin or more likely to be spurious.
dada_forward <- dada(derep_forward, err=err_forward_reads, multithread=TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, multithread=TRUE)
saveRDS(dada_forward, "~/fishingbats/many_species/dada_forward.RData")
saveRDS(dada_reverse, "~/fishingbats/many_species/dada_reverse.RData")

#Merging forward and reverse reads
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                               derep_reverse, minOverlap=5, verbose=TRUE)
saveRDS(merged_amplicons, "~/fishingbats/many_species/merged_amplicons.RData")
#maxMismatch is by default 0 but you can make it say 1 or 2 if you’re finding that a lot of your forward and reverse reads are not merging
#minOverlap option has been change because thhey were sequenced at 250 bp

#Generating a count table
seqtab <- makeSequenceTable(merged_amplicons)
table(nchar(getSequences(seqtab)))
#to see the breakdown of the size of these amplicons. On the top row is the size of the merged reads, and on the botton is the frequency

#Chimera identification
seqtab.nochim <- removeBimeraDenovo(seqtab, multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, "~/fishingbats/Final_copy/seqtab.nochim.RData")

# this is one quick way to look at sequences that have been lost, to know whether they held a lot in terms of abundance
sum(seqtab.nochim)/sum(seqtab)

#Overview of counts throughout:quick way to pull out how many reads were dropped at various points of the pipeline
# set a little function
getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

summary_tab

write.table(summary_tab, "summary_reads_table.txt", sep="\t", quote=F)
#Assigning taxonomy(before this step, download the database)
#There are different DADA2-formatted databases available in DADA2 website
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", tryRC=T)
saveRDS(taxa, "taxa.RData")

#Extracting the standard goods from DADA2
#giving to seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

#making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

#count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

#tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)

#save the working environment
save.image("all_species.RData")
