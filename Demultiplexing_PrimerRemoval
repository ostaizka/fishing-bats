#### Demultiplexing ####

mkdir ${workdir}/1-QualityFiltered_DADA2
module load AdapterRemoval/2.2.2

fastqpath=
fastqfile=

cd Mads
AdapterRemoval --file1 ${fastqpath}/${fastqfile}_1.fq.gz \
--file2 ${fastqpath}/${fastqfile}_2.fq.gz \
--basename $WORK/"$a"/ \
--barcode-list $WORK/"$a".txt --barcode-mm-r1 2 --barcode-mm-r2 2

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

#### PRIMER REMOVAL ####
#Primers:341F/806R
for sample in $(cat samples)
do
echo "On sample: $sample"
cutadapt -e 0.15 -g ^CTANGGGNNGCANCAG -G ^GACTACNNGGGTATCTAAT \
--discard-untrimmed \
-o ${sample}_1a_trimmed.fq.gz -p ${sample}_2a_trimmed.fq.gz \
${sample}.1.fq.gz ${sample}.2.fq.gz \
>> cutadapt_primer_trimming_statsa.txt 2>&1

echo "On sample: $sample"
cutadapt -e 0.15 -g ^GACTACNNGGGTATCTAAT -G ^CTANGGGNNGCANCAG \
--discard-untrimmed \
-o ${sample}_1b_trimmed.fq.gz -p ${sample}_2b_trimmed.fq.gz \
${sample}.1.fq.gz ${sample}.2.fq.gz \
>> cutadapt_primer_trimming_statsb.txt 2>&1

#Merge both files
cat ${sample}_1a_trimmed.fq.gz ${sample}_2b_trimmed.fq.gz > ${sample}_1_trimmed.fq.gz
cat ${sample}_2a_trimmed.fq.gz ${sample}_1b_trimmed.fq.gz > ${sample}_2_trimmed.fq.gz

#Remove intermediate files
rm ${sample}_1a_trimmed.fq.gz ${sample}_2b_trimmed.fq.gz ${sample}_2a_trimmed.fq.gz ${sample}_1b_trimmed.fq.gz
done
