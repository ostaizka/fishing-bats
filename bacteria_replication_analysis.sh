workdir="/home/projects/ku-cbd/people/antalb/aeromonas_irep"

#####
# Index Aeromonas salmonicida genome
#####

#jobs/index_referece_Asa.sh
module load tools samtools/1.11 bowtie2/2.4.2 pigz/2.3.4
bowtie2-build ${workdir}/reference/GCF_012931585.1_ASM1293158v1_genomic.fna ${workdir}/reference/aeromonas_salmonicida

qsub -V -A ku-cbd -W group_list=ku-cbd -v "workdir=${workdir}"  -d `pwd` -e ${workdir}/jobs/index_referece.err -o ${workdir}/jobs/index_referece.out -l nodes=1:ppn=40,mem=180gb,walltime=1:00:00:00 -N index_referece ${workdir}/jobs/index_referece_Asa.sh

#jobs/index_referece_Ave.sh
module load tools samtools/1.11 bowtie2/2.4.2 pigz/2.3.4
bowtie2-build ${workdir}/reference/GCF_008693705.1_ASM869370v1_genomic.fna ${workdir}/reference/aeromonas_veronii

qsub -V -A ku-cbd -W group_list=ku-cbd -v "workdir=${workdir}"  -d `pwd` -e ${workdir}/jobs/index_referece.err -o ${workdir}/jobs/index_referece.out -l nodes=1:ppn=40,mem=180gb,walltime=1:00:00:00 -N index_referece ${workdir}/jobs/index_referece_Ave.sh

#####
# Map reads against Aeromonas salmonicida genome
#####

#jobs/map_to_referece.sh
module load tools samtools/1.11 bowtie2/2.4.2 pigz/2.3.4
#Aeromonas salmonicida
bowtie2 --threads 40 --reorder --rg-id ${sample} -x ${workdir}/reference/aeromonas_salmonicida -q ${workdir}/reads/${sample}.fastq -S ${workdir}/map/${sample}.Asa.sam
#Aeromonas veronii
bowtie2 --threads 40 --reorder --rg-id ${sample} -x ${workdir}/reference/aeromonas_veronii -q ${workdir}/reads/${sample}.fastq -S ${workdir}/map/${sample}.Ave.sam

sample=ES437
sample=ES470
sample=ES466
sample=ES166
sample=ES172
sample=ES174
sample=ES190
sample=ES244
sample=ES264
sample=ES354
sample=ES365
sample=ES441
sample=ES468
sample=ES472
sample=ES610
sample=ES612
sample=ES619
sample=ES621

qsub -V -A ku-cbd -W group_list=ku-cbd -v "workdir=${workdir},sample=${sample}"  -d `pwd` -e ${workdir}/jobs/map_to_referece.err -o ${workdir}/jobs/map_to_referece.out -l nodes=1:ppn=40,mem=180gb,walltime=1:00:00:00 -N map_to_referece_${sample} ${workdir}/jobs/map_to_referece.sh

#####
# Run iRep bPTR
#####

#jobs/bptr.sh
module load irep/1.10
#Aeromonas salmonicida
bPTR -f ${workdir}/reference/GCF_012931585.1_ASM1293158v1_genomic.fna -s ${workdir}/map/${sample}.Asa.sam -o ${workdir}/bptr/${sample}.Asa.tsv -plot ${workdir}/bptr/${sample}.Asa.pdf -m coverage
#Aeromonas veronii
bPTR -f ${workdir}/reference/GCF_008693705.1_ASM869370v1_genomic.fna -s ${workdir}/map/${sample}.Ave.sam -o ${workdir}/bptr/${sample}.Ave.tsv -plot ${workdir}/bptr/${sample}.Ave.pdf -m coverage

sample=ES437
sample=ES470
sample=ES466
sample=ES166
sample=ES172
sample=ES174
sample=ES190
sample=ES244
sample=ES264
sample=ES354
sample=ES365

sample=ES441
sample=ES468
sample=ES472
sample=ES610
sample=ES612
sample=ES619
sample=ES621
qsub -V -A ku-cbd -W group_list=ku-cbd -v "workdir=${workdir},sample=${sample}"  -d `pwd` -e ${workdir}/jobs/bptr.err -o ${workdir}/jobs/bptr.out -l nodes=1:ppn=40,mem=180gb,walltime=1:00:00:00 -N bptr_${sample} ${workdir}/jobs/bptr.sh


cat ${workdir}/bptr/${sample}.Ave.tsv
