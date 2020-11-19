### R

#Before running it, change files path and decide the level of taxonomic assigment:

library(dada2)
packageVersion("dada2") # 1.11.5 when this was put together

setwd("~/fishingbats/many_species/fishing")###Change it according to the study
setwd("~/fishingbats/many_species/others_capaccinii")
setwd("~/fishingbats/many_species/evie")

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
#ls *_1_filtered.fq.gz | sed "s/_1_filtered\.fq\.gz//g" > samples1
#samples1 <- scan("samples1", what="character")
#filtered_forward_reads <- paste0(samples1, "_1_filtered.fq.gz")
#filtered_reverse_reads <- paste0(samples1, "_2_filtered.fq.gz")


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
saveRDS(filtered_out, "filtered_out.RData")

#look at the filtered reads
#plotQualityProfile(filtered_forward_reads)
#plotQualityProfile(filtered_reverse_reads)
#plotQualityProfile(filtered_reverse_reads[17:20])
#plotQualityProfile(filtered_forward_reads[1:2])

#Generating an error model of our data by learning the specific error-signature of our dataset
err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE)
saveRDS(err_forward_reads, "err_forward_reads.RData")
saveRDS(err_reverse_reads, "err_reverse_reads.RData")

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
saveRDS(derep_forward, "derep_forward.RData")
saveRDS(derep_reverse, "derep_reverse.RData")

#Inferring ASVs
#It does this by incorporating the consensus quality profiles and abundances of each unique sequence, and then figuring out if each sequence more likely to be of biological origin or more likely to be spurious.
dada_forward <- dada(derep_forward, err=err_forward_reads, multithread=TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, multithread=TRUE)
saveRDS(dada_forward, "dada_forward.RData")
saveRDS(dada_reverse, "dada_reverse.RData")

#Merging forward and reverse reads
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                               derep_reverse, minOverlap=5, verbose=TRUE)
saveRDS(merged_amplicons, "merged_amplicons.RData")
#maxMismatch is by default 0 but you can make it say 1 or 2 if you’re finding that a lot of your forward and reverse reads are not merging
#minOverlap option has been change because thhey were sequenced at 250 bp
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(merged_amplicons)
saveRDS(seqtab, "seqtab.rds") 

#Generating a count table
table(nchar(getSequences(seqtab)))

#### MADS ####
#### NOW YOU SHOULD FIRST OPEN EACH COUNT TABLE AND MERGE THE RUNS ####
# Merge multiple runs (if necessary)
st1 <- readRDS("~/fishingbats/many_species/others_capaccinii/seqtab.rds")
st2 <- readRDS("~/fishingbats/many_species/fishing/seqtab.rds")
st3 <- readRDS("~/fishingbats/many_species/evie/seqtab.rds")
st.all <- mergeSequenceTables(st1, st2, st3)

########

#Chimera identification
seqtab.nochim <- removeBimeraDenovo(st.all, multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, "seqtab.nochim.RData")

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
taxa <- assignTaxonomy(seqtab.nochim, "~/fishingbats/many_species/silva_nr_v132_train_set.fa.gz", tryRC=T)
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
save.image("all_enviroment.RData")


