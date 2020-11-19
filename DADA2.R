### R

library(dada2)
packageVersion("dada2") # 1.11.5 when this was put together

## setting variables ##
samples <- scan("samples", what="character")
forward_reads <- paste0(samples, "_1_trimmed.fq.gz")
reverse_reads <- paste0(samples, "_2_trimmed.fq.gz")
filtered_forward_reads <- paste0(samples, "_1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples, "_2_filtered.fq.gz")

#Quality trimming/filtering
plotQualityProfile(reverse_reads[1:4])
plotQualityProfile(forward_reads[1:4])

#Filtering
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=c(223,223), compress=TRUE, multithread=TRUE)
saveRDS(filtered_out, "filtered_out.RData")
plotQualityProfile(filtered_reverse_reads[1:4])
plotQualityProfile(filtered_forward_reads[1:4])

#Generating an error model of the data by learning the specific error-signature of our dataset
err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE)
saveRDS(err_forward_reads, "err_forward_reads.RData")
saveRDS(err_reverse_reads, "err_reverse_reads.RData")

#Visualize how well the estimated error rates match up with the observed
plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)

#Dereplication
derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples 
saveRDS(derep_forward, "derep_forward.RData")
saveRDS(derep_reverse, "derep_reverse.RData")

#Inferring ASVs
dada_forward <- dada(derep_forward, err=err_forward_reads, multithread=TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, multithread=TRUE)
saveRDS(dada_forward, "dada_forward.RData")
saveRDS(dada_reverse, "dada_reverse.RData")

#Merging forward and reverse reads
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                               derep_reverse, minOverlap=5, verbose=TRUE)
saveRDS(merged_amplicons, "merged_amplicons.RData")

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(merged_amplicons)
saveRDS(seqtab, "seqtab.rds") 

#Generating a count table
table(nchar(getSequences(seqtab)))

# Merge multiple runs (if necessary)
st1 <- readRDS("~/fishingbats/many_species/others_capaccinii/seqtab.rds")
st2 <- readRDS("~/fishingbats/many_species/fishing/seqtab.rds")
st3 <- readRDS("~/fishingbats/many_species/evie/seqtab.rds")
st.all <- mergeSequenceTables(st1, st2, st3)

#Chimera identification
seqtab.nochim <- removeBimeraDenovo(st.all, multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, "seqtab.nochim.RData")

#Overview of counts throughout
# set a little function
getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))
write.table(summary_tab, "summary_reads_table.txt", sep="\t", quote=F)

#Taxonomy assigment
taxa <- assignTaxonomy(seqtab.nochim, "~/fishingbats/many_species/silva_nr_v132_train_set.fa.gz", tryRC=T)
saveRDS(taxa, "taxa.RData")

#Extracting the standard goods from DADA2
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)


