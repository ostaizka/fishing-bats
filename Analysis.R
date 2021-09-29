asv_count <- read.table("~/ASVs_counts.tsv",row.names=1)
asv_taxonomy <- read.table("~/ASVs_taxonomy.tsv",row.names=1)
metadata <- read.csv("~/samples.csv",sep=";",header=TRUE,row.names=1)
hierarchy <- read.csv("~/samples.csv",sep=";",header=TRUE)

library(phyloseq)
asv_taxonomy.matrix <- as.matrix(asv_taxonomy)

count_phy <- otu_table(asv_count, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(metadata)
TAX = tax_table(asv_taxonomy.matrix)
physeq = phyloseq(count_phy, TAX, sample_info_tab_phy)

############# DECONTAM ################
library(decontam)
physeq1 <- subset_samples(physeq, Run=="First")
physeq1 <- prune_taxa(taxa_sums(physeq1)>0, physeq1)
sample_data(physeq1)$is.neg <- sample_data(physeq1)$Sample_type == "Control" 
contamdf.prev1 <- isContaminant(physeq1, method="prevalence", neg="is.neg")
table(contamdf.prev1$contaminant)
head(which(contamdf.prev1$contaminant))
# getting IDs of those identified as likely contaminants
contam_asvs1 <- row.names(contamdf.prev1[contamdf.prev1$contaminant == TRUE, ])
# We can see this by looking at their taxonomic designations in our tax table
asv_taxonomy.matrix[row.names(asv_taxonomy.matrix) %in% contam_asvs1, ]
# making new phyloseq object without ASV contaminats
goodTaxa1 <- setdiff(taxa_names(physeq1), contam_asvs1)
physeq_no_contam1 <- prune_taxa(goodTaxa1, physeq1)
# making new phyloseq object without blanks/controls samples 
physeq_no_blanks1 <- subset_samples(physeq_no_contam1, Sample_type %in% c("Sample"))
physeq_clean1 <- prune_taxa(taxa_sums(physeq_no_blanks1)>0, physeq_no_blanks1)
head(physeq_clean1@otu_table)
# making new fasta file
contam_indices1 <- which(asv_fasta %in% paste0(">", contam_asvs1))
dont_want1 <- sort(c(contam_indices1, contam_indices1 + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want1]


#Second run
physeq2 <- subset_samples(physeq, Run=="Second")
physeq2 <- prune_taxa(taxa_sums(physeq2)>0, physeq2)
sample_data(physeq2)$is.neg <- sample_data(physeq2)$Sample_type == "Control" 
contamdf.prev2 <- isContaminant(physeq2, method="prevalence", neg="is.neg")
table(contamdf.prev2$contaminant)
head(which(contamdf.prev2$contaminant))
# getting IDs of those identified as likely contaminants
contam_asvs2 <- row.names(contamdf.prev2[contamdf.prev2$contaminant == TRUE, ])
# We can see this by looking at their taxonomic designations in our tax table: taxmat
asv_taxonomy.matrix[row.names(asv_taxonomy.matrix) %in% contam_asvs2, ]
# making new phyloseq object without ASV contaminats
goodTaxa2 <- setdiff(taxa_names(physeq2), contam_asvs2)
physeq_no_contam2 <- prune_taxa(goodTaxa2, physeq2)
# making new phyloseq object without blanks/controls samples
physeq_no_blanks2 <- subset_samples(physeq_no_contam2, Sample_type %in% c("Sample"))
physeq_clean2 <- prune_taxa(taxa_sums(physeq_no_blanks2)>0, physeq_no_blanks2)
head(physeq_clean2@otu_table)
# new fasta file
contam_indices2 <- which(asv_fasta_no_contam %in% paste0(">", contam_asvs2))
dont_want2 <- sort(c(contam_indices2, contam_indices2 + 1))
asv_fasta_clean <- asv_fasta_no_contam[- dont_want2]

#Merge phyloseq
physeqall <- merge_phyloseq(physeq_clean1, physeq_clean2)

#extract the asv table from the phyloseq
asv.table <- data.frame(physeqall@otu_table)

########### DIVERSITY ANALYSIS (hilldiv) #########
##Remove samples with less than 10000 reads for the count table
library(hilldiv)
counts_filtdepth <- depth_filt(asv.table,10000)
col_asv.table <- colnames(asv.table)
col_counts_filtdepth <- colnames(counts_filtdepth)
#colnames of samples remove for filtering, check for equality
samples_out <- setdiff(col_asv.table, col_counts_filtdepth)

##Remove all OTUs with less copies than 0.01% of the total number of reads of each sample
counts_filtcopy <- copy_filt(counts_filtdepth, 0.0001)
row.filtered <- row.names(counts_filtdepth)
row.filtcopy <- row.names(counts_filtcopy)
#asvs remove after filtering
asv_out <- setdiff(row.filtered, row.filtcopy)

#new phyloseq
count_phy <- otu_table(counts_filtcopy, taxa_are_rows=T)
physeq = phyloseq(count_phy, TAX, sample_info_tab_phy)

###TAXONOMY filtering
##Remove ASVs identified as Eukaryota or Archaea, and Bacteria lacking Phylum info in the taxonomy file
physeqall_bacteria <- subset_taxa(physeq, Kingdom =="Bacteria")
physeqall_phylum <- subset_taxa(physeqall_bacteria, Phylum != "NA")
physeqall_class <- subset_taxa(physeqall_phylum, Class != "NA")



########### Bar plot: gut microbiota diversity of all bats #########
pseq.rel <- microbiome::transform(physeqall_class, "compositional")
physeqall.comp.rare <- aggregate_rare(pseq.rel, level = "Phylum", detection = 0.1/100, prevalence = 20/100, include.lowest = TRUE)

top_phylum <- plot_composition(physeqall.comp.rare, sample.sort = NULL, 
                               otu.sort = NULL,
                               x.label = "Species", 
                               plot.type = "barplot", 
                               verbose = FALSE) + 
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Top phyla") + 
  theme_bw() + scale_fill_brewer("Phylum", palette = "Paired")
print(top_phylum + theme(axis.text.x = element_text(angle = 90,size=5)))


########### DIVERSITY ANALYSIS (hilldiv) #########
#get the tables
asv.table.clean <- data.frame(physeqall_class@otu_table)
metadata.clean <- data.frame(physeqall_class@sam_data)
taxonomy.clean <- data.frame(physeqall_class@tax_table)
library(tibble)
library(dplyr)
hierarchy <- rownames_to_column(metadata.clean, var = "Sample")
identical(sort(colnames(asv.table.clean)),sort(as.character(hierarchy[,1])))
#Filter hierarchy table
hierarchy <- hierarchy[which(hierarchy[,1] %in% colnames(asv.table.clean)),]
#Tree
tree <- read.tree("ASVs.align.raxml.bestTree")
tree <- force.ultrametric(tree, method = "extend")
#Show basic statistics
tree
str(tree)
#hilldiv incorporates the function match_data(), which enables checking in a straightforward way whether the OTU/ASV names in the count table and the tree match, and subsetting count tables or trees accordingly.
tree.clean <- match_data(asv.table.clean,tree,output="tree")
match_data(asv.table.clean,tree.clean)

library(hilldiv)
#q=0
divq0 <- div_test(asv.table.clean,qvalue=0,hierarchy=hierarchy[,c(1,2)],posthoc=TRUE)
#q=1
divq1 <- div_test(asv.table.clean,qvalue=1,hierarchy=hierarchy[,c(1,2)],posthoc=TRUE)
#Phylogenies
divq0phylo <- div_test(asv.table.clean,qvalue=0,hierarchy=hierarchy[,c(1,2)],tree=tree.clean,posthoc=TRUE)
divq1phylo <- div_test(asv.table.clean,qvalue=1,hierarchy=hierarchy[,c(1,2)],tree=tree.clean,posthoc=TRUE)

################## Dissimilarity ##################
pairdis.q0 <- pair_dis(asv.table.clean,qvalue=0,hierarchy=hierarchy[,c(1,2)])
dis_nmds(pairdis.q0$L1_UqN,hierarchy=hierarchy[,c(1,2)], centroids = TRUE)
pairdis.q1 <- pair_dis(asv.table.clean,qvalue=1,hierarchy=hierarchy[,c(1,2)])
dis_nmds(pairdis.q1$L1_UqN,hierarchy=hierarchy[,c(1,2)], centroids = TRUE)
pairdis.q0.Phy <- pair_dis(asv.table.clean,qvalue=0,hierarchy=hierarchy[,c(1,2)],tree=tree.clean)
dis_nmds(pairdis.q0.Phy$L1_UqN,hierarchy=hierarchy[,c(1,2)], centroids = TRUE)
pairdis.q1.Phy <- pair_dis(asv.table.clean,qvalue=1,hierarchy=hierarchy[,c(1,2)],tree=tree.clean)
dis_nmds(pairdis.q1.Phy$L1_UqN,hierarchy=hierarchy[,c(1,2)], centroids = TRUE)

################## Permanova ##################
library(vegan)
# extract the values from pairdis
u0n <- pairdis.q0$L1_UqN
# coerce a object to a distance matrix
u0n.dist <- as.dist(u0n)
# level of homogeneity of dispersion
ps.disper.u0n.diet <- betadisper(u0n.dist, hierarchy$diet)
ps.disper.u0n.habitat <- betadisper(u0n.dist, hierarchy$habitat)
ps.disper.u0n.Species <- betadisper(u0n.dist, hierarchy$Species)
ps.disper.u0n.Family <- betadisper(u0n.dist, hierarchy$Family)
permutest(ps.disper.u0n.diet, pairwise = TRUE)
permutest(ps.disper.u0n.habitat, pairwise = TRUE)
permutest(ps.disper.u0n.Species, pairwise = TRUE)
permutest(ps.disper.u0n.Family, pairwise = TRUE)
# adonis
adonis(u0n.dist ~ diet + Family + habitat + Species, data =metadata.clean, permutations = 999)
library(mctoolsr)
calc_pairwise_permanovas(u0n.dist, metadata.clean, "habitat")
calc_pairwise_permanovas(u0n.dist, metadata.clean, "Species")
calc_pairwise_permanovas(u0n.dist, metadata.clean, "Family")
#q1
u1n <- pairdis.q1$L1_UqN
u1n.dist <- as.dist(u1n)
ps.disper.u1n.diet <- betadisper(u1n.dist, hierarchy$diet)
ps.disper.u1n.habitat <- betadisper(u1n.dist, hierarchy$habitat)
ps.disper.u1n.Species <- betadisper(u1n.dist, hierarchy$Species)
ps.disper.u1n.Family <- betadisper(u1n.dist, hierarchy$Family)
permutest(ps.disper.u1n.diet, pairwise = TRUE)
permutest(ps.disper.u1n.habitat, pairwise = TRUE)
permutest(ps.disper.u1n.Species, pairwise = TRUE)
permutest(ps.disper.u1n.Family, pairwise = TRUE)
# adonis
adonis(u1n.dist ~ diet + Family + habitat + Species, data =metadata.clean, permutations = 999)
library(mctoolsr)
calc_pairwise_permanovas(u1n.dist, metadata.clean, "habitat")
calc_pairwise_permanovas(u1n.dist, metadata.clean, "Species")
calc_pairwise_permanovas(u1n.dist, metadata.clean, "Family")

################## Bubble plot ##################

Mric <- subset_samples(pseq.rel, Group == "Myotis_ricketti")
Mric <- prune_taxa(taxa_sums(Mric)>0, Mric)#remove ASV==0
Mric.abund <- subset_taxa(Mric, Genus %in% c("Aeromonas", "Cetobacterium", "Vibrio", "Plesiomonas", "Paraclostridium", "Photobacterium", "Candidatus_Rhabdochlamydia", "Mydepth_filtcoplasma", "Romboutsia", "Chlamydia", "Orbus"))
Mric.genus <- aggregate_taxa(Mric.abund, 'Genus')
Mric.genus.prop <- rowMeans(Mric.genus@otu_table)*100
Mric.genus.prop.table <- as.data.frame(Mric.genus.prop)

Mcap <- subset_samples(pseq.rel, Group == "Myotis_capaccinii")
Mcap.abund <- subset_taxa(Mcap, Genus %in% c("Aeromonas", "Cetobacterium", "Vibrio", "Plesiomonas", "Paraclostridium", "Photobacterium", "Candidatus_Rhabdochlamydia", "Mycoplasma", "Romboutsia", "Chlamydia", "Orbus"))
Mcap.genus <- aggregate_taxa(Mcap.abund, 'Genus')
Mcap.genus.prop <- rowMeans(Mcap.genus@otu_table)*100
Mcap.genus.prop.table <- as.data.frame(Mcap.genus.prop)

#Merge both tables
Mric.genus.prop.table$Bacteria1 <- rownames(Mric.genus.prop.table)
Mcap.genus.prop.table$Bacteria2 <- rownames(Mcap.genus.prop.table)
bats_merge <- merge(Mric.genus.prop.table, Mcap.genus.prop.table[,c("Bacteria2","Mcap.genus.prop")], by.x="Bacteria1", by.y="Bacteria2", all.x=TRUE, all.y=TRUE) 

FMcap <- subset_samples(pseq.rel, Group == "Fishing_capaccinii")
FMcap <- prune_taxa(taxa_sums(FMcap)>0, FMcap)#remove ASV==0
FMcap.abund <- subset_taxa(FMcap, Genus %in% c("Aeromonas", "Cetobacterium", "Vibrio", "Plesiomonas", "Paraclostridium", "Photobacterium", "Candidatus_Rhabdochlamydia", "Mycoplasma", "Romboutsia", "Chlamydia", "Orbus"))
FMcap.genus <- aggregate_taxa(FMcap.abund, 'Genus')
FMcap.genus.prop <- rowMeans(FMcap.genus@otu_table)*100
FMcap.genus.prop.table <- as.data.frame(FMcap.genus.prop)
#Merge tables
FMcap.genus.prop.table$Bacteria2 <- rownames(FMcap.genus.prop.table)
bats_merge <- merge(bats_merge, FMcap.genus.prop.table[,c("Bacteria2","FMcap.genus.prop")], by.x="Bacteria1", by.y="Bacteria2", all.x=TRUE, all.y=TRUE) 

Mvi <- subset_samples(pseq.rel, Group == "Myotis_vivesi")
Mvi <- prune_taxa(taxa_sums(Mvi)>0, Mvi)#remove ASV==0
Mvi.abund <- subset_taxa(Mvi, Genus %in% c("Aeromonas", "Cetobacterium", "Vibrio", "Plesiomonas", "Paraclostridium", "Photobacterium", "Candidatus_Rhabdochlamydia", "Mycoplasma", "Romboutsia", "Chlamydia", "Orbus"))
Mvi.genus <- aggregate_taxa(Mvi.abund, 'Genus')
Mvi.genus.prop <- rowMeans(Mvi.genus@otu_table)*100
Mvi.genus.prop.table <- as.data.frame(Mvi.genus.prop)
#Merge tables
Mvi.genus.prop.table$Bacteria2 <- rownames(Mvi.genus.prop.table)
bats_merge <- merge(bats_merge, Mvi.genus.prop.table[,c("Bacteria2","Mvi.genus.prop")], by.x="Bacteria1", by.y="Bacteria2", all.x=TRUE, all.y=TRUE) 

Nlep <- subset_samples(pseq.rel, Group == "Noctilio_leporinus")
Nlep <- prune_taxa(taxa_sums(Nlep)>0, Nlep)#remove ASV==0
Nlep.abund <- subset_taxa(Nlep, Genus %in% c("Aeromonas", "Cetobacterium", "Vibrio", "Plesiomonas", "Paraclostridium", "Photobacterium", "Candidatus_Rhabdochlamydia", "Mycoplasma", "Romboutsia", "Chlamydia", "Orbus"))
Nlep.genus <- aggregate_taxa(Nlep.abund, 'Genus')
Nlep.genus.prop <- rowMeans(Nlep.genus@otu_table)*100
Nlep.genus.prop.table <- as.data.frame(Nlep.genus.prop)
#Merge tables
Nlep.genus.prop.table$Bacteria2 <- rownames(Nlep.genus.prop.table)
bats_merge <- merge(bats_merge, Nlep.genus.prop.table[,c("Bacteria2","Nlep.genus.prop")], by.x="Bacteria1", by.y="Bacteria2", all.x=TRUE, all.y=TRUE) 

Insect <- subset_samples(pseq.rel, Group == "Insectivorous")
Insect <- prune_taxa(taxa_sums(Insect)>0, Insect)#remove ASV==0
Insect.abund <- subset_taxa(Insect, Genus %in% c("Aeromonas", "Cetobacterium", "Vibrio", "Plesiomonas", "Paraclostridium", "Photobacterium", "Candidatus_Rhabdochlamydia", "Mycoplasma", "Romboutsia", "Chlamydia", "Orbus"))
Insect.genus <- aggregate_taxa(Insect.abund, 'Genus')
Insect.genus.prop <- rowMeans(Insect.genus@otu_table)*100
Insect.genus.prop.table <- as.data.frame(Insect.genus.prop)
#Merge tables
Insect.genus.prop.table$Bacteria2 <- rownames(Insect.genus.prop.table)
bats_merge <- merge(bats_merge, Insect.genus.prop.table[,c("Bacteria2","Insect.genus.prop")], by.x="Bacteria1", by.y="Bacteria2", all.x=TRUE, all.y=TRUE) 

#fish samples
sample <- read.csv("~/Sampleinfo_fish.csv",sep=";", header=TRUE, row.names = 1)
asv <- read.csv("~/ASVtable_fish.csv", header=TRUE, row.names = 1)
taxonomy <- read.csv("~/taxonomy_fish.csv", header=TRUE, row.names = 1)
count_phy <- otu_table(asv, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample)
taxonomy.matrix <- as.matrix(taxonomy)

TAX = tax_table(taxonomy.matrix)
physeq.fish = phyloseq(count_phy, TAX, sample_info_tab_phy)
pseq.rel.fish <- microbiome::transform(physeq.fish, "compositional")

gambusia <- subset_samples(pseq.rel.fish, Group == "Gamb_content")
gambusia.abund <- subset_taxa(gambusia, Genus %in% c("Aeromonas", "Cetobacterium", "Vibrio", "Plesiomonas", "Paraclostridium", "Photobacterium", "Candidatus_Rhabdochlamydia", "Mycoplasma", "Romboutsia", "Chlamydia", "Orbus"))
gambusia.genus <- aggregate_taxa(gambusia.abund, 'Genus')
gambusia.genus.prop <- rowMeans(gambusia.genus@otu_table)*100
gambusia.genus.prop.table <- as.data.frame(gambusia.genus.prop)
#Merge tables
gambusia.genus.prop.table$Bacteria2 <- rownames(gambusia.genus.prop.table)
bats_merge <- merge(bats_merge, gambusia.genus.prop.table[,c("Bacteria2","gambusia.genus.prop")], by.x="Bacteria1", by.y="Bacteria2", all.x=TRUE, all.y=TRUE) 

#Plot it
x <- c("Aeromonas", "Cetobacterium", "Plesiomonas", "Paraclostridium", "Photobacterium", "Romboutsia", "Chlamydia", "Vibrio", "Candidatus_Rhabdochlamydia", "Mycoplasma", "Orbus")
bats_merge1 <- bats_merge %>%
  slice(match(x, Bacteria1))

library(reshape)
pcm = melt(bats_merge1, id = c("Bacteria1"))
pcm %>%
  slice(match(x, Bacteria1))
colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
             "#F6AE2D","#86BBD8")
pcm$Bacteria1 <- factor(pcm$Bacteria1,levels=unique(pcm$Bacteria1))

pdf("~/bubble_plot.pdf",width=10, height=6)
ggplot(pcm, aes(x = Bacteria1, y = variable)) + 
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 60), range = c(1,16), breaks = c(0.1,10,30,60)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 10, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", size = 10), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") +  
  scale_fill_manual(values = colours, guide = FALSE) + 
  scale_y_discrete(limits = rev(levels(pcm$variable))) 
dev.off() 

################## Enrichment ##################
library(DESeq2)
###Arthropodivorous vs. piscivorous
diagdds = phyloseq_to_deseq2(physeqall_class, ~ diet)
diagdds <- estimateSizeFactors(diagdds, type="poscounts",locfunc=genefilter::shorth)#count table is sparse with many zeroes. In that case you'd need to use a specific function first that is equipped to deal with that.
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
diagdds.ins.fis <- results(diagdds, alpha=0.01, contrast=c("fish_gm", "Insectivorous", "Fishing"))
sigtab_diagdds.ins.fis <- diagdds.ins.fis[which(diagdds.ins.fis$padj < 0.01), ]
sigtab_diagdds.ins.fis_with_tax <- cbind(as(sigtab_diagdds.ins.fis, "data.frame"), as(tax_table(TAX)[row.names(sigtab_diagdds.ins.fis), ], "matrix"))
sigtab_diagdds.ins.fis_with_tax[order(sigtab_diagdds.ins.fis_with_tax$baseMean, decreasing=T), ]
deseq2_ins_fis <- as.data.frame(sigtab_diagdds.ins.fis_with_tax)

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

x = tapply(sigtab_diagdds.ins.fis_with_tax$log2FoldChange, sigtab_diagdds.ins.fis_with_tax$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_diagdds.ins.fis_with_tax$Phylum = factor(as.character(sigtab_diagdds.ins.fis_with_tax$Phylum), levels=names(x))
x = tapply(sigtab_diagdds.ins.fis_with_tax$log2FoldChange, sigtab_diagdds.ins.fis_with_tax$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_diagdds.ins.fis_with_tax$Genus = factor(as.character(sigtab_diagdds.ins.fis_with_tax$Genus), levels=names(x))
ggplot(sigtab_diagdds.ins.fis_with_tax, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + geom_hline(yintercept=0) + coord_flip()
