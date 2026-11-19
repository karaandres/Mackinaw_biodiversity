### This code creates phyloseq objects for the 12S, COI, & 18S sequences for the full study
### Includes data filtering, rarefaction, and basic data exploration
### Inputs: ASV table, taxonomy, and phylogeny from qiime2
### Outputs: final phyloseq objects (rarefied datasets), figures
### Last updated 8 May 2026 by Kara Andres (kjandr2@ilstu.edu)

### Clear work environment and load packages
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
library(ggplot2)
library(phyloseq) ## NOTE: must be installed BiocManager::install("phyloseq"); 'igraph' must be installed first
library(qiime2R) ## devtools::install_github("jbisanz/qiime2R")
library(reshape2)
library(decontam) ## BiocManager::install("decontam")
library(Biostrings)
library(stringr)
library(phylosmith) ## devtools::install_github("schuyler-smith/phylosmith")
library(RColorBrewer)
library(MiscMetabar)
library(colorspace)
library(vegan)
library(forcats)
library(microbiome) ## BiocManager::install("microbiome")
library(dplyr)
library(tidyverse)

##############################################################################
################ Starting with 12S MiFish ####################################
##############################################################################

### Produce phyloseq object
ASVs <- read_qza("/Users/kbja10/Library/CloudStorage/OneDrive-ILStateUniversity/ISU/Research/Mackinaw_biodiversity/data/12S/output_29Apr2026/table.qza")
taxonomy <- read_qza("/Users/kbja10/Library/CloudStorage/OneDrive-ILStateUniversity/ISU/Research/Mackinaw_biodiversity/data/12S/output_29Apr2026/taxonomy.qza")
R_tree <- read_qza("/Users/kbja10/Library/CloudStorage/OneDrive-ILStateUniversity/ISU/Research/Mackinaw_biodiversity/data/12S/output_29Apr2026/rooted-tree.qza")
seqs <- read_qza("/Users/kbja10/Library/CloudStorage/OneDrive-ILStateUniversity/ISU/Research/Mackinaw_biodiversity/data/12S/output_29Apr2026/rep-seqs-dada2.qza")
metadata <- read_q2metadata("/Users/kbja10/Library/CloudStorage/OneDrive-ILStateUniversity/ISU/Research/Mackinaw_biodiversity/data/Mackinaw_eDNA_metadata_R.txt")
rownames(metadata) <- metadata$Sample_name_12S
names(metadata)[names(metadata) == "SampleID"] <- "Site_name"
metadata$Sample_name <- as.factor(metadata$Sample_name)
metadata$Sample_name_12S <- as.factor(metadata$Sample_name_12S)
metadata$Site_name <- as.factor(metadata$Site_name)
metadata$Site_abbreviation <- as.factor(metadata$Site_abbreviation)
metadata$Sample_type <- as.factor(metadata$Sample_type)
phyloseq_dat <- phyloseq(otu_table(ASVs$data,taxa_are_rows = TRUE), 
                         tax_table(as.matrix(parse_taxonomy(taxonomy$data))), 
                         phy_tree(R_tree$data),
                         refseq(DNAStringSet(seqs$data)),
                         sample_data(metadata))
phyloseq_dat # 444 taxa and 70 samples
sum(sample_sums(phyloseq_dat)) # 4,840,823 reads

### Remove ASVs that are not of target length
seq_lengths <- width(refseq(phyloseq_dat))
hist(seq_lengths) # Remove above 220
phyloseq_dat <- prune_taxa(seq_lengths <= 220, phyloseq_dat)
phyloseq_dat # 301 taxa and 70 samples

### Remove sequences not in target classes (Actinopteri & Hyperoartia)
unique(tax_table(phyloseq_dat)[, "Class"])
phyloseq_dat_offtarget <- subset_taxa(phyloseq_dat,
                                      (Class %in% c("Amphibia", "Sarcopterygii", "Insecta", "Aves", "Mammalia")))
View(phyloseq_dat_offtarget@tax_table)
# Off-target seqs include bats, rodents, minks, frogs, turtles, birds, deer, otter
# Subset to just target taxa
phyloseq_dat <- subset_taxa(phyloseq_dat,
                            Class %in% c("Actinopteri", "Hyperoartia"))
phyloseq_dat # 228 taxa and 70 samples

### Inspect negative controls 
phyloseq_blank <- subset_samples(phyloseq_dat, Sample_type=="blank")
phyloseq_blank <- prune_taxa(taxa_sums(phyloseq_blank) > 0, phyloseq_blank)
dim(phyloseq_blank@otu_table) # 58 ASVs present in 16 negative controls
sum(rowSums(unique(phyloseq_blank@otu_table))) # 33,716 reads total 
phyloseq_blank_table <- phyloseq_blank %>% psmelt()
p <- ggplot(phyloseq_blank_table, aes(x=Sample_name, y=Abundance, fill=Genus)) + 
  geom_col(width=0.9) +
  ylab("Read Count") + xlab("Sample name") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text=element_text(size=8), axis.title=element_text(size=14),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(size=8))
p
# Almost all reads are in one particularly bad extraction blank 

### Identify and remove contaminants by decontam
sample_data(phyloseq_dat)$is.neg <- sample_data(phyloseq_dat)$Sample_type == "blank"
contam <- isContaminant(phyloseq_dat, method="prevalence", neg="is.neg", threshold=0.1)
contaminants <- rownames(contam)[which(contam$contaminant)] 
contaminants # 0 contaminant sequences identified
phyloseq_decontam <- prune_taxa(!rownames(phyloseq_dat@otu_table) %in% contaminants, phyloseq_dat)

### Remove negative controls
phyloseq_clean <- subset_samples(phyloseq_decontam, Sample_type!="blank")
phyloseq_clean # 228 taxa and 54 samples

### Inspect sequencing depth per sample
hist(sample_sums(phyloseq_clean@otu_table), xlab="Number of reads", main="Histogram of number of reads per sample", breaks=50)
read_depth <- data.frame(Sample_names=sample_data(phyloseq_clean)$Sample_name,
                         Site_abbreviation=sample_data(phyloseq_clean)$Site_abbreviation,
                         read_depth=sample_sums(phyloseq_clean))
p <- ggplot(read_depth, aes(x=Site_abbreviation, y=read_depth, fill=Site_abbreviation)) + 
  geom_boxplot(outlier.shape=NA, outlier.colour=NA) + 
  geom_point(shape=16, size=2, color="gray") + 
  labs(x="Sample name", y="Read Depth") + 
  theme_bw() +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        legend.position = "none")
p

### Determine rarefaction depth with ASV accumulation curves by Site
accum_curves <- accu_plot(phyloseq_clean, "Site_abbreviation", add_nb_seq=TRUE, step=1000) + 
  theme_bw()
accum_curves

### Most samples plateau by 1,000 reads -- target rarefaction depth ###

### Identify samples with reads < 1,000
sample_sums(prune_samples(sample_sums(phyloseq_clean)<1000, phyloseq_clean))
length(sort(sample_sums(prune_samples(sample_sums(phyloseq_clean)<1000, phyloseq_clean)))) # 4 samples below 1,000 reads
sample_filter <- colnames(phyloseq_clean@otu_table)[sample_sums(phyloseq_clean@otu_table)<1000]
phyloseq_1000 <- subset_samples(phyloseq_clean,!sample_names(phyloseq_clean) %in% sample_filter)
phyloseq_1000 # 228 taxa and 50 samples

### Produce final phyloseq after rarefaction
set.seed(1000)
phyloseq_rarefied <- rarefy_even_depth(phyloseq_1000, sample.size=min(sample_sums(phyloseq_1000)))
phyloseq_rarefied <- prune_taxa(taxa_sums(phyloseq_rarefied)>0, phyloseq_rarefied)
phyloseq_rarefied # 173 taxa and 50 samples
sum(sample_sums(phyloseq_rarefied)) # 52,250 reads
phyloseq_rarefied_12S <- phyloseq_rarefied
saveRDS(phyloseq_rarefied_12S,file = "phyloseq_rarefied_12S.rds")

### Sample size per site 
sample_size <- as.data.frame(with(sample_data(phyloseq_clean), table(Site_abbreviation)))
clean_n <- as.data.frame(with(sample_data(phyloseq_rarefied), table(Site_abbreviation)))
merge(sample_size, clean_n, by="Site_abbreviation", all=TRUE)

### Pie chart of Family: read abundance across all samples
ps_family <- tax_glom(phyloseq_rarefied_12S, taxrank = "Family")
my.cols <- c(brewer.pal(8,"Dark2"), "#C46254","bisque","darkorange", rainbow_hcl(12), rainbow(13), "azure3")
data <- ps_family %>%
  psmelt() %>%
  as_tibble() %>%
  group_by(Family) %>%
  summarise(Reads = sum(Abundance), .groups = "drop") %>%
  mutate(RelAb = Reads / sum(Reads)) %>%
  mutate(Family = ifelse(RelAb < 0.01, "Other", Family)) %>% # Combine families < 1% into "Other"
  group_by(Family) %>%
  summarise(RelAb = sum(RelAb), .groups = "drop") %>%
  arrange(desc(RelAb)) %>%
  mutate(Family = factor(Family, levels = Family))
pie <- ggplot(data, aes(x = "", y = RelAb, fill = Family)) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Family", values = my.cols) +
  coord_polar("y") +
  theme_void() +
  theme(legend.title = element_text(size = 18),legend.text = element_text(size = 16),
        legend.box.spacing = unit(0, "cm"))
pie

### Closer look at taxonomic proportions (Family) per site
familyabundance <- phyloseq_rarefied_12S %>%
  tax_glom(taxrank = "Family") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>%
  as_tibble()
familyabundance_avg <- familyabundance %>%
  group_by(Family, Site_name, Land_use_index) %>%
  summarize(avg_abundance = mean(Abundance), .groups = "drop") %>%
  mutate(Family = ifelse(avg_abundance < 0.01, "Other", Family)) %>%
  group_by(Family, Site_name, Land_use_index) %>%
  summarize(avg_abundance = sum(avg_abundance), .groups = "drop")
p <- ggplot(familyabundance_avg) +
  geom_col(mapping = aes(x=fct_reorder(Site_name, Land_use_index), y=avg_abundance, fill=fct_reorder(Family, avg_abundance, .desc=TRUE)), position="fill", show.legend=TRUE) +
  ylab("Proportion of 12S sequences") +
  scale_fill_manual(name="Family", values = my.cols) +
  xlab(NULL) + theme_bw() +
  theme(axis.text.y.left = element_text(size=15),
        axis.text.x = element_text(size=10, angle=45, hjust=1),
        axis.title.y = element_text(size=15))
p

### Diversity stats per Site
sample_data(phyloseq_rarefied_12S)$Site_name <- fct_reorder(sample_data(phyloseq_rarefied_12S)$Site_name,
                                                        sample_data(phyloseq_rarefied_12S)$Land_use_index)
plot_richness(phyloseq_rarefied_12S, x="Site_name", measures=c("Observed","Shannon"), color="Site_name") +
  geom_boxplot(alpha=0.6) + 
  scale_color_manual(values=my.cols) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=10))

### PCOA of Site
ord.PCoA.bray <- ordinate(phyloseq_rarefied_12S, method="PCoA", distance="bray")
plot_ordination(phyloseq_rarefied_12S, ord.PCoA.bray, type="samples", color="Site_name") +
  geom_point(size=3, stroke=1) +
  scale_color_manual(name="Site_name", values=my.cols) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=14))
plot_ordination(phyloseq_rarefied_12S, ord.PCoA.bray, type = "samples", color = "Land_use_index") +
  geom_point(size = 3, stroke = 1) +
  scale_color_gradientn(colours = c("darkgreen", "gold", "red"), name = "Land Use Index") +
  theme_bw() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20),
        legend.title = element_text(size = 20), legend.text = element_text(size = 14))

##############################################################################
################ Next with 18S ##############################################
##############################################################################

### Produce phyloseq object
ASVs <- read_qza("/Users/kbja10/Library/CloudStorage/OneDrive-ILStateUniversity/ISU/Research/Mackinaw_biodiversity/data/18S/output_29Apr2026/table.qza")
taxonomy <- read_qza("/Users/kbja10/Library/CloudStorage/OneDrive-ILStateUniversity/ISU/Research/Mackinaw_biodiversity/data/18S/output_29Apr2026/taxonomy.qza")
R_tree <- read_qza("/Users/kbja10/Library/CloudStorage/OneDrive-ILStateUniversity/ISU/Research/Mackinaw_biodiversity/data/18S/output_29Apr2026/rooted-tree.qza")
seqs <- read_qza("/Users/kbja10/Library/CloudStorage/OneDrive-ILStateUniversity/ISU/Research/Mackinaw_biodiversity/data/18S/output_29Apr2026/rep-seqs-dada2.qza")
metadata <- read_q2metadata("/Users/kbja10/Library/CloudStorage/OneDrive-ILStateUniversity/ISU/Research/Mackinaw_biodiversity/data/Mackinaw_eDNA_metadata_R.txt")
rownames(metadata) <- metadata$Sample_name_18S
names(metadata)[names(metadata) == "SampleID"] <- "Site_name"
metadata$Sample_name <- as.factor(metadata$Sample_name)
metadata$Sample_name_18S <- as.factor(metadata$Sample_name_18S)
metadata$Site_name <- as.factor(metadata$Site_name)
metadata$Site_abbreviation <- as.factor(metadata$Site_abbreviation)
metadata$Sample_type <- as.factor(metadata$Sample_type)
phyloseq_dat <- phyloseq(otu_table(ASVs$data,taxa_are_rows = TRUE), 
                         tax_table(as.matrix(parse_taxonomy(taxonomy$data))), 
                         phy_tree(R_tree$data),
                         refseq(DNAStringSet(seqs$data)),
                         sample_data(metadata))
phyloseq_dat # 5085 taxa and 70 samples
sum(sample_sums(phyloseq_dat)) # 4,763,084 reads

### Remove ASVs that are not of target length
seq_lengths <- width(refseq(phyloseq_dat))
hist(seq_lengths) # Remove above 180
phyloseq_dat <- prune_taxa(seq_lengths <= 180, phyloseq_dat)
phyloseq_dat # 5078 taxa and 70 samples

### Remove sequences not in target classes (Actinopteri & Hyperoartia)
unique(tax_table(phyloseq_dat)[, "Kingdom"])
phyloseq_dat <- subset_taxa(phyloseq_dat,
                            !(Kingdom %in% "Unassigned"))
phyloseq_dat # 4172 taxa and 70 samples

### Inspect negative controls 
phyloseq_blank <- subset_samples(phyloseq_dat, Sample_type=="blank")
phyloseq_blank <- prune_taxa(taxa_sums(phyloseq_blank) > 0, phyloseq_blank)
dim(phyloseq_blank@otu_table) # 715 ASVs present in 16 negative controls
sum(rowSums(unique(phyloseq_blank@otu_table))) # 1,287,370 reads total 
phyloseq_blank_table <- phyloseq_blank %>% psmelt()
p <- ggplot(phyloseq_blank_table, aes(x=Sample_name, y=Abundance, fill=Kingdom)) + 
  geom_col(width=0.9) +
  ylab("Read Count") + xlab("Sample name") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text=element_text(size=8), axis.title=element_text(size=14),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(size=8))
p
# Very high numbers of reads across all types of negative controls 

### Identify and remove contaminants by decontam
sample_data(phyloseq_dat)$is.neg <- sample_data(phyloseq_dat)$Sample_type == "blank"
contam <- isContaminant(phyloseq_dat, method="prevalence", neg="is.neg", threshold=0.1)
contaminants <- rownames(contam)[which(contam$contaminant)] 
length(contaminants) # 82 contaminant sequences identified
phyloseq_decontam <- prune_taxa(!rownames(phyloseq_dat@otu_table) %in% contaminants, phyloseq_dat)

### Remove negative controls
phyloseq_clean <- subset_samples(phyloseq_decontam, Sample_type!="blank")
phyloseq_clean # 4090 taxa and 54 samples

### Inspect sequencing depth per sample
hist(sample_sums(phyloseq_clean@otu_table), xlab="Number of reads", main="Histogram of number of reads per sample", breaks=50)
read_depth <- data.frame(Sample_names=sample_data(phyloseq_clean)$Sample_name,
                         Site_abbreviation=sample_data(phyloseq_clean)$Site_abbreviation,
                         read_depth=sample_sums(phyloseq_clean))
p <- ggplot(read_depth, aes(x=Site_abbreviation, y=read_depth, fill=Site_abbreviation)) + 
  geom_boxplot(outlier.shape=NA, outlier.colour=NA) + 
  geom_point(shape=16, size=2, color="gray") + 
  labs(x="Sample name", y="Read Depth") + 
  theme_bw() +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        legend.position = "none")
p

### Determine rarefaction depth with ASV accumulation curves by Site
accum_curves <- accu_plot(phyloseq_clean, "Site_abbreviation", add_nb_seq=TRUE, step=1000) + 
  theme_bw()
accum_curves

### Most samples plateau by 1,000 reads -- target rarefaction depth ###

### Identify samples with reads < 1,000
sample_sums(prune_samples(sample_sums(phyloseq_clean)<1000, phyloseq_clean))
length(sort(sample_sums(prune_samples(sample_sums(phyloseq_clean)<1000, phyloseq_clean)))) # 4 samples below 1,000 reads
sample_filter <- colnames(phyloseq_clean@otu_table)[sample_sums(phyloseq_clean@otu_table)<1000]
phyloseq_1000 <- subset_samples(phyloseq_clean,!sample_names(phyloseq_clean) %in% sample_filter)
phyloseq_1000 # 4090 taxa and 50 samples

### Produce final phyloseq after rarefaction
set.seed(1000)
phyloseq_rarefied <- rarefy_even_depth(phyloseq_1000, sample.size=min(sample_sums(phyloseq_1000)))
phyloseq_rarefied <- prune_taxa(taxa_sums(phyloseq_rarefied)>0, phyloseq_rarefied)
phyloseq_rarefied # 2784 taxa and 50 samples
sum(sample_sums(phyloseq_rarefied)) # 282,350 reads
phyloseq_rarefied_18S <- phyloseq_rarefied
saveRDS(phyloseq_rarefied_18S,file = "phyloseq_rarefied_18S.rds")

### Sample size per site 
sample_size <- as.data.frame(with(sample_data(phyloseq_clean), table(Site_abbreviation)))
clean_n <- as.data.frame(with(sample_data(phyloseq_rarefied), table(Site_abbreviation)))
merge(sample_size, clean_n, by="Site_abbreviation", all=TRUE)

### Pie chart of Class: read abundance across all samples
ps_class <- tax_glom(phyloseq_rarefied_18S, taxrank = "Class")
my.cols <- c(brewer.pal(8,"Dark2"), "#C46254","bisque","darkorange", rainbow_hcl(12), rainbow(13), "azure3")
data <- ps_class %>%
  psmelt() %>%
  as_tibble() %>%
  group_by(Class) %>%
  summarise(Reads = sum(Abundance), .groups = "drop") %>%
  mutate(RelAb = Reads / sum(Reads)) %>%
  mutate(Class = ifelse(RelAb < 0.01, "Other", Class)) %>% # Combine families < 1% into "Other"
  group_by(Class) %>%
  summarise(RelAb = sum(RelAb), .groups = "drop") %>%
  arrange(desc(RelAb)) %>%
  mutate(Class = factor(Class, levels = Class))
pie <- ggplot(data, aes(x = "", y = RelAb, fill = Class)) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Class", values = my.cols) +
  coord_polar("y") +
  theme_void() +
  theme(legend.title = element_text(size = 18),legend.text = element_text(size = 16),
        legend.box.spacing = unit(0, "cm"))
pie

### Closer look at taxonomic proportions (Phylum) per Site
phylumabundance <- phyloseq_rarefied_18S %>%
  tax_glom(taxrank = "Phylum") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>%
  as_tibble()
phylumabundance_avg <- phylumabundance %>%
  group_by(Phylum, Site_name, Land_use_index) %>%
  summarize(avg_abundance = mean(Abundance), .groups = "drop") %>%
  mutate(Phylum = ifelse(avg_abundance < 0.01, "Other", Phylum)) %>%
  group_by(Phylum, Site_name, Land_use_index) %>%
  summarize(avg_abundance = sum(avg_abundance), .groups = "drop")
p <- ggplot(phylumabundance_avg) +
  geom_col(mapping = aes(x=fct_reorder(Site_name, Land_use_index), y=avg_abundance, fill=fct_reorder(Phylum, avg_abundance, .desc=TRUE)), position="fill", show.legend=TRUE) +
  ylab("Proportion of 18S sequences") +
  scale_fill_manual(name="Phylum", values = my.cols) +
  xlab(NULL) + theme_bw() +
  theme(axis.text.y.left = element_text(size=15),
        axis.text.x = element_text(size=10, angle=45, hjust=1),
        axis.title.y = element_text(size=15))
p

### Diversity stats per Site
sample_data(phyloseq_rarefied_18S)$Site_name <- fct_reorder(sample_data(phyloseq_rarefied_18S)$Site_name,
                                                        sample_data(phyloseq_rarefied_18S)$Land_use_index)
plot_richness(phyloseq_rarefied_18S, x="Site_name", measures=c("Observed","Shannon"), color="Site_name") +
  geom_boxplot(alpha=0.6) + 
  scale_color_manual(values=my.cols) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=10))

### PCOA of Site
site_order <- sample_data(phyloseq_rarefied_18S) %>%
  data.frame() %>% distinct(Site_name, Land_use_index) %>%
  arrange(Land_use_index) %>% pull(Site_name)
sample_data(phyloseq_rarefied_18S)$Site_name <- factor(sample_data(phyloseq_rarefied_18S)$Site_name, levels = site_order)
ord.PCoA.bray <- ordinate(phyloseq_rarefied_18S, method="PCoA", distance="bray")
plot_ordination(phyloseq_rarefied_18S, ord.PCoA.bray, type="samples", color="Site_name") +
  geom_point(size=3, stroke=1) +
  scale_color_manual(name="Site_name", values=my.cols) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=14))
plot_ordination(phyloseq_rarefied_18S, ord.PCoA.bray, type = "samples", color = "Land_use_index") +
  geom_point(size = 3, stroke = 1) +
  scale_color_gradientn(colours = c("darkgreen", "gold", "red"), name = "Land Use Index") +
  theme_bw() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20),
    legend.title = element_text(size = 20), legend.text = element_text(size = 14))

##############################################################################
################ Next with COI ##############################################
##############################################################################

### Produce phyloseq object
ASVs <- read_qza("/Users/kbja10/Library/CloudStorage/OneDrive-ILStateUniversity/ISU/Research/Mackinaw_biodiversity/data/COI/output_29Apr2026/table.qza")
taxonomy <- read_qza("/Users/kbja10/Library/CloudStorage/OneDrive-ILStateUniversity/ISU/Research/Mackinaw_biodiversity/data/COI/output_29Apr2026/taxonomy.qza")
R_tree <- read_qza("/Users/kbja10/Library/CloudStorage/OneDrive-ILStateUniversity/ISU/Research/Mackinaw_biodiversity/data/COI/output_29Apr2026/rooted-tree.qza")
seqs <- read_qza("/Users/kbja10/Library/CloudStorage/OneDrive-ILStateUniversity/ISU/Research/Mackinaw_biodiversity/data/COI/output_29Apr2026/rep-seqs-dada2.qza")
metadata <- read_q2metadata("/Users/kbja10/Library/CloudStorage/OneDrive-ILStateUniversity/ISU/Research/Mackinaw_biodiversity/data/Mackinaw_eDNA_metadata_R.txt")
rownames(metadata) <- metadata$Sample_name_COI
names(metadata)[names(metadata) == "SampleID"] <- "Site_name"
metadata$Sample_name <- as.factor(metadata$Sample_name)
metadata$Sample_name_COI <- as.factor(metadata$Sample_name_COI)
metadata$Site_name <- as.factor(metadata$Site_name)
metadata$Site_abbreviation <- as.factor(metadata$Site_abbreviation)
metadata$Sample_type <- as.factor(metadata$Sample_type)
phyloseq_dat <- phyloseq(otu_table(ASVs$data,taxa_are_rows = TRUE), 
                         tax_table(as.matrix(parse_taxonomy(taxonomy$data))), 
                         phy_tree(R_tree$data),
                         refseq(DNAStringSet(seqs$data)),
                         sample_data(metadata))
phyloseq_dat # 7769 taxa and 70 samples
sum(sample_sums(phyloseq_dat)) # 1,783,945 reads

### Remove ASVs that are not of target length
seq_lengths <- width(refseq(phyloseq_dat))
hist(seq_lengths) # Remove below 180
phyloseq_dat <- prune_taxa(seq_lengths >= 180, phyloseq_dat)
phyloseq_dat # 7741 taxa and 70 samples

### Remove sequences with no Phylum assignment
unique(tax_table(phyloseq_dat)[, "Phylum"])
phyloseq_dat <- subset_taxa(phyloseq_dat,
                            !is.na(Phylum) & !is.na(Class))
phyloseq_dat # 4938 taxa and 70 samples

### Inspect negative controls 
phyloseq_blank <- subset_samples(phyloseq_dat, Sample_type=="blank")
phyloseq_blank <- prune_taxa(taxa_sums(phyloseq_blank) > 0, phyloseq_blank)
dim(phyloseq_blank@otu_table) # 87 ASVs present in 16 negative controls
sum(rowSums(unique(phyloseq_blank@otu_table))) # 19,953 reads total 
phyloseq_blank_table <- phyloseq_blank %>% psmelt()
p <- ggplot(phyloseq_blank_table, aes(x=Sample_name, y=Abundance, fill=Class)) + 
  geom_col(width=0.9) +
  ylab("Read Count") + xlab("Sample name") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text=element_text(size=8), axis.title=element_text(size=14),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(size=8))
p

### Identify and remove contaminants by decontam
sample_data(phyloseq_dat)$is.neg <- sample_data(phyloseq_dat)$Sample_type == "blank"
contam <- isContaminant(phyloseq_dat, method="prevalence", neg="is.neg", threshold=0.1)
contaminants <- rownames(contam)[which(contam$contaminant)] 
length(contaminants) # 4 contaminant sequences identified
phyloseq_decontam <- prune_taxa(!rownames(phyloseq_dat@otu_table) %in% contaminants, phyloseq_dat)

### Remove negative controls
phyloseq_clean <- subset_samples(phyloseq_decontam, Sample_type!="blank")
phyloseq_clean # 4934 taxa and 54 samples

### Inspect sequencing depth per sample
hist(sample_sums(phyloseq_clean@otu_table), xlab="Number of reads", main="Histogram of number of reads per sample", breaks=50)
read_depth <- data.frame(Sample_names=sample_data(phyloseq_clean)$Sample_name,
                         Site_abbreviation=sample_data(phyloseq_clean)$Site_abbreviation,
                         read_depth=sample_sums(phyloseq_clean))
p <- ggplot(read_depth, aes(x=Site_abbreviation, y=read_depth, fill=Site_abbreviation)) + 
  geom_boxplot(outlier.shape=NA, outlier.colour=NA) + 
  geom_point(shape=16, size=2, color="gray") + 
  labs(x="Sample name", y="Read Depth") + 
  theme_bw() +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        legend.position = "none")
p

### Determine rarefaction depth with ASV accumulation curves by Site
accum_curves <- accu_plot(phyloseq_clean, "Site_abbreviation", add_nb_seq=TRUE, step=1000) + 
   theme_bw()
accum_curves

### Most samples plateau by 1,000 reads -- target rarefaction depth ###

### Identify samples with reads < 1,000
sample_sums(prune_samples(sample_sums(phyloseq_clean)<1000, phyloseq_clean))
length(sort(sample_sums(prune_samples(sample_sums(phyloseq_clean)<1000, phyloseq_clean)))) # 36 samples below 1,000 reads
sample_filter <- colnames(phyloseq_clean@otu_table)[sample_sums(phyloseq_clean@otu_table)<1000]
phyloseq_1000 <- subset_samples(phyloseq_clean,!sample_names(phyloseq_clean) %in% sample_filter)
phyloseq_1000 # 4934 taxa and 18 samples

### Produce final phyloseq after rarefaction
set.seed(1000)
phyloseq_rarefied <- rarefy_even_depth(phyloseq_1000, sample.size=min(sample_sums(phyloseq_1000)))
phyloseq_rarefied <- prune_taxa(taxa_sums(phyloseq_rarefied)>0, phyloseq_rarefied)
phyloseq_rarefied # 1903 taxa and 18 samples
sum(sample_sums(phyloseq_rarefied)) # 23,544 reads
phyloseq_rarefied_COI <- phyloseq_rarefied
saveRDS(phyloseq_rarefied_COI,file = "phyloseq_rarefied_COI.rds")

### Sample size per site 
sample_size <- as.data.frame(with(sample_data(phyloseq_clean), table(Site_abbreviation)))
clean_n <- as.data.frame(with(sample_data(phyloseq_rarefied), table(Site_abbreviation)))
merge(sample_size, clean_n, by="Site_abbreviation", all=TRUE)

### Pie chart of Class: read abundance across all samples
ps_class <- tax_glom(phyloseq_rarefied_COI, taxrank = "Class")
my.cols <- c(brewer.pal(8,"Dark2"), "#C46254","bisque","darkorange", rainbow_hcl(12), rainbow(13), "azure3")
data <- ps_class %>%
  psmelt() %>%
  as_tibble() %>%
  group_by(Class) %>%
  summarise(Reads = sum(Abundance), .groups = "drop") %>%
  mutate(RelAb = Reads / sum(Reads)) %>%
  mutate(Class = ifelse(RelAb < 0.01, "Other", Class)) %>% # Combine families < 1% into "Other"
  group_by(Class) %>%
  summarise(RelAb = sum(RelAb), .groups = "drop") %>%
  arrange(desc(RelAb)) %>%
  mutate(Class = factor(Class, levels = Class))
pie <- ggplot(data, aes(x = "", y = RelAb, fill = Class)) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Class", values = my.cols) +
  coord_polar("y") +
  theme_void() +
  theme(legend.title = element_text(size = 18),legend.text = element_text(size = 16),
        legend.box.spacing = unit(0, "cm"))
pie

### Closer look at taxonomic proportions (Class) per Site
classabundance <- phyloseq_rarefied_18S %>%
  tax_glom(taxrank = "Class") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>%
  as_tibble()
classabundance_avg <- classabundance %>%
  group_by(Class, Site_name, Land_use_index) %>%
  summarize(avg_abundance = mean(Abundance), .groups = "drop") %>%
  mutate(Class = ifelse(avg_abundance < 0.01, "Other", Class)) %>%
  group_by(Class, Site_name, Land_use_index) %>%
  summarize(avg_abundance = sum(avg_abundance), .groups = "drop")
p <- ggplot(classabundance_avg) +
  geom_col(mapping = aes(x=fct_reorder(Site_name, Land_use_index), y=avg_abundance, fill=fct_reorder(Class, avg_abundance, .desc=TRUE)), position="fill", show.legend=TRUE) +
  ylab("Proportion of 18S sequences") +
  scale_fill_manual(name="Phylum", values = my.cols) +
  xlab(NULL) + theme_bw() +
  theme(axis.text.y.left = element_text(size=15),
        axis.text.x = element_text(size=10, angle=45, hjust=1),
        axis.title.y = element_text(size=15))
p

### Diversity stats per Site
sample_data(phyloseq_rarefied_COI)$Site_name <- fct_reorder(sample_data(phyloseq_rarefied_COI)$Site_name,
                                                            sample_data(phyloseq_rarefied_COI)$Land_use_index)
plot_richness(phyloseq_rarefied_COI, x="Site_name", measures=c("Observed","Shannon"), color="Site_name") +
  geom_boxplot(alpha=0.6) + 
  scale_color_manual(values=my.cols) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=10))

### PCOA of Site
site_order <- sample_data(phyloseq_rarefied_COI) %>%
  data.frame() %>% distinct(Site_name, Distance_IL_river_km) %>%
  arrange(Distance_IL_river_km) %>% pull(Site_name)
sample_data(phyloseq_rarefied_COI)$Site_name <- factor(sample_data(phyloseq_rarefied_COI)$Site_name, levels = site_order)
ord.PCoA.bray <- ordinate(phyloseq_rarefied_COI, method="PCoA", distance="bray")
plot_ordination(phyloseq_rarefied_COI, ord.PCoA.bray, type="samples", color="Site_name") +
  geom_point(size=3, stroke=1) +
  scale_color_manual(name="Site_name", values=my.cols) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=14))
plot_ordination(phyloseq_rarefied_COI, ord.PCoA.bray, type = "samples", color = "Land_use_index") +
  geom_point(size = 3, stroke = 1) +
  scale_color_gradientn(colours = c("darkgreen", "gold", "red"), name = "Land Use Index") +
  theme_bw() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20),
        legend.title = element_text(size = 20), legend.text = element_text(size = 14))

#######################################################
##### Shannon diversity across sites for all three loci
#######################################################

# Extract Shannon and collapse to site
get_site_shannon <- function(ps, locus_name) {
  shannon <- estimate_richness(ps, measures="Shannon")
  meta <- sample_data(ps) %>%
    data.frame()
  df <- cbind(meta, shannon) %>%
    group_by(Site_name, Land_use_index, Distance_IL_river_km) %>%
    summarise(Shannon = mean(Shannon, na.rm = TRUE), .groups = "drop") %>%
    mutate(Locus=locus_name, Land_use_index=as.numeric(as.character(Land_use_index)))
  return(df)
}

# Combine all loci at site level
shannon_sites <- bind_rows(
  get_site_shannon(phyloseq_rarefied_12S, "12S"),
  get_site_shannon(phyloseq_rarefied_COI, "COI"),
  get_site_shannon(phyloseq_rarefied_18S, "18S"))

# Plot
locus.cols <- c("12S"="dodgerblue4", "COI"="darkorange2", "18S"="forestgreen")

# Set facet order
shannon_sites$Locus <- factor(shannon_sites$Locus, levels = c("12S", "COI", "18S"))
p <- ggplot(shannon_sites, aes(x=Land_use_index, y=Shannon, color=Locus)) +
  geom_point(aes(shape=Locus), size=3) +
  geom_smooth(method="lm", se=FALSE, linewidth=1.2) +
  facet_wrap(~ Locus, nrow=1) +
  scale_color_manual(values=locus.cols) +
  theme_bw() +
  labs(x="Land use intensity", y="Site-level Shannon diversity") +
  theme(axis.text=element_text(size=14),
    axis.title=element_text(size=18),
    strip.text=element_text(size=18, face="bold"),
    legend.position="none")
p
