###Load Packages
library(phyloseq)
library(ggplot2)
library(scales)
library(grid)
library(vegan)
library(RColorBrewer)
library(DESeq2)
library(reshape2)
library(dplyr)
library(rio)
library(data.table)
library(VennDiagram)
library(microbiome)
library(goeveg)
library(venneuler)

##DADA2
load(file="Namrata_DADA2.RData")
samdf <- read.table("Namrata_mappingfile.txt")
colnames(samdf)<- c("SampleID","MouseID","Organ","Location","Genotype")
mappingfile <- samdf[,-1]
rownames(mappingfile) <- samdf[,1]
physeq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(mappingfile), 
                   tax_table(taxa.plus))
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
physeq1 = merge_phyloseq(physeq, mappingfile, random_tree)
save.image(file="Namrata_physeq1.RData")

###Load Data
load(file="Namrata_physeq1.RData")
theme_set(theme_bw())

###Remove Singletons
physeq1_trim<-prune_taxa(taxa_sums(physeq1) > 1, physeq1)
physeq1_trim

###Subset Data
Fecal<- physeq1_trim %>%
  subset_samples(Location == "Fecal")
Fecal

##PCoA Fecal Weighted Unifrac PCoA
Analysis<-Fecal
GP1=transform_sample_counts(Analysis, function(x) 1E6 * x/sum(x))
orduW = ordinate(GP1, "PCoA","unifrac",weighted=TRUE)
pW = plot_ordination(GP1, orduW, color = "Genotype", title = "Weighted Unifrac")+ geom_point(size = 2, shape = 1, colour = "black", stroke = 1)
pW + theme(plot.title = element_text(size=18))+ theme(legend.text=element_text(size=14)) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + scale_colour_manual(values=c("red","white"))

###Permanova by Genotype
set.seed(42)
analysis_unifrac_weighted<-phyloseq::distance(Analysis,method = "unifrac", weighted=TRUE)
sampledf<- data.frame(sample_data(Analysis))
adonis(analysis_unifrac_weighted ~ Genotype, data = sampledf)

##PCoA Fecal Unweighted Unifrac PCoA
GP1=transform_sample_counts(Analysis, function(x) 1E6 * x/sum(x))
orduU = ordinate(GP1, "PCoA","unifrac",weighted=FALSE)
pU = plot_ordination(GP1, orduU, color = "Genotype", title = "Unweighted Unifrac")+ geom_point(size = 2, shape = 1, colour = "black", stroke = 1)
pU + theme(plot.title = element_text(size=18))+ theme(legend.text=element_text(size=14)) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + scale_colour_manual(values=c("red","white"))

###Permanova by Genotype
set.seed(42)
analysis_unifrac_unweighted<-phyloseq::distance(Analysis,method = "unifrac", weighted=FALSE)
sampledf<- data.frame(sample_data(Analysis))
adonis(analysis_unifrac_unweighted ~ Genotype, data = sampledf)

##Normalize dataset
subset=transform_sample_counts(Analysis,function(x)x/sum(x))

###Condense by taxonomic level
ps.class<-tax_glom(Analysis,taxrank="Class")

###Percent Abundance
otu.table<-otu_table(ps.class)
write.csv(otu.table,"class.otu.table.csv")
tax.table<-tax_table(ps.class)
write.csv(tax.table,"class.tax.table.csv")

