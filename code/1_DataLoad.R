
############################################################################################
   ##  Sweet spheres -- Bunse et al.  ##
############################################################################################

## This script: load and format DADA2 ASVs & metadata
## taxdata: NA/uncultured/unclassified replaced w/ last tax-level + "unclassified" 

# Load packages and colors
library(gtools)
library(phyloseq)
library(ampvis2)
library(dplyr)
library(reshape2)
source("6_Colors.R")


############################################################################################
   ###  RAW DATA -- LOAD AND FORMATTING ###
############################################################################################

# Load ASVs
ASV <- read.table(
  "ASV_seqtab.txt",
  h = T,
  sep = "\t")

# Load taxonomy 
TAX <- read.table(
  "ASV_tax.txt",
  h = T, 
  sep = "\t", 
  stringsAsFactors = F, 
  row.names = 1)

##############################################

# Import metadata 
ENV <- read.table("metadata.txt", 
  h=T, sep = "\t", stringsAsFactors=F)
row.names(ENV) = ENV$SampleID

# Generate factors for parameter of interest
ENV$condition <- factor(ENV$condition, 
  levels = c("start","CTR","FL","AlgP","PecP"))
ENV$sample <- factor(ENV$sample, 
  levels = c("Start","Control","Filter","Particle"))
ENV$time_h <- factor(ENV$time_h, 
   levels = c("0","24","48","60"))
ENV$timepoint <- factor(ENV$timepoint, 
   levels = c("start","early","late"))

###############################################

## FORMAT TABLES

# rename BAC 
ASV <- ASV[,mixedsort(names(ASV))]
ENV <- ENV[mixedorder(ENV$seq_clip),]
colnames(ASV) = ENV$SampleID

# Check: reorder so ENV and ASV tables match
ENV <- ENV[mixedorder(ENV$SampleID),]
ASV <- ASV[, rownames(ENV)]


############################################################################################
   ###  PHYLOSEQ LOAD  ###
############################################################################################

asv = otu_table(ASV, taxa_are_rows=T)
tax = tax_table(TAX)
rownames(tax) <- rownames(asv)

pseq <- phyloseq(
  otu_table(asv, taxa_are_rows=F), 
   sample_data(ENV), 
    tax_table(tax))

# Fix tax-IDs
colnames(pseq@tax_table)<- c(
  "Kingdom","Phylum","Class","Order","Family","Genus")

####################################

## Filter: ASV >3 counts in >3% of samples; to rel. ab
pseq.abs = filter_taxa(
  pseq, function(x) sum(x > 3) > (0.03 * length(x)), T)
pseq.rel = transform_sample_counts(
  pseq.abs, function(x) x / sum(x) * 100) 

## Agglomerate taxranks
genus.abs <- tax_glom(
  pseq.abs, taxrank="Genus")
genus.rel <- tax_glom(
  pseq.rel, taxrank="Genus")


############################################################################################
   ###  AMPVIS LOAD  ###
############################################################################################

ampvis <- data.frame(OTU = rownames(
   phyloseq::otu_table(pseq.abs)@.Data),
   phyloseq::otu_table(pseq.abs)@.Data,
   phyloseq::tax_table(pseq.abs)@.Data,
   check.names=F)

# Add dummy "Species" column (required >> copy "Genus")
ampvis$Species <- ampvis$Genus

# Extract metadata from phyloseq; format and combine
ampvis.env <- data.frame(
  phyloseq::sample_data(pseq.abs), check.names=F)
ampvis.env <- cbind(
  Sample = rownames(ampvis.env), ampvis.env) 
ampvis <- amp_load(ampvis, ampvis.env)


############################################################################################

# remove temporary datasets
rm(asv, tax, ampvis.env)
save.image("BAH_particles.Rdata")
