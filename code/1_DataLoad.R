
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
   ###  SET COLORS FOR PLOTTING  ###
############################################################################################

col.taxon <- c(
  "Alphaproteobacteria"="#4769d1",
  "Gammaproteobacteria"="#77d6dd",
  "Deltaproteobacteria"="#f2de4f",
  "Bacteroidia"= "#c93636", 
  "Planctomycetacia"="#965ec4", 
  "Verrucomicrobiae"= "#3d5e34",
  "NS9_marine_group_unclassified" = "#dec50b",
  "NS2b_marine_group" = "#ff7a7a", 
  "NS4_marine_group" = "#82eff5", 
  "NS5_marine_group" = "#82eff5",
  "Cryomorphaceae_unclassified" = "#000000", 
  "Aurantivirga" = "#2e7ec9", 
  "Formosa" = "#117777", 
  "Polaribacter_1" = "#22c952", 
  "Polaribacter_3" = "#8a6ca3", 
  "Amylibacter" = "#e5e6a1", 
  "Glaciecola" = "#419cba",
  "Colwellia" = "#64d974",
  "Psychrobium" = "#ffb463",
  "Psychromonas" = "#f5e90c",
  "Tenacibaculum" = "#ab5572",
  "Thalassotalea" = "#f5e5c6",  
  "Wenyingzhuangia" = "#2f6e32",  
  "SAR92_clade" = "#070087", 
  "Acinetobacter" = "#626f99",
  "Persicirhabdus" = "#c9bc93",
  "Rubritalea" = "#e7dfe8", 
  "Planktomarina" = "#ed7dff", 
  "Pseudoalteromonas" = "#f73982",
  "Ulvibacter"= "#a1c2bd",
  "Flavicella" = "#ff8000",
  "Fluviicola" = "#ff7a7a",
  "Lentimonas" = "#2e7ec9")

col.treatment = c(
  "CTR" = "#5f6161", 
  "ALG" = "#d9582e",
  "PEC" = "#EBCC2A",
  "AlgP" = "#d9582e", 
  "PecP" = "#EBCC2A",
  "FL" = "#78B7C5",
  "T0" = "#7900f2",
  "start" = "#7900f2")
   
############################################################################################

# remove temporary datasets
rm(asv, tax, ampvis.env)
save.image("BAH_particles.Rdata")
