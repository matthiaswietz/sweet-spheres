############################################################################################
   ## Sweet spheres -- Bunse et al.  ##
############################################################################################

# This script: KEGG / dbCan analysis metagenomes

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)


############################################################################################

## Load genes, translations and expression ##
# These large inputs are stored at https://zenodo.org/record/4171148

# Translations
meta_faa <- read.csv(
  "./meta-omics/MetaG_translations.txt", 
  h=T, sep="\t", stringsAsFactors=F)

# Gene/contig info
meta_info <- read.csv(
  "./meta-omics/MetaG_genes.txt", 
  h=T, sep="\t", stringsAsFactors=F)

# Expression
meta_exp <- read.csv2(
  "./meta-omics/MetaT_TPM.csv", 
  h=T, sep=";", row.names=1, stringsAsFactors=F)

# rename 
names(meta_exp) = c(
  "locus_tag","FL-1","FL-2","FL-3",
  "AlgP-1","AlgP-2","AlgP-3",
  "FL-4","FL-5","FL-6",
  "PecP-1","PecP-2","PecP-3")

# Merge with contig info
meta_exp = inner_join(
  meta_exp, meta_info, by="locus_tag", all.x=T) %>%
  left_join(meta_faa, by='locus_tag') %>%
  mutate_if(is.numeric, round, 3)

# sort naturally
meta_exp <- meta_exp[
  mixedorder(meta_exp$locus_tag), ]

# Export
write.table(meta_exp, file = "MetaOmics_Exp.txt", 
  sep ="\t", row.names=F)

############################################################################################

## KEGG annotation ##

# based on GhostKoala on MetaG_genes.faa (on Zenodo)
# output saved as KEGG_anno.txt
# Add locus_tag/KEGG header  

# KOs  "human diseases", "organismal systems", "aging" removed
# All following inputs are on https://github.com/matthiaswietz/sweet-spheres/tree/main/data

# Import results
KEGG.anno <- read.csv(
  "./meta-omics/KEGG_anno.txt", h=T, sep="\t")

# Import KO-Table
KEGG.ortho <- read.csv(
  "./meta-omics/KEGG_ortho.txt", h=T, sep="\t")

# Merge 
meta_KEGG = merge(
  KEGG.anno, KEGG.ortho, by="KEGG", all.x=T) %>%
  tidyr::drop_na() %>%
  left_join(meta_exp, by="locus_tag") %>%
  mutate_if(is.numeric, round, 3) %>%
  dplyr::select(-c(
    "Category1","Contig_length","Gene",
    "Gene_start","Gene_length","Strand",
    "COG","EC_number","product"))

# sort naturally
meta_KEGG <- meta_KEGG[
  mixedorder(meta_KEGG$locus_tag), ]

# Export full table for supplement
write.table(
  meta_KEGG, file = "MetaOmics_KEGG.txt", 
   sep ="\t", row.names=F)


############################################################################################

## CAZYme annotation ##

# based on MetaG_genes.faa (on Zenodo)
# separated into five chunks (on linux-server: split -d -b 18M)
# Sequence headers modified so only PROKKA ID remains
# Submitted to dbCan2 
# Downloaded HMMER output; keep only evalue < 1e-15, coverage > 0.80

# This file dbcan_coverage80.csv is available on Github
# further processed as follows:

# Load  
meta_dbcan <- read.csv(
  "./meta-omics/dbcan_coverage80.csv", h=T, sep=";")

# Remove everything except CAZY-ID 
meta_dbcan[] <- lapply(meta_dbcan, 
  gsub, pattern="_Glycos_transf_[1-9].hmm", replacement="")   
meta_dbcan[] <- lapply(meta_dbcan, 
  gsub, pattern="_Glyco_tranf_[1-9].hmm", replacement="") 
meta_dbcan[] <- lapply(meta_dbcan, 
  gsub, pattern="_Glyco_tranf_[1-9]_[1-9].hmm", replacement="") 
meta_dbcan[] <- lapply(meta_dbcan, 
  gsub, pattern="_Glyco_trans_[1-9]_[1-9]", replacement="") 
meta_dbcan[] <- lapply(meta_dbcan, 
  gsub, pattern=".hmm", replacement ="") 

# Merge with expression
meta_dbcan <- merge(
  meta_exp, meta_dbcan, by="locus_tag", all.x=F) %>%
  mutate_if(is.numeric, round, 3) %>%
  dplyr::select(-c("HMM_Profile_simple","HK"))

# sort naturally
meta_dbcan <- meta_dbcan[
  mixedorder(meta_dbcan$locus_tag), ]

# Export full table for supplement
write.table(
  meta_dbcan, file = "MetaOmics_dbcan.txt", 
  sep ="\t", row.names=F)

######################################################

## Top CAZymes for SignalP analysis

# Calculate mean by Alg / Pec / FL
dbcan_avg <- meta_dbcan %>%
  mutate(AlgP = rowMeans(
    dplyr::select(., starts_with("Alg")), na.rm=T)) %>%
  mutate(PecP = rowMeans(
    dplyr::select(., starts_with("Pec")), na.rm=T)) %>%
  mutate(FL = rowMeans(
    dplyr::select(., starts_with("FL")), na.rm=T))

# filter + export top X in PA
dbcan_avg %>%
  filter(AlgP & PecP > 1 & FL < 1) %>%
write.table(., file="MetaOmics_dbcan_topPA.txt", 
  sep ="\t", row.names=F)

# filter + export top X in FL
dbcan_avg %>%
  filter(AlgP & PecP < 1 & FL > 1) %>%
write.table(., file="MetaOmics_dbcan_topFL.txt", 
  sep ="\t", row.names=F)

############################################################################################

## STATS ##
## Higher expression of PLs on AlgP?

dbcan_stats <- meta_dbcan %>%
  mutate_at(
    .vars = "HMM.Profile", 
    .funs = gsub,
    pattern = "_[0-9]",
   replacement = "") %>% 
mutate(AlgP = rowMeans(dplyr::select(
  ., starts_with("Alg")), na.rm=T)) %>%
mutate(PecP = rowMeans(dplyr::select(
  ., starts_with("Pec")), na.rm=T)) %>%
mutate(FL = rowMeans(dplyr::select(
  ., starts_with("FL")), na.rm=T)) %>%
  dplyr::select(c("HMM.Profile","AlgP","PecP")) %>%
filter(HMM.Profile=="PL18")  %>%
reshape2::melt(id.vars=c("HMM.Profile")) 

wilcox.test(
  value ~ variable, 
  data = dbcan_stats, paired=T, exact=F)

# PL6 PL7 PL17 PL18 all significant; PL15 not

######################################################

## Gln-Synthase + NH4 transport higher in PA?
## Methylcitrate cycle and glyox shunt higher in FL?

AmmTrans <- read.csv(
  "./meta-omics/Ammonium.txt", 
  h=T, sep="\t", stringsAsFactors=F) %>%
  dplyr::select(c("locus_tag"))  %>%
  left_join(meta_exp, by="locus_tag")

# Add P-II N regulator (glnB, glnK) 
# Glutamin synthetase (6.3.1.2 -- glnA, glnII, glnT)
# Glutamate synthase (1.4.1.13)
# Glutamate dehydrogenase (1.4.1.2, 1.4.1.3)
# Methylcitrate (acn=COG1048, prpB=COG2513, prpC=COG0372, prpF=COG2828)
# Glyox cycle (4.1.3.1, 2.3.3.9)
AmmMeta <- meta_exp %>%
  filter(EC_number %in% c(
    "2.3.3.9","4.1.3.1",
    "6.3.1.2","1.4.1.13",
    "1.4.1.2","1.4.1.3") |
  COG %in% c("COG2828","COG0372","COG2513",
    "COG2224","COG2225","COG1048") |
  grepl("glnB|glnK", Gene)) 

# Reformat for plotting
AmmExpr <- bind_rows(AmmTrans, AmmMeta) %>%
  mutate(product = case_when(
    grepl("Ammonia", product) ~ "amtB", 
    grepl("Ammonium", product) ~ "nrgA",
    TRUE ~ product)) %>% 
  mutate(fct = case_when(
    EC_number=="6.3.1.2" ~ "GlnSynth",
    EC_number=="1.4.1.13" ~ "GlcSynth",
    EC_number %in% c("1.4.1.2", "1.4.1.3") ~ "GlcDehyd",
    #product=="amtB" ~ "amtB",
    product %in% c(
      "nrgA","amtB","hypothetical protein") ~ "AmmTrans",
    product %in% c("Nitrogen regulatory protein P-II",  
      "Nitrogen regulatory protein P-II 2") ~ "P-II",
    COG %in% c(
      "COG2828","COG0372",
      "COG2513","COG1048") ~ "Citrate",
    COG %in% c("COG2224","COG2225") ~ "Glyox")) %>%
  mutate(type = case_when(
    fct %in% c(
      "GlnSynth","GlcSynth","GlcDehyd",
      "AmmTrans","P-II") ~ "Ammonium",
    fct %in% c(
      "Citrate","Glyox") ~ "Fatty acids")) %>%
  mutate(PA = rowMeans(dplyr::select(
    ., starts_with("Alg"), starts_with("Pec")))) %>%
  mutate(FL = rowMeans(dplyr::select(
    ., starts_with("FL")))) %>%
  dplyr::select(c(
    "PA","FL","fct","type",
    "locus_tag","product")) %>%
  drop_na() %>%
  reshape2::melt() 

# Reorder for plotting
AmmExpr$fct <- factor(AmmExpr$fct, levels=c(
  "GlnSynth","GlcSynth","GlcDehyd",
  "P-II","AmmTrans","Citrate","Glyox")) 

# Plot without amtB (same in PA+FL)
ggplot(data=subset(AmmExpr, !fct %in% 
   c("GlcSynth","GlcDehyd") 
   #& product!="amtB"
   )) +  
  aes(x=variable, y=value, fill=fct) +
geom_bar(stat="identity") +
scale_fill_viridis(
 option = "D",
  begin=0.99, end=0.12,
  discrete=T) +
facet_wrap(~type, ncol=1, scales="free") +
theme_bw() 

# nrgA*   P-II* Citrate**
wilcox.test(
  value ~ variable, 
  data = subset(AmmExpr, fct %in% c(
    "nrgA")), paired=T)

# GlcDehyd + GlnSynth ***
wilcox.test(
  value ~ variable, 
  data = subset(AmmExpr, fct %in% c(
    "GlnSynth")), paired=T)

######################################################

# Val/Ile/Leu differences PA-FL

stats <- meta_KEGG %>%
  mutate_if(is.factor,as.character) %>%
  filter(grepl("00280|00290", Category3)) %>% 
  mutate(fct = case_when(
    grepl("00280", Category3) ~ "Degrade", 
    grepl("00290", Category3) ~ "Biosynth",
    TRUE ~ Category3)) %>% 
  mutate(PA = rowMeans(
    dplyr::select(
      ., starts_with("Alg"),starts_with("Pec")))) %>%
  mutate(FL = rowMeans(
    dplyr::select(., starts_with("FL")))) %>%
  dplyr::select(c("fct","PA","FL")) %>%
  reshape2::melt(id.vars=c("fct")) 

# Biosynth** Degrade********
wilcox.test(
  value ~ variable, 
  data = subset(stats, fct %in% c(
    "Biosynth")), paired=T)

############################################################################################

# remove temporary datasets
rm(KEGG.anno, KEGG.ortho, meta_info, AmmMeta, AmmTrans)

save.image("BAH_particles.Rdata") 
