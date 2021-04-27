############################################################################################
   ##  Sweet spheres - Bunse et al  ##
############################################################################################

# This script: KEGG / dbCan analysis MAGs

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(PNWColors)
library(PMCMR)


############################################################################################

## CAZYmes --- dbCan ##
# filter for MAGs of interest

bins_dbcan <- meta_dbcan %>%
  filter(auto_bin %in% c(
    "21","24","26","32","57","73")) %>%
  dplyr::select(-c(
    "Gene_start","Gene_length","Strand",
    "Contig","Contig_length","Gene",
    "EC_number","COG"))

# sort naturally
bins_dbcan <- bins_dbcan[
  mixedorder(bins_dbcan$locus_tag), ]

# Export table 
write.table(bins_dbcan, file="MAGs_dbcan.txt", 
    sep ="\t", row.names=F)

############################################################################################

## CAZymes - MAG21

# Append highly expressed GHs with >70% coverage
# rename NA to GH for simplicity
bin21 <- meta_exp %>%
  filter(locus_tag %in% c(
    "PROKKA_108763","PROKKA_108765")) %>%
  bind_rows(meta_dbcan) %>%  
  replace(is.na(.), "GH")

# calculate mean per treatment  
bin21 <- bin21 %>%
  mutate(AlgP = rowMeans(
    dplyr::select(., starts_with("Alg")))) %>%
  mutate(PecP = rowMeans(
    dplyr::select(., starts_with("Pec")))) %>%
  mutate(FL = rowMeans(
    dplyr::select(., starts_with("FL"))))

# reshape
bin21 <- bin21 %>%
  filter(auto_bin=="21") %>%
  mutate(CAZyme = str_extract(HMM.Profile,"^.{2}")) %>%
  #filter(locus_tag !="PROKKA_118003") %>%
  reshape2::melt(measure.vars=c("AlgP","PecP","FL"))

# order for plotting
bin21$CAZyme <- factor(bin21$CAZyme, levels=c(
    "GH","GT","CE","PL","CB"))
bin21$variable <- factor(bin21$variable, levels=c(
  "PecP","AlgP","FL"))

ggplot(bin21, 
  aes(x=CAZyme, y=value, fill=CAZyme)) +
geom_bar(
  stat="identity", 
  position=position_dodge(width=0.9)) +
scale_x_discrete(
  expand = c(0.25,0)) +  
scale_y_continuous(
  limits = c(0,130),
  expand = c(0.00,0.1),
  breaks = c(0,50,100)) + 
scale_fill_manual(values = rev(pnw_palette("Bay",5))) +
facet_wrap(variable~., ncol=1) +
theme_bw() +
theme(panel.grid.minor = element_blank()) 
        
kruskal.test(
  value ~ variable, 
  data = subset(bin21, CAZyme=="GH" &
    variable %in% c("AlgP","PecP")))
posthoc.kruskal.dunn.test(
  value ~ variable, 
  data = subset(
    bin21, CAZyme=="GH", 
  p.adjust.method = "bonf"))

############################################################################################

# Is diff. expression in PecP due to MAG21?
# Inout files: see https://github.com/matthiaswietz/sweet-spheres/tree/main/data

# load DESEQ data
deseq_AlgPec <- read.csv(
  "./meta-omics/DESEq_AlgP-PecP.csv", h=T, sep=",")

# load locus tags 
bins_loci <- read.csv(
  "./MAGs/bins_locus-tags.txt", h=T, sep="\t")

# reformat KEGG data
bins_KEGG <- meta_KEGG %>%
  dplyr::select(c("KEGG","locus_tag","Category2","Category3"))

# merge everything
bins_deseq <- merge(merge(
  deseq_AlgPec, bins_loci, by="locus_tag"),
  bins_KEGG, by="locus_tag", all=T) %>%
filter(auto_bin %in% c("21","24","26","57","73")) %>%
distinct(KEGG, locus_tag, .keep_all=T) %>%
mutate(condition = case_when(
  log2FoldChange > 2 ~ "PEC", 
  log2FoldChange < -2 ~ "ALG")) 

# count by MAG
count <- plyr::count(bins_deseq, c(
  "auto_bin","condition")) %>%
  tidyr::drop_na()

ggplot(count, 
  aes(x=factor(auto_bin), y=condition, 
      size=freq, group=1)) +
geom_point(color="gray33") +   
scale_size(range = c(1,12)) +
theme_bw()+
theme(axis.text = element_text(size=10)) 

# remove temp-files
rm(count)


############################################################################################

save.image("BAH_particles.Rdata")