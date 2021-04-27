
###############################################################
   ## Sweet spheres -- Bunse et al. ##
###############################################################

## This script: ASV-based community composition ##

# load packages 
library(phyloseq)
library(ampvis2)
library(vegan)
library(tibble)
library(PMCMRplus)
library(ggplot2)
library(cowplot)
library(dplyr)


###############################################################
   ###  ORDINATION -- Fig. 1a ###
###############################################################

amp_ordinate(
  ampvis, 
  type = "NMDS",
  distmeasure = "bray",
  transform = "hellinger", 
  sample_color_by = "condition", sample_colorframe = T,
  sample_shape_by = "time_h",
  sample_point_size = 3) +
scale_color_manual(values=col.treatment) +
scale_fill_manual(values=col.treatment) +
scale_shape_manual(values=c(18,15,17,19)) + 
theme_bw() + 
theme(axis.text = element_blank(),
      #axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right")


###################################################################################
   ###  SIGNIFICANT DISTANCES? -- PERMANOVA ###
###################################################################################

#  PA vs FL ***
dist.sample <- subset_samples(
  pseq.rel, sample %in% c("Particle","Filter"))
otu_table(dist.sample) = otu_table(decostand(
  otu_table(dist.sample), method="hellinger"), 
  taxa_are_rows=T)
dist.sample <- phyloseq::distance(
  dist.sample, method = "bray")

adonis2(
  dist.sample ~ sample, 
  data = subset(ENV, sample %in% c(
    "Particle","Filter")), 
  sqrt.dist = T)

#########################

# PA vs Control ***
dist.sample <- subset_samples(
  pseq.rel, sample %in% c("Particle","Control"))
otu_table(dist.sample) = otu_table(decostand(
  otu_table(dist.sample), method="hellinger"), 
  taxa_are_rows=T)
dist.sample <- phyloseq::distance(
  dist.sample, method = "bray")

adonis2(
  dist.sample ~ sample, 
  data = subset(ENV, sample %in% c(
    "Particle","Control")), 
  sqrt.dist = T)

#########################

# FL vs Control ***
dist.sample <- subset_samples(
  pseq.rel, sample %in% c("Filter","Control"))
otu_table(dist.sample) = otu_table(decostand(
  otu_table(dist.sample), method="hellinger"), 
  taxa_are_rows=T)
dist.sample <- phyloseq::distance(
  dist.sample, method = "bray")

adonis2(
  dist.sample ~ sample, 
  data = subset(ENV, sample %in% c(
    "Filter","Control")), 
  sqrt.dist = T)

#########################

# Rreplicate congruency - n.s. as expected (=identical)
dist.sample <- subset_samples(
  pseq.rel, condition %in% c("PecP"))
otu_table(dist.sample) = otu_table(decostand(
  otu_table(dist.sample), method="hellinger"), 
  taxa_are_rows=T)
dist.sample <- phyloseq::distance(
  dist.sample, method = "bray")

adonis2(
  dist.sample ~ replicate, 
  data = subset(ENV, condition %in% c(
    "PecP")), 
  sqrt.dist = T)


############################################################################################
   ###  TOP GENERA -- Fig. 1c ###
############################################################################################

amp_heatmap(
  data = ampvis,
  tax_aggregate = "Genus",
  tax_show = c(
    "Colwellia","Tenacibaculum","Psychrobium",
    "Catenovulum","Psychromonas",
     "Pseudoalteromonas","Glaciecola"),
  normalise = T,
  plot_values = T,
  round = 0,
  max_abundance = 26,
  min_abundance = 1,
  plot_colorscale = "sqrt",
  plot_legendbreaks = c(5, 10, 20),
  facet_by = "condition",
  group_by = "time_h",
  color_vector = c(
    "white","#AEDBE7","#93CFE0",
    "#43ABC9","#1287a8","#3c6478",
    "#0c374d","#093145"))

############################################################################################
   ###  GENUS + DIVERSITY DETAILS -- Fig. S2  ###
############################################################################################
 
ASV.avg <- psmelt(genus.rel) %>% 
  as_tibble %>% 
  group_by(
    Genus, condition, time_h, 
    sample, timepoint) %>%
  summarize_at(c("Abundance"), mean) %>%
  ungroup

# Reorder Abundances
ASV.avg$condition <- factor(
ASV.avg$condition, levels = c(
  "start","AlgP","PecP","FL","CTR"))

# Reorder AlphaDiv 
Diversity$condition <- factor(
Diversity$condition, levels = c(
  "start","AlgP","PecP","FL","CTR"))

# Plot
plot1 <- ggplot(data = subset(
    ASV.avg, Abundance > 4), 
    aes(x=time_h, y=Abundance, fill=Genus)) + 
geom_bar(aes(), stat = "identity", position = "stack") +
scale_fill_manual(values=col.taxon) +
scale_y_continuous(
  limits = c(0, 90), 
  expand = c(0, 0)) +
facet_grid(~condition, scales="free") +
theme_bw() +
theme(legend.position="top",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()) +
guides(fill = guide_legend(ncol=5))

plot2 <- ggplot(
  Diversity, aes(
    x=time_h, y=Simpson, fill=condition)) +
geom_boxplot() + 
scale_fill_manual(values=col.treatment) +
facet_grid(~condition, scales = "free") +
theme_bw() +
theme(legend.position="none",
      panel.grid.minor = element_blank())

plot_grid(
  plot1,
  plot2,
  ncol=1,
  rel_heights = c(2,1))


###################################################################################
   ###  STATS: TREATMENT + TEMPORAL PATTERNS  ###
###################################################################################

# Melt + avg  
ASV.all <- psmelt(genus.rel) %>% 
  as_tibble %>% 
  group_by(
    OTU, Genus, condition, time_h, 
    sample, timepoint) %>%
  ungroup

# significant enrichment overall?
kruskal.test(
  Abundance ~ sample, 
  data = subset(ASV.all, Genus=="Tenacibaculum"))
posthoc.kruskal.dunn.test(
  Abundance ~ sample, 
  data = subset(ASV.all, Genus=="Tenacibaculum"), 
  p.adjust.method = "bonf")

# Aureispira at 60h***
kruskal.test(
  Abundance ~ time_h, 
  data = subset(ASV.all, 
  Genus=="Aureispira" &
  time_h %in% c("24","60") &
  sample %in% c("Particle","Filter")))

# Catenovulum 24 vs 60h**
kruskal.test(
  Abundance ~ time_h, 
  data = subset(ASV.all, 
  Genus=="Catenovulum" &
  time_h %in% c("24","48","60") &
  condition %in% c("PecP","AlgP")))
posthoc.kruskal.dunn.test(
  Abundance ~ time_h, 
  data = subset(
    ASV.all, Genus=="Catenovulum" &
    time_h %in% c("24","60") &
    sample %in% c("Particle")))

# Tenaci -- PecP 24 vs 60*
kruskal.test(
  Abundance ~ time_h, 
  data = subset(
    ASV.all, Genus=="Pseudoalteromonas" &
    time_h %in% c("24","48","60") &
    condition %in% c("PecP")))
posthoc.kruskal.dunn.test(
  Abundance ~ time_h, 
  data = subset(ASV.all, 
                Genus=="Pseudoalteromonas" &
                  time_h %in% c("24","48","60") &
                  condition %in% c("PecP")))

# Ps.alt -- PA 24 vs 60** 
kruskal.test(
  Abundance ~ time_h, 
  data = subset(
    ASV.all, 
    Genus=="Pseudoalteromonas" &
    time_h %in% c("24","60") &
    sample %in% c("Particle")))
posthoc.kruskal.dunn.test(
  Abundance ~ time_h, 
  data = subset(ASV.all, 
  Genus=="Pseudoalteromonas" &
  time_h %in% c("24","48","60") &
  sample %in% c("Particle")))


###################################################################################
   ###  MICRODIVERSITY -- Fig. S3 ###
###################################################################################

## subset taxa
asv_col <- subset_taxa(pseq.rel, Genus == "Colwellia")
asv_col = filter_taxa(asv_col, function(x) sum(x)>2, T)
asv_col = as(otu_table(asv_col), "matrix") 
asv_col <- reshape2::melt(asv_col) %>%
  add_column(tax="Colwellia")

asv_glac <- subset_taxa(
  pseq.rel, Genus == "Glaciecola")
asv_glac = filter_taxa(asv_glac, function(x) sum(x)>3, T)
asv_glac = as(otu_table(asv_glac), "matrix") 
asv_glac <- reshape2::melt(asv_glac) %>%
  add_column(tax="Glaciecola")

## Merge data
microdiv <- rbind(
  asv_col, asv_glac) %>% 
  dplyr::rename(SampleID = Var2, ASV = Var1) %>%
  left_join(dplyr::select(
    ENV, SampleID, condition, time_h, timepoint), 
    by=c("SampleID"="SampleID")) %>%
  group_by(ASV, time_h, condition, tax) %>%
  summarize_at(c("value"), mean, na.rm=T) %>% 
  ungroup %>%
  mutate(condition=factor(condition, levels=c(
    "start","AlgP","PecP","FL","CTR")))

## Plot 
ggplot(data=subset(
  microdiv, value > 0.41),
  aes(x = as.factor(time_h), y = value)) +
  geom_bar(aes(
    group=ASV, fill=ASV), 
    position ="dodge", stat="identity")+ 
  scale_fill_manual(values = c(col.25)) +
  facet_grid(tax~condition, scales="free") +
  guides(fill = guide_legend(ncol=12))+
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    axis.title = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom")

# save ASV data if desired
pivot_wider(
    data = microdiv,
    id_cols = c("ASV","tax"), 
    names_from = c("condition","time_h"),
    values_from = c("value"),
    names_sort = T) %>%
  mutate_if(is.numeric, round, 2) %>%
write.table(., "microdiv.txt", sep="\t")            

###################################################################################

save.image("BAH_particles.Rdata")

