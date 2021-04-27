############################################################################################
   ##  Sweet spheres -- Bunse et al. ##
############################################################################################

## This script: rarefaction and diversity analyses ##

library(iNEXT)
library(olsrr)
library(cowplot)
library(phyloseq)

############################################################################################

# Load raw counts
iNEXT = otu_table(
  ASV, taxa_are_rows = F)

# Calculate results
iNEXT = iNEXT(
  as.data.frame(otu_table(iNEXT)), q=c(0),
  datatype="abundance", conf = 0.95, nboot = 100)

# Rarefaction
rarefac <- fortify(iNEXT, type=1)
rarefac.point <- rarefac[
  which(rarefac$method == "observed"),]
rarefac.line <- rarefac[
  which(rarefac$method != "observed"),]
rarefac.line$method <- factor(rarefac.line$method,
     c("interpolated", "extrapolated"),
     c("interpolation", "extrapolation"))

# Species richness
richness <- ggplot(rarefac, aes(x=x, y=y, colour=site)) +
  geom_line(aes(linetype=method), lwd=0.5, 
            data=rarefac.line) +
  scale_colour_discrete(guide=F) +
  scale_x_continuous(limits = c(0,1e+5)) +
  labs(x="Sample size", y="Species richness") +
  theme_classic(base_size = 12) + 
  theme(legend.position="bottom")

# Sample coverage
cover <- fortify(iNEXT, type=2)
cover.point <- cover[which(
  cover$method == "observed"),]
cover.line <- cover[which(
  cover$method != "observed"),]
cover.line$method <- factor(
  cover.line$method,
  c("interpolated", "extrapolated"),
  c("interpolation", "extrapolation"))

coverage <- ggplot(cover, aes(x=x, y=y, colour=site))+ 
  geom_line(aes(linetype=method), lwd=0.5, data=cover.line) +
  scale_colour_discrete(guide=F) +
  scale_x_continuous(limits = c(0,1e+5)) +
  scale_y_continuous(
    breaks = seq(0.9,1,0.05), 
    limits = c(0.9,1)) +
  labs(x="Sample size", y="Sample coverage") +
  theme_classic(base_size = 12) + 
  theme(legend.position="bottom")

############################################################################################

## PLOTS AND TABLES ##

rich <- iNEXT$AsyEst[
  iNEXT$AsyEst$Diversity == "Species richness",]
simp <- iNEXT$AsyEst[
  iNEXT$AsyEst$Diversity == "Simpson diversity",]

Diversity <- data.frame(
  SampleID = iNEXT$DataInfo$site, 
  timepoint = ENV$timepoint, 
  time_h = ENV$time_h,
  condition = ENV$condition,
  sample = ENV$sample,
  Richness = rich$Observed,
  Simpson = simp$Observed)
  
# Export
write.table(
  Diversity, file="Diversity.txt", 
 sep="\t", row.names=F)

# plot Rarefaction and coverage 
plot_grid(
  richness, 
  coverage,  
  label_x = 0.1,
  label_size = 9,
  ncol = 2)

############################################################################################

# Delete temporary datasets
rm(rich, simp, cover)
rm(list=ls(pattern=
  "cover.point.*|cover.line.*|rarefac.point.*|rarefac.line.*"))

############################################################################################

save.image("BAH_particles.Rdata")
