#!/usr/bin/env Rscript

### TreeCertainty: Exploring the phylogenetic signal amongst orthologs of the *Podospora* complex
#############################################################################
# Part of the Snakemake pipeline TreeCertainty.smk
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019-12-30
# Version 1
# =======================================
library(dplyr, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE)
library(cowplot) 
library(reshape2) # for melt

# ============================
# Reading the data
# ============================

# orthos <- read.delim("/Users/Lorena/Dropbox/PhD_UU/Analyses/JohannessonsServer/4_PhylogenyPaper/4b_TreeCertainty_noOutgr/results/TCgenes.txt", header = TRUE, sep="\t")
orthos <- read.delim(snakemake@input[[1]], header = TRUE, sep="\t")

## Genes
# gff <- read.table("/Users/Lorena/Dropbox/PhD_UU/Analyses/JohannessonsServer/4_PhylogenyPaper/4b_TreeCertainty_noOutgr/results/GenesInfo.txt", header = FALSE)
gff <- read.table(snakemake@input[[2]], header = FALSE) 
names(gff) <- c("chromosome", "start", "end", "sense", "panames")

# Names of the orthologs in terms of the P. anserina genes 
# panames <- read.delim("/Users/Lorena/Dropbox/PhD_UU/Analyses/JohannessonsServer/4_PhylogenyPaper/1b_OrthoTrees_noOutgr/filtering/Podan2_1n.clean.txt", header = FALSE, sep="\t")
panames <- read.delim(snakemake@input[[3]], header = FALSE, sep="\t")

allorthos <- cbind(orthos, panames)
names(allorthos)[6] <- "panames"

alldata <- merge(allorthos, gff, by = "panames")
alldata <- alldata %>% mutate(len = end - start) %>% mutate(mid = (end - start)/2 + start)

# ============================
# Plotting TC
# ============================

# Change it into long format for ggplot
lorthos <- orthos %>% melt(id = c("Ortholog")) 

# Plot the distributions
# raw <- ggplot(lorthos %>% filter(variable == "TC" | variable == "TCA") , aes(x = value, colour = variable, fill = variable)) +
#   geom_density(alpha = 0.4) +
#   xlab("TC/TCA") +
#   theme_cowplot() +
#   theme(legend.title=element_blank(), # remove title of legend
#         legend.position=c(0.75,0.3)) 

relatives <- ggplot(lorthos %>% filter(variable == "RelTC" | variable == "RelTCA") , aes(x = value, colour = variable, fill = variable)) +
  geom_density(alpha = 0.4) +
  xlab("Relative TC or TCA") +
  theme_cowplot() +
  theme(legend.title=element_blank(), # remove title of legend
        legend.position=c(0.75,0.3)) +
  geom_vline(xintercept = mean(orthos$RelTC), colour = "coral3") +
  geom_vline(xintercept = mean(orthos$RelTCA), colour = "cadetblue4") +
  # geom_vline(xintercept = 0.5, colour = "white") +
  scale_linetype_manual(values=c(5,5,4))

ggsave(plot = relatives, snakemake@output[[1]], width = 4, height = 4)

# Save it into a table
write.table(orthos %>% filter(RelTC >= 0.5) %>% .$Ortholog, snakemake@output[[2]], sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE) 
write.table(orthos %>% filter(RelTC >= 0.70) %>% .$Ortholog, snakemake@output[[3]], sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE) 
write.table(orthos %>% filter(RelTC >= 0.75) %>% .$Ortholog, snakemake@output[[4]], sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE) 

# ============================
# Where are the genes
# ============================
#### --- Points ----
centromeres <- data.frame(chromosome = c("chromosome_1", "chromosome_2", "chromosome_3", "chromosome_4", "chromosome_5", "chromosome_6", "chromosome_7"),
                          mid = c(4479205, 236468.5, 675115.5, 1236388.5, 2062808.5, 238150.0, 3562141.5), Locus = "Centromere", RelTC = 0 )


# Make a dataframe with the lengths of the Podan2 chromosomes
chrlines <- data.frame(mid =c(0) , endchr = c(8813524, 5165621, 4137471, 3808395, 4734292, 4264132, 4087160), chromosome = c("chromosome_1", "chromosome_2", "chromosome_3", "chromosome_4", "chromosome_5", "chromosome_6", "chromosome_7"))

# Plot all chromosomes on top
genesdist <- ggplot(alldata, aes(x = mid, y = RelTC)) +
  theme_cowplot() + facet_grid(chromosome ~ .) +
  geom_segment(data = chrlines, 
                 aes(x = mid, y = 0, xend = endchr, yend = 0), 
                 colour = "gray",
                 size = 3) +
  geom_point(aes(colour = RelTC), alpha = 0.5) + 
  xlab("Chromosomal location (bp)") +
  scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0.5) +
  geom_point(data = centromeres, aes(x = mid, shape = Locus), colour = "black", size = 5) +
  geom_point(data = alldata %>% filter(RelTC >= 0.75), colour = "purple") # Make them more obvious

ggsave(plot = genesdist, snakemake@output[[5]], width = 9, height = 10)

# ============================
# Exploring the genes
# ============================

# ggplot(alldata, aes(x = RelTC, y = len)) + geom_point() +
#   theme_cowplot()
# 
# ggplot(alldata, aes(x = RelTC, y = log(len))) + geom_point() +
#   theme_cowplot()
# 
# library(ggpubr) # compare means of two or multiple groups; ii) and to automatically add p-values and significance levels to a ggplot 
# 
# ggplot(alldata, aes(x = chromosome, y = RelTC)) + geom_boxplot() +
#   theme_cowplot() +
#   stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.") # In this case, each of the grouping variable levels is compared to all (i.e. basemean).
# 
# ggscatter(alldata, x = "RelTC", y = "len", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "pearson",
#           xlab = "RelTC", ylab = "Length (bp)")
# 
# # Are they even normal?
# shapiro.test(alldata$RelTC[1:5000]) # Not normal
# shapiro.test(alldata$len[1:5000]) # Not normal
# 
# shapiro.test(alldata %>% filter(RelTC >= 0) %>% .$RelTC %>% .[1:5000] %>% log) # Not normal
# shapiro.test(alldata %>% filter(RelTC >= 0) %>% .$len %>% .[1:5000] %>% log) # Not normal
# # But the shape is not bad...
# 
# ggqqplot(alldata$RelTC, ylab = "MPG")
# ggqqplot(alldata$RelTC %>% log, ylab = "MPG")

# library(MASS)
# # run a linear model
# m <- lm(alldata$RelTC ~ alldata$len)
# 
# # run the box-cox transformation
# bc <- boxcox(alldata$RelTC ~ alldata$len)
