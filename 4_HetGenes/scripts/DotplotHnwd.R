#!/usr/bin/env Rscript

# Dotplot of hnwd genes
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020-02-07
# =======================================
Version <- 1.0
# ============================
# Load the necessary libraries
# ============================
library(ggplot2, quietly = TRUE)
library(dplyr, warn.conflicts = FALSE)
library(cowplot)

# ============================
# Data
# ============================
# The query is Podan2 in this case
alignment <- snakemake@input[[1]]
gff_ref <- snakemake@input[[2]]
gff_query <- snakemake@input[[3]]

selfalign <- snakemake@input[[4]]
gff_self <- snakemake@input[[5]]

# alignment <- "/Users/Lorena/Dropbox/PhD_UU/Analyses/JohannessonsServer/4_PhylogenyPaper/5_HetGenes/mummer/CBS415.72m_hnwd-pc13.fa-vs-Podan2_hnwd-pc13.fa.coords"
# gff_ref <- "/Users/Lorena/Dropbox/PhD_UU/Analyses/JohannessonsServer/4_PhylogenyPaper/5_HetGenes/data/dotplots/Podan2_hnwd-pc13.gff3"
# gff_query <- "/Users/Lorena/Dropbox/PhD_UU/Analyses/JohannessonsServer/4_PhylogenyPaper/5_HetGenes/data/dotplots/CBS415.72m_hnwd-pc13.gff3"
# 
# selfalign <- "/Users/Lorena/Dropbox/PhD_UU/Analyses/JohannessonsServer/4_PhylogenyPaper/5_HetGenes/mummer/CBS415.72m_hnwd-pc13_insertion.fa-vs-CBS415.72m_hnwd-pc13_insertion.fa.coords"
# gff_self <- "/Users/Lorena/Dropbox/PhD_UU/Analyses/JohannessonsServer/4_PhylogenyPaper/5_HetGenes/data/dotplots/CBS415.72m_hnwd-pc13_insertion.gff3"

# ============================
# Functions
# ============================
gff2genenames <- function(gff){
  genenames <- c()
  for(i in 1:nrow(gff)){
    attrib <- gff[i,9] %>% as.character() %>% strsplit(.,";") %>% .[[1]]
    name <- attrib[pmatch("Name=", attrib)] %>% strsplit(.,"=") %>% .[[1]] %>% .[2]
    genenames <- c(genenames, name)
  }
  gff <- cbind(gff, genenames) %>% select(V1, V4, V5, genenames)
  return(gff)
}

plotmummers <- function(dotplot, coords, linesize = 1, pointsize= 1){
  # Plot the alignments 
  for (i in 1:dim(coords)[1]){
    currentrow <- coords[i,]
    if (currentrow$S2 > currentrow$E2){ # Reverse
      colorsito <- "coral2"
    }
    else{ # Forward
      colorsito <- "black"
    }
    dotplot <- dotplot + geom_segment(data = currentrow, aes(x = S1, xend = E1, y = S2, yend = E2), colour = colorsito, size = linesize) +
      geom_point(data = currentrow, aes(x=S1, y = S2), colour = colorsito, size = pointsize) + geom_point(data = currentrow, aes(x=E1, y = E2), colour = colorsito, size = pointsize)
  }
  return(dotplot)
}

# ============================
# Processing
# ============================
coords <- read.delim(alignment, header = FALSE)
# S1 and E1 are reference, S2 and E2 are query
names(coords) <- c("S1", "E1", "S2", "E2", "LEN_R", "LEN_Q", "IDY", "COV_R", "COV_Q", "REF", "QUERY")

gffr <- read.delim(gff_ref, header = FALSE) %>% filter(V3 == "gene") %>% gff2genenames() %>% cbind(., S1 = 0)
names(gffr) <- c("contig", "LEN_Q", "end", "Locus", "S1")
gffq <- read.delim(gff_query, header = FALSE) %>% filter(V3 == "gene") %>% gff2genenames() %>% cbind(., LEN_Q = 0)
names(gffq) <- c("contig", "S1", "end", "Locus", "LEN_Q")

## self alignment
selfcoords <- read.delim(selfalign, header = FALSE)
names(selfcoords) <- c("S1", "E1", "S2", "E2", "LEN_R", "LEN_Q", "IDY", "COV_R", "COV_Q", "REF", "QUERY")
gffs <- read.delim(gff_self, header = FALSE) %>% filter(V3 == "gene") %>% gff2genenames()
names(gffs) <- c("contig", "KEY", "end", "Locus")

# ============================
# Plotting
# ============================
# Make base plot
dotplot <- ggplot(coords, aes(x=S1, y=LEN_Q)) + geom_blank() +
  theme_bw() +
  theme(strip.text.y = element_blank(), panel.border=element_blank()) + # Remove the gray label on the side and border of plot
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  ylab(expression(paste(italic("P. anserina")~'(Podan2)'))) + 
  xlab(expression(paste(italic("P. pseudocomata")~'(CBS415.72m)')))

# Add annotation
dotplot <- dotplot +
  geom_segment(data = gffr, aes(x = S1, xend = 1, yend = end, colour = Locus), size = 5) +
  geom_segment(data = gffq, aes(x = S1, xend = end, yend = 1, colour = Locus), size = 5) +
  scale_color_manual(values=c("#8dd35fff", "#aa87deff", "#d35f8d2f", "azure3", "azure4")) 

dotplot <- plotmummers(dotplot, coords, linesize = 1, pointsize = 0.5)

# Save plot
ggsave(dotplot, file = snakemake@output[[1]], width = 6, height = 3.5)
# ggsave(dotplot, file = "/Users/Lorena/Dropbox/PhD_UU/Analyses/JohannessonsServer/4_PhylogenyPaper/5_HetGenes/results/dotplot.pdf", width = 6, height = 3.5)

## Self-alignment
# Make base plot
selfdotplot <- ggplot(selfcoords, aes(x=S1, y=LEN_Q)) + geom_blank() +
  theme_bw() +
  theme(strip.text.y = element_blank(), panel.border=element_blank()) + # Remove the gray label on the side and border of plot
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  xlab(expression(paste(italic("hnwd-pc13")~' insertion'))) + 
  ylab(expression(paste(italic("hnwd-pc13")~' insertion')))

# Add annotation
selfdotplot <- selfdotplot +
  geom_segment(data = gffs %>% mutate(LEN_Q = KEY, S1 = 0), aes(x = S1, xend = 1, yend = end, colour = Locus), size = 5) +
  geom_segment(data = gffs %>% mutate(LEN_Q = 0, S1 = KEY), aes(x = S1, xend = end, yend = 1, colour = Locus), size = 5) +
  scale_color_manual(values=c("#8dd35fff", "#aa87deff", "#d35f8d2f"))

selfdotplot <- plotmummers(selfdotplot, selfcoords, linesize = 1, pointsize = 0.5)

# Save plot
ggsave(selfdotplot, file = snakemake@output[[2]], width = 5, height = 3.5)



