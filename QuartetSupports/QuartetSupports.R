#!/usr/bin/env Rscript

### Quartet distribution in the phylogeny of the *Podospora* complex
#############################################################################
# The raw data I collected from the output of the Snakemake pipelines TreeCertainty.smk and Phylogeny.smk
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020-01-13
# Version 1
# =======================================
library(dplyr, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE)
library(cowplot) 
library(reshape2) # for melt

rawqq <- read.csv("/Users/Lorena/Dropbox/PhD_UU/Writing/PaperPhylogeny/QuartetSupports.csv", header = TRUE)

longqq <- rawqq %>% melt(id = c("Taxa", "Strategy", "Opinionated", "n", "MinUFBoot"))
names(longqq) <- c("Taxa", "Strategy", "Opinionated", "n", "MinUFBoot", "clade", "quartet_support")

ggplot(longqq, aes(x = MinUFBoot, y = quartet_support, colour = as.factor(Opinionated))) + 
  scale_y_continuous(limits = c(0, 1)) + 
  scale_x_continuous(limits = c(0, 100)) + 
  scale_size_continuous(breaks = c(20, 100, 1000, 2000, 4000, 6000, 7000)) +
  geom_point(aes(size = n)) + 
  facet_grid(Taxa ~ clade) +
  geom_line(aes(linetype = as.factor(Opinionated))) +
  ylab("Quartet support") +
  xlab("Minimum UFBoot support") +
  guides(colour=guide_legend(title="Minimum\nrelative TC"), linetype=guide_legend(title="Minimum\nrelative TC")) +
  theme_bw()

ggsave("/Users/Lorena/Dropbox/PhD_UU/Writing/PaperPhylogeny/FigS2_QuartetSupports.pdf", width = 12, height = 5)


# plotqq <- function(data = rawqq, clade){
#   p <- ggplot(data, aes(x = MinUFBoot, y = get(clade), colour = as.factor(Opinionated))) + # The get() is used so I can retrive it as a variable
#     scale_y_continuous(limits = c(0, 1)) + 
#     geom_point(aes(size = n)) + 
#     facet_grid(Taxa ~ .) +
#     geom_line() +
#     ylab("Quartet support") +
#     guides(colour=guide_legend(title="Fraction of genes")) +
#     ggtitle(clade)
#   return(p)
# }
# 
# plotqq(clade = "ans.pauci")
# plotqq(clade = "pseudos")
