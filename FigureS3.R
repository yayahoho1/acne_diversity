# Figure S3 --------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggsignif)
library(ggpubr)
library(ggplot2)
library(egg)

data = read.table("pairwise_phylogenetic_distance_with_groupinfo.txt",header = T,sep = "\t")

# color palette
color = c("#90A6BB","#FFD579","#9D4E3F")
#color2 = c("#75944F","#C75A68")

#
my_comparisons = list(c("NC","AD"),c("NC","Acne"),c("AD","Acne"))
#my_comparisons = list(c("In","Bewteen"))

# plot
p = ggboxplot(data, x="Group", y="distance",color="black",
              fill="Group", 
              palette = "npg",
              bxp.errorbar=T,
              bxp.errorbar.width=0.05,
              outlier.shape=NA,
              short.panel.labs = FALSE,
              order = c("NC","AD","Acne"),
              width = 0.7) +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test") +
  #stat_compare_means(label.y = .7) +
  ylab("Pairwise phylogenetic distance") +
  scale_fill_manual(values=color)

# amplification
#ylim1 = boxplot.stats(data$distance)$stats[c(1, 5)] 
#p1 = p + coord_cartesian(ylim = ylim1*1.5) #

# plot
#p1 + ylab("Pairwise phylogenetic distance") + scale_fill_manual(values=color2)
