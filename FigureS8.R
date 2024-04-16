# Figure S8 --------------------------------------------------------------------

#setwd("E:/project/08.STU004_transcriptom/202212/04.featureCounts/unique_reads/Distance/Bary_distance")
#setwd("E:/project/02.metabolim/非靶向/metabolim_tree/Distance/Bray_ditance")
setwd("E:/project/06.STU004_piyansuo_补充_20220711/NC_AD_Acne/202210/all/review/20240226_回复讨论/reviewer1/point3/pangenome_81strains_for_multi-omics/site_vs_group")

# Load required packages and data ----------------------------------------------

library(ggsignif)
library(ggpubr)
library(ggplot2)
library(egg)

data=read.table("pairwise_phylogenetic_distance.txt",header = T, sep = "\t")
#data=read.table("Bray-Curtis_distance.txt",header = T, sep = "\t")

# color palette
color1=c("#90A6BB","#FFD579","#9D4E3F") # for Group
#color2=c("#ab6f26","#d6b474")

#
my_comparisons = list(c("NC","AD"),c("NC","Acne"),c("AD","Acne"))
#my_comparisons = list(c("Bewteen_sites","Bewteen_groups"))
#my_comparisons = list(c("In_individual","Between_individual"))

#
p=ggboxplot(data, x="Group", y="distance",
            color="black",
            fill="Group", 
            palette = "npg",
            #add = "jitter",
            bxp.errorbar=T,
            bxp.errorbar.width=0.05,
            outlier.shape=NA, 
            short.panel.labs = FALSE,
            width = 0.7) +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test")

p1=p+stat_compare_means(label.y = .6) +
  #ylab("Pairwise Bray-Curtis distance") + 
  ylab("Pairwise phylogenetic distance") + 
  scale_fill_manual(values=color6) + 
  scale_y_continuous(limits = c(0,0.5))

p1



# Local amplification
#ylim1<-boxplot.stats(data$distance)$stats[c(1, 5)] 

#p1=p + coord_cartesian(ylim = ylim1*1.5) #

#p1 + ylab("Pairwise phylogenetic distance") + scale_fill_manual(values=color2)

















