#Figure S7A,B,H-J --------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggplot2)

data=read.table("dNdS_range_stat_for_group.txt",header = T,sep = "\t")

# color palette
color1=c("#90A6BB","#FFD579","#9D4E3F") # for Group
#color2=c("#baced9","#8b9aa2","#d1b19d","#de8387","#bcb9ad") # for site

#
data$Group=factor(data$Group, levels = c("HC","AD","Acne"))
#data$Group=factor(data$Group, levels = c("AF","AM","Face"))

# plot
p = ggplot(data,aes(x=position,y=percentage,fill=Group)) + 
  geom_bar(stat='identity',position = "dodge") + #
  geom_smooth(se = FALSE,aes(linetype=as.factor(Group), colour = Group), size=1) +
  scale_linetype_manual(values=c("AF"=1,"AM"=4,"face"=3,"foot"=6)) +
  scale_color_manual(values=color2) +
  ylab("Proportion of genes") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.text = element_text(size=9,color='black'), 
    plot.title = element_text(size=15,color="black"), 
    legend.key = element_blank(), 
    axis.text = element_text(size=10,color='black'), 
    strip.text = element_text(size=12,color="black"), 
    strip.background = element_blank()) +
  scale_fill_manual(values=color2) + 
  scale_x_continuous(expand = c(0,0.01),breaks=seq(0, 1.3, 0.1)) + # Y axis gap
  scale_y_continuous(expand = c(0,0),limits = c(0,0.25)) + # X axis gap
  xlab("dN/dS")
p




# Figure S7C-E -----------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggsignif)
library(ggpubr)
library(ggplot2)
library(egg)

data = read.table("accessory_dNdS_summary.txt",header = T,sep = "\t")

# color palette
color1=c("#9D4E3F","#FFD579","#90A6BB")

#
my_comparisons = list(c("HC","AD"),c("HC","Acne"),c("AD","Acne"))

# log
data[,2] = log(data[,2] + 0.01)

# plot
p = ggboxplot(data, x="Group", y="dNdS",color="black",
            fill="Group", 
            palette = "npg",
            #add = "jitter",
            bxp.errorbar=T,
            bxp.errorbar.width=0.05,
            outlier.shape=NA,
            short.panel.labs = FALSE,
            order = c("HC","AD","Acne"),
            width = 0.7) +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test") +
  stat_compare_means(label.y = 5.5) +
  scale_fill_manual(values = color)
p


# Figure S7F,G -----------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggplot2)
library(tidyverse)

KEGG = read.csv("pathway_enrichment_result.csv",header = T)

#
kegg = KEGG[KEGG$pvalue <0.05,]

#
kegg$Group = factor(kegg$Group,levels = c("HC","AD","Acne"))

# plot
p3 = ggplot(kegg,aes(x=reorder(Description,Count),y=Count)) + 
  geom_bar(stat='identity',aes(fill=pvalue)) +
  ylab("Gene Counts") +
  xlab("KEGG Pathway")+
  coord_flip()+
  theme(axis.text.y=element_text(color='black', size=12),
        axis.text.x = element_text(size = 10),
        legend.text = element_text(size=10))+
  facet_grid(Group~.,space="free_y",scales="free_y") +
  scale_fill_gradient(high = "#C6C6C6",low ="#DD604B" )





