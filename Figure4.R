# Figure 4A --------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggsignif)
library(ggpubr)
library(ggplot2)

data = read.table("hgt_number_stat_with_groupinfo.txt",header = T,row.names = 1,check.names = F)

# log
data[,1] = log(data[,1] + 1)

# color palette
color = c("#90A6BB","#FFD579","#9D4E3F")

#
my_comparisons = list(c("NC","AD"),c("NC","Acne"),c("AD","Acne"))

# plot
p = ggboxplot(data, x="Group", y="hgtector_number",color="black",fill="Group", palette = "npg",
            #add = "jitter",
            bxp.errorbar=T,
            bxp.errorbar.width=0.05,
            outlier.shape=NA,
            facet.by = "Site", 
            short.panel.labs = FALSE,
            order = c("NC","AD","Acne"),
            width = 0.7) +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test") +
  stat_compare_means(label.y = 5.5) +
  ylab("log(HGT number + 1)") + 
  scale_fill_manual(values=color)

# zoom
##ylim1<-boxplot.stats(data$hgtector_number)$stats[c(1, 5)] 
#p2=p + coord_cartesian(ylim = ylim1*0.999)
#p2+ylab("log(HGT number + 1)") + scale_fill_manual(values=color)



# Figure 4B --------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(reshape2)
library(ggplot2)
library(dplyr)

data = read.csv("all_singleton_melt.csv",header = T)

#
data$Group = factor(data$Group,levels = c("NC_singleton_gene","AD_singleton_gene","Acne_singleton_gene"))

#
empty_bar = 0
# id
data$id = 1:nrow(data)

# inner annotation
base_anno = group_by(data, COG) %>%
  summarise(start = min(id), end = max(id) - empty_bar) %>%
  mutate(mid = (start + end) / 2)
base_anno

# color palette
color = c("#90A6BB","#FFD579","#9D4E3F")

# plot
p = ggplot() +
  geom_col(data=data, aes(x=id,y=value,fill=Group), position = position_dodge2()) +
  # inner
  geom_segment(data = base_anno, aes(x = start, y = -5, xend = end, yend = -5),
               colour = "grey40") +
  geom_text(data = base_anno, aes(x = mid, y = -10, label = COG), 
            #angle = c(-13, -50, -25, 13), 
            colour = "grey40") +
  coord_polar() +
  ylim(-40,50) +
  theme(
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
  ) +
  scale_fill_manual(values=color)



# Figure 4C --------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggsignif)
library(ggpubr)
library(ggplot2)

data = read.table("accessory_stat_boxplot.txt",header = T,sep = "\t")

# color palette
color=c("#ab6f26","#d6b474")

#
my_comparisons = list(c("high","low"))

#plot
p=ggboxplot(data, x="Group", y="Accessory.gene.number",color="black",
            fill="Group", 
            palette = "npg",
            #add = "jitter",
            bxp.errorbar=T,
            bxp.errorbar.width=0.05,
            outlier.shape=NA,
            short.panel.labs = FALSE,
            width = 0.7)

p1=p+stat_compare_means(comparisons = my_comparisons,method="wilcox.test") +
  ylab("The number of HGT beloging to accessory genes") +
  scale_fill_manual(values=color)
p1


# Figure 4D --------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggplot2)
library(tidyverse)

KEGG = read.table("pathway_enriched_result.txt",header = T,sep = "\t")

#
kegg = KEGG[KEGG$pvalue <0.05,]

#
kegg$Group = factor(kegg$Group,levels = c("HC","AD","Acne"))

#plot
p = ggplot(kegg,aes(x=reorder(Description,Count),y=Count)) + 
  geom_bar(stat='identity',aes(fill=pvalue)) +
  ylab("Gene Counts") +
  xlab("KEGG Pathway")+
  coord_flip()+
  theme(axis.text.y=element_text(color='black', size=12),
        axis.text.x = element_text(size = 10),
        legend.text = element_text(size=10))+
  facet_grid(Group~.,space="free_y",scales="free_y") + 
  #facet_grid(Site~.,space="free_y",scales="free_y") + 
  scale_fill_gradient(high = "#C6C6C6",low ="#DD604B" )

p


# Figure 4E,F ------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggplot2)

data = read.table("dNdS_range_stat_for_group.txt",header = T,sep = "\t")

#
data$Group = factor(data$Group, levels = c("NC","AD","Acne"))

# color palette
color1 = c("#90A6BB","#FFD579","#9D4E3F") # for Group
#color2=c("#8ECFC9","#BEB8DC","#E7DAD2","#999999","#F7E1ED") # for site

#plot
p = ggplot(data,aes(x=position,y=percentage,fill=Group)) + 
  geom_bar(stat='identity',position = "dodge") + #
  geom_smooth(se = FALSE,aes(linetype=as.factor(Group), colour = Group), size=1) +
  scale_linetype_manual(values=c("AF"=1,"AM"=4,"face"=3,"foot"=6)) +
  scale_color_manual(values=color2) +
  ylab("Proportion of genes") +
  xlab("dN/dS") +
  theme_bw() +
  theme(#panel.border = element_blank(), 
    #panel.grid.major = element_line(linetype = 'dashed'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.text = element_text(size=9,color='black'), 
    plot.title = element_text(size=15,color="black"), 
    legend.key = element_blank(), 
    axis.text = element_text(size=10,color='black'), 
    strip.text = element_text(size=12,color="black"), 
    strip.background = element_blank()) +
  scale_fill_manual(values=color1)


p1=p + scale_x_continuous(expand = c(0,0.01),breaks=seq(0, 1.3, 0.1)) + # Y axis gap
  scale_y_continuous(expand = c(0,0),limits = c(0,0.25)) # X axis gap



