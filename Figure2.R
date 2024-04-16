# Figure 2A,B,C ----------------------------------------------------------------
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


# Figure 2D --------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggplot2)
library(cowplot)

NC = read.table("bewteen_and_in_for_single_subject_with_group_NC.txt",header = T)
AD = read.table("bewteen_and_in_for_single_subject_with_group_AD.txt",header = T)
Acne = read.table("bewteen_and_in_for_single_subject_with_group_Acne.txt",header = T)

#
NC$Site = factor(NC$Site,levels = c("Bewteen","AF","AM","Face","Foot"))
AD$Site = factor(AD$Site,levels = c("Bewteen","AF","AM","Face"))
Acne$Site = factor(Acne$Site,levels = c("Bewteen","AF","AM","Face"))


# color palette
color=c("#818e48","#626C83","#D5A551","#C3735C","#999999")

# plot
p1 = ggplot(NC,aes(x=Site,y=Distance)) +
  geom_violin(aes(fill=Site),show.legend = F,colour = "black") +
  geom_boxplot(width=0.02,size=0.1,alpha=0.5,outlier.shape = NA,position=position_dodge(0.9),show.legend = F)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Pairwise phylogenetic distance") +
  theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5))+
  facet_grid(Indivadual~.,scales = "free_x") +
  scale_fill_manual(values = color)

p2 = ggplot(AD,aes(x=Site,y=Distance)) +
  geom_violin(aes(fill=Site),show.legend = F,colour = "black") +
  geom_boxplot(width=0.02,size=0.1,alpha=0.5,outlier.shape = NA,position=position_dodge(0.9),show.legend = F)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Pairwise phylogenetic distance") +
  theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5))+
  facet_grid(Indivadual~.,scales = "free_x") +
  scale_fill_manual(values = color)

p3 = ggplot(Acne,aes(x=Site,y=Distance)) +
  geom_violin(aes(fill=Site),show.legend = F,colour = "black") +
  geom_boxplot(width=0.02,size=0.1,alpha=0.5,outlier.shape = NA,position=position_dodge(0.9),show.legend = F)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Pairwise phylogenetic distance") +
  theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5))+
  facet_grid(Indivadual~.,scales = "free_x") +
  scale_fill_manual(values = color)

# merge plot
plot_grid(p1, p2,p3,align = "v",nrow = 1)



# Figure 2E --------------------------------------------------------------------


# Figure 2F,G ------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(reshape2)
library(ggplot2)

#data=read.table("core_gene_compare_stat_with_Groupinfo.txt",header = T,sep="\t",check.names = F)
data = read.table("Pan_gene_compare_stat_with_Groupinfo.txt",header = T,sep="\t",check.names = F)

#
df_group = melt(data)

#
df_group$variable = factor(df_group$variable,levels = c("Unique to the subject", "Shared with at least one other subject", "Shared with all other subjects"))
df_group$symbol = factor(df_group$symbol,levels = c("NC","AD","Acne"))

# color palette
color = c("#ffd579","#90a6bb","#305757")

#plot
p1 = ggplot(df_group,aes(x=ID,y=value,fill=variable)) + geom_bar(stat='identity',width = 0.5) +
  ylab("Number of genes")+
  xlab("")+
  facet_wrap(~symbol,scales = "free_x") +
  labs(fill = "") +
  #mytheme+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size=14,color='black'), 
        plot.title = element_text(size=15,color="black"), 
        legend.key = element_blank(), 
        axis.text = element_text(size=10,color='black'), 
        strip.text = element_text(size=12,color="black"), 
        strip.background = element_blank(),
        axis.text.x = element_text(size = 14,angle = 90),
        axis.text.y = element_text(size = 12),
        legend.position="bottom") +
  scale_y_continuous(expand = c(0,0.1),limits = c(0,4500)) +
  scale_fill_manual(values = color)


# Figure 2H,I,J ----------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggplot2)
library(tidyverse)

KEGG = read.csv("pathway_enrichment_result.csv",header = T)

#
kegg = KEGG[KEGG$pvalue <0.05,]

#
kegg$Group = factor(kegg$Group,levels = c("AD","Acne"))
#kegg$Group = factor(kegg$Group,levels = c("NC","Acne"))
#kegg$Group = factor(kegg$Group,levels = c("NC","AD"))

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





