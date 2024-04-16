# Figure S10A ------------------------------------------------------------------
# Load required packages and data ----------------------------------------------
library(ggplot2)
library(vegan)
library(ape)

data=read.table("transcriptome_gene_counts_vsd.txt",row.names = 1,header = T,check.names = F)
groupinfo = read.table("group_info.txt",row.names = 1,header = T)

# Determine whether the sample order is consistent
all(rownames(groupinfo) %in% colnames(data))
all(rownames(groupinfo) == colnames(data))

# Reorder the sample
data = data[, rownames(groupinfo)]
all(rownames(groupinfo) == colnames(data))

# color palette
color = c("#9D4E3F","#FFD579","#90A6BB")

# bray-curtis distance
dune_dist = vegdist(t(data), method="bray", binary=F)

# PcoA
res = pcoa(dune_dist,correction="none")

res_data = res$vectors[,1:2]
res_data = data.frame(res_data)
colnames(res_data) = c('PCoA1','PCoA2')
eig = as.numeric(res$value[,1])

dune_pcoa_result = cbind(res_data,groupinfo)
#head(dune_pcoa_result)
dune_pcoa_result$Group = factor(dune_pcoa_result$Group)

# PERMANOVA
# set seed
set.seed(321)

# Row represent sample
dune.div = adonis2(t(data) ~ Group, data = groupinfo, permutations = 999, method="bray")
dune_adonis = paste0("adonis R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)

# plot 
ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, color=Group,shape=Group)) +  
  geom_point(alpha=.9, size=5) +
  
  labs(x=paste("PCoA 1 (",format(100* eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100*eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=dune_adonis)+ 
  stat_ellipse(level=0.95) +
  theme_classic() +
  scale_color_manual(values = color)



# Figure S10B,D ----------------------------------------------------------------
# Load required packages and data ----------------------------------------------
library(ggVolcano)
library(ggrepel)

deg_data = read.csv("AD_vs_NC_All_WilcoxonTest.csv",row.names = 1,header = T)

#
#data = add_regulate(deg_data, log2FC_name = "log2FoldChange",
#                     fdr_name = "padj",log2FC = 1, fdr = 0.05)
#
#data$row = rownames(data)
#head(data)
#p = ggvolcano(data, x = "log2FoldChange", y = "padj",
#                fills = c("#CBD289","#D2D2D2","#C47A9B"),
#                colors = c("#CBD289","#D2D2D2","#C47A9B"),
#                pointSize = 4,
#                label = "label", label_number = 500, legend_position = "UL", output = FALSE)

# # color palette
color = c(Down="#CBD289",Normal="#D2D2D2",Up="#C47A9B")

#
deg_data$Regulation = ifelse(deg_data$padj<0.05 & abs(deg_data$log2FoldChange)>= 1,ifelse(deg_data$log2FoldChange >= 1,'Up','Down'),'Normal')

# plot
p1 = ggplot(deg_data,aes(log2FoldChange, -log10(padj),col=Regulation)) +
  geom_hline(yintercept = -log10(0.05), lty=5,col="black",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=5,col="black",lwd=0.6) +
  geom_point(size=4,alpha=1) +
  scale_colour_manual(values = color) +
  theme_bw() + 
  xlab("Log2FC") +
  ylab("-Log10FDR")

p1

# add label
p2 = p1+ggrepel::geom_text_repel(
  aes(label=label),deg_data,
  show.legend=F,
  max.overlaps = 500,
  box.padding=unit(0.5, "lines"),
  segment.color = "black", segment.size = 0.5,
  arrow = arrow(length=unit(0.01, "npc"))
)

p2


# Figure S10C,E ----------------------------------------------------------------
# Load required packages and data ----------------------------------------------
library(ggplot2)

COG = read.table("result.COG.matrix",sep = "\t")

# plot
p = ggplot(COG,aes(x=V1,y=V2)) + geom_bar(stat='identity',fill="#DD604B") +
  ylab("Gene Counts") +
  #xlab("KEGG Pathway")+
  coord_flip()+
  theme(axis.text.y=element_text(color='black', size=10)) + 
  scale_y_continuous(expand = c(0,0),breaks  = seq(0,14,2))
p


# Figure S10F ------------------------------------------------------------------
# Load required packages and data ----------------------------------------------
library(ggplot2)
library(tidyverse)
library("vegan")

KEGG = read.csv("diffgene_pathway_enrichment_AD_vs_Acne.csv",header = T)

#
kegg = KEGG[KEGG$pvalue <0.05,]

# plot
p3 = ggplot(kegg,aes(x=Site,y=Description)) + geom_point(aes(size=GeneRatio,color=pvalue)) +
  theme(
    axis.text.y=element_text(color='black', size=18),
    axis.text.x = element_text(size = 10),
    legend.text = element_text(size=10))+
  facet_grid(Symbol~Group,space="free_y",scales="free_y") +
  scale_color_gradient(high = "#C6C6C6",low ="#DD604B" )



# Figure S10G ------------------------------------------------------------------
# Load required packages and data ----------------------------------------------
library(ggplot2)
library(vegan)
library(pairwiseAdonis)
library(ggpubr)
library(patchwork)

data = read.table("transcriptome_gene_counts_vsd.txt",row.names = 1,header = T,check.names = F)
groupinfo = read.table("group_info.txt",row.names = 1,header = T)


# Determine whether the sample order is consistent
all(rownames(groupinfo) %in% colnames(data))
all(rownames(groupinfo) == colnames(data))

# Reorder the sample
data = data[, rownames(groupinfo)]
all(rownames(groupinfo) == colnames(data))

# color palette
color=c("#336274","#17bb7f","#eba957","#a54d41")

# bray-curtis distance
dune_dist = vegdist(t(data), method="bray", binary=F)

# PcoA
res = pcoa(dune_dist,correction="none")

res_data = res$vectors[,1:2]
res_data = data.frame(res_data)
colnames(res_data) = c('PCoA1','PCoA2')
eig = as.numeric(res$value[,1])

dune_pcoa_result = cbind(res_data,groupinfo)
#head(dune_pcoa_result)
dune_pcoa_result$Group = factor(dune_pcoa_result$Group)

# PERMANOVA
# set seed
set.seed(321)

# Row represent sample
dune.div = adonis2(t(data) ~ Group, data = groupinfo, permutations = 999, method="bray")
dune_adonis = paste0("adonis R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)
dune_adonis
# plot 
p = ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, color=Group,shape=Group)) +  
  geom_point(alpha=.8, size=5) +
  
  labs(x=paste("PCoA 1 (",format(100* eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100*eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=dune_adonis)+ 
  stat_ellipse(level=0.95) +
  theme_classic() +
  scale_color_manual(values = color)

p
# This is a wrapper function for multilevel pairwise comparison 
# using adonis() from package 'vegan'. 
# The function returns adjusted p-values using p.adjust().
dune.pairwise.adonis = pairwise.adonis(x=t(data), factors=groupinfo$Group, sim.function = "vegdist",
                                        sim.method = "bray",
                                        p.adjust.m = "BH",
                                        reduce = NULL,
                                        perm = 999)

dune.pairwise.adonis

#
#dune_adonis_pairs = paste0("pairwise.adonis R2: ",round(dune.pairwise.adonis$R2,3), "; P-value: ", dune.pairwise.adonis$p.value)

#
tab2 = ggtexttable(dune.pairwise.adonis[,c("pairs","R2","p.value","p.adjusted")], rows = NULL, 
                    theme = ttheme("blank")) %>% 
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1)  %>% 
  tab_add_hline(at.row = nrow(dune.pairwise.adonis)+1, row.side = "bottom", linewidth = 1)  

p2 = p + tab2 
p2

















