# Figure S12B ------------------------------------------------------------------
# Load required packages and data ----------------------------------------------
library(ggplot2)
library(reshape2)

data = read.csv("target_metabolome.csv",header = T,row.names = 1,check.names = F)
groupinfo = read.csv("group_info.csv",header = T,row.names = 1)

# normalization
df = scale(data,center = T,scale =T )
df = as.data.frame(df)
df$Group = groupinfo$Group
df_long = melt(df)

#
df_long$Group = factor(df_long$Group,levels=c("NC","AD","Acne"))

# color palette
color = c("#90A6BB","#FFD579","#9D4E3F")

p = ggplot(df_long, aes(variable, value, shape=Group,color=Group)) + geom_point(size=5,alpha=.6) +
  coord_flip()+
  theme_bw()+
  scale_color_manual(values = color)

# Figure S12C,D ----------------------------------------------------------------
# Load required packages and data ----------------------------------------------
library(ggsignif)
library(ggpubr)
library(ggplot2)

data = read.csv("Propionic_acid.csv",header = T,row.names = 1,check.names = F)
#data = read.csv("Butyric_acid.csv",header = T,row.names = 1,check.names = F)

# color palette
color=c("#90A6BB","#FFD579","#9D4E3F")

#
my_comparisons = list(c("NC","AD"),c("NC","Acne"),c("AD","Acne"))

# plot
p=ggboxplot(data, x="Group", y="Propionic_acid",color="black",fill="Group", palette = "npg",
#p=ggboxplot(data, x="Group", y="Butyric_acid",color="black",fill="Group", palette = "npg",            
            #add = "jitter",
            bxp.errorbar=T,
            bxp.errorbar.width=0.05,
            short.panel.labs = FALSE,
            order = c("NC","AD","Acne"),
            width = 0.7) + 
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test") +
  stat_compare_means(label.y = 3) +
  ylab("Normalized intensity") + 
  scale_fill_manual(values=color)
  
p


# Figure S12E -------------------------------------------------------------------
# Load required packages and data ----------------------------------------------
library(ggplot2)
library(ropls)
library(ggforce)
library(factoextra)

data=read.csv("target_metabolome.csv",header = T,row.names = 1,check.names = F)
groupinfo=read.csv("group_info.csv",header = T,row.names = 1)

# normalization
data1=scale(data)

# 进行PCA分析
res.pca <- prcomp(data1, scale = TRUE)
#res.pca
summary(res.pca)

#
screeplot(res.pca, npcs = 10, type = "lines")

# opls
plsda <- opls(data1, groupinfo$Group, predI = NA, orthoI = 0)
#plsda <- opls(data, groupinfo$Group, predI = 7, orthoI = 0)
#
scoreMN <- plsda@scoreMN
scoreMN = data.frame(scoreMN)
scoreMN <- cbind(scoreMN, groupinfo$Group)
colnames(scoreMN)[3]="Group"
scoreMN$Group=factor(scoreMN$Group,levels = c("NC","AD","Acne"))

#
x_lab <- plsda@modelDF[1, "R2X"] * 100
y_lab <- plsda@modelDF[2, "R2X"] * 100

# color palette
color=c("#90A6BB","#FFD579","#9D4E3F")

# plot
p=ggplot(scoreMN, aes(x=p1,y=p2,color=Group,shape=Group)) +
  geom_point(alpha=.9, size=4) +
  #stat_ellipse(level=0.95) +
  #theme_classic()
  theme_bw() +
  labs(x = paste0("PC1(", x_lab, ")"), y = paste0("PC2(", y_lab, ")")) +
  geom_ellipse(aes(x0 = 0, y0 = 0, a = 4, b = 3, angle = 0),
               linetype="dashed",
               color="grey") +
  coord_fixed() +
  scale_color_manual(values = color)
p


# Figure S12F -------------------------------------------------------------------

setwd("E:/project/06.STU004_piyansuo_补充_20220711/NC_AD_Acne/202210/all/文章图片/后期/C. acnes文章修改_V10_20231126_figure修改/Figure 7/Figure7C")

# Load required packages and data ----------------------------------------------
library(ggplot2)
library(tidyverse)

KEGG=read.csv("pathway_enrichment_top20.csv",header = T)

# plot
p3 <- ggplot(KEGG,aes(x=reorder(pathway_name,-Pvalue),y=Impact)) + geom_point(aes(size=Hits,color=Pvalue)) +
  ylab("Impact") +
  coord_flip()+
  theme(
    axis.text.y=element_text(color='black', size=14),
    axis.text.x = element_text(size = 12),
    legend.text = element_text(size=12)) + 
  scale_color_gradient(high = "#C6C6C6",low ="#DD604B" )

p3


# Figure S12G,H ----------------------------------------------------------------
# Load required packages and data ----------------------------------------------
library(ggplot2)
library(ggrepel)

#data=read.csv("differential_metabolites_AD_vs_HC.csv",header = T,check.names = F)
data=read.csv("differential_metabolites_Acne_vs_HC.csv",header = T,check.names = F)

# color palette
color = c(Down="#CBD289",Normal="#D2D2D2",Up="#C47A9B")

#
data$Regulation <- ifelse(data$Pvalue<0.05 & abs(data$log2FC)>= 0,ifelse(data$log2FC >= 0,'Up','Down'),'Normal')

# plot
p=ggplot(data,aes(log2FC, -log10(Pvalue),col=Regulation)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_point(aes(size=VIP),alpha=.6) +
  scale_colour_manual(values = color) +
  theme_bw() + 
  xlab("Log2FC") +
  ylab("-Log10(Pvalue)")

p

# add label
p1=p+ggrepel::geom_text_repel(
  aes(label=label),data,
  show.legend=F,
  max.overlaps = 50,
  box.padding=unit(0.5, "lines"),
  segment.color = "black", segment.size = 0.5,
  arrow = arrow(length=unit(0.01, "npc"))
)

p1



