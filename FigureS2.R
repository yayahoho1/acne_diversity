# Figure S2
# Load required packages and data ----------------------------------------------

library(vegan)
library(ggplot2)
library(patchwork)

# gene content
data = read.table("gene_presence_absence.Rtab",row.names = 1,header = T, check.names = F)

# sample info
groupinfo = read.table("sample_info.txt",row.names = 1,header = T)

# Determine whether the sample order is consistent
all(rownames(groupinfo) %in% colnames(data))
all(rownames(groupinfo) == colnames(data))

# Reorder the sample
data = data[, rownames(groupinfo)]
all(rownames(groupinfo) == colnames(data))

data = as.matrix(data)

# PCA
pca1 = prcomp(t(data))
df1 = pca1$x # PC score
df1 = as.data.frame(df1) 

#
summ1 = summary(pca1)
xlab1 = paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 = paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")

# group info
Site = groupinfo$Site
Group = groupinfo$Individual

# color pelette
colorE = c("#C79757FF", "#8D41C4FF", "#6DC954FF", "#CA5896FF", "#CACF4FFF", "#7374C8FF", "#D24F3BFF", "#87CFACFF", 
         "#442F56FF", "#4A603BFF", "#A0A9C0FF", "#7C403AFF" )  

colorF = c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF","#8C564BFF","#E377C2FF","#7F7F7FFF",
         "#BCBD22FF","#17BECFFF","#AEC7E8FF","#FFBB78FF", "#98DF8AFF","#FF9896FF","#C5B0D5FF","#C49C94FF","#C7C7C7FF",
         "#F7B6D2FF","#DBDB8DFF","#9EDAE5FF")

colorG = c("#A6CEE3FF","#1F78B4FF","#B2DF8AFF", "#33A02CFF", "#FB9A99FF", "#E31A1CFF", "#FDBF6FFF","#FF7F00FF","#CAB2D6FF",
         "#6A3D9AFF","#FFFF99FF","#B15928FF")



# PERMANOVA
# set seed
set.seed(321)

dune.div1 = adonis2(t(data) ~ Individual, data = groupinfo, permutations = 999, method="bray")
dune.div2 = adonis2(t(data) ~ Site, data = groupinfo, permutations = 999, method="bray")

dune_adonis3 = paste0("Individual adonis R2: ",round(dune.div1$R2,2), ", P-value: ", dune.div1$`Pr(>F)`,"; ",
                       "Site adonis R2: ",round(dune.div2$R2,2), ", P-value: ", dune.div2$`Pr(>F)`)

# plot
p = ggplot(data = df1,aes(x = PC1,y = PC2,color = Group,shape=Site))+
  geom_point(size = 3,alpha=.5)+
  labs(x = xlab1,y = ylab1,color = "Individual",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw() +
  labs(title=dune_adonis3) +
  scale_color_manual(values = colorE)

# Local amplification
ppp = p + xlim(17.5,19) + ylim(-7.5,-6) + theme(legend.position = 'none') +
  xlab("") + ylab("") + labs(title = "")
ppp

#p + inset_element(ppp, 0.01, 0.6, 0.6, 0.95, on_top = TRUE)

