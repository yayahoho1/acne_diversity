# Figure 5A --------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(reshape2)
library(ggplot2)

data = read.csv("COG_abundance_summary.csv",header = T)

# color palette
color11 = c("#D36511", "#5A7BBC", "#E96C7B", "#50D4C8", "#8DC862", "#FDD347", "#F3A76D","#6D90D5", "#B64C58", "#0C6961", "#467128", "#917001", "#977810", "#254380","#E54C5E", "#30C0B4", "#75BD42", "#F2BA02", "#EE822F", "#4874CB")

# 
data_long = melt(data)

#plot
p = ggplot() + geom_bar(data =data_long, aes(x = variable, y = value, fill = factor(COG)),
                      stat = "identity",
                      position = "fill",width = 0.7)+ #
  scale_fill_manual(values=color11)+
  theme_bw()+
  theme(#panel.border = element_blank(), 
    #panel.grid.major = element_line(linetype = 'dashed'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.text = element_text(size=9,color='black'), 
    plot.title = element_text(size=15,color="black"), 
    legend.key = element_blank(), 
    axis.text = element_text(size=10,color='black'), 
    strip.text = element_text(size=12,color="black"), 
    strip.background = element_blank())

p
p1 = p + scale_y_continuous(expand = c(0,0))
p1



# Figure 5B --------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggplot2)
library(tidyverse)

KEGG = read.csv("pathway_enrichment.csv",header = T)

#
kegg = KEGG[KEGG$pvalue <0.05,]

#
kegg$Group = factor(kegg$Group,levels = c("AD_vs_NC","Acne_vs_NC"))

# plot
p = ggplot(kegg,aes(x=Site,y=Description)) + 
  geom_point(aes(size=GeneRatio,color=pvalue)) +
  theme(
    axis.text.y=element_text(color='black', size=18),
    axis.text.x = element_text(size = 10),
    legend.text = element_text(size=10))+
  facet_grid(Symbol~Group,space="free_y",scales="free_y") +
  scale_color_gradient(high = "#C6C6C6",low ="#DD604B" )

p


# Figure 5C --------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ComplexHeatmap)
library(circlize)

data = read.table("gene_counts_matrix_vsd.txt",header = T,row.names = 1,sep = "\t",check.names = F)

#
sampleDists <- dist(t(data))
sampleDistMatrix <- as.matrix(sampleDists)

# color palette
col_fun = colorRamp2(c(0,80),c("#DD604B","#f2f2f2"))

# plot
ComplexHeatmap::Heatmap(sampleDistMatrix, 
                        col = col_fun,
                        row_names_side = "right",
                        show_row_names = T,
                        column_names_gp = gpar(fontsize = 10),
                        row_names_gp = gpar(fontsize = 10),
                        show_column_names = T,
                        cluster_rows = T,cluster_columns = T,
                        cluster_column_slices = F,
                        row_names_max_width = unit(10,"cm"),
                        rect_gp = gpar(col = "grey60", lwd = 0.000001),
                        heatmap_legend_param = list(title= "Euclidean distance",
                                                    title_position = "topcenter", 
                                                    legend_direction= "horizontal",
                                                    legend_height=unit(5,"cm"))
)




# Figure 5D --------------------------------------------------------------------

# Load required packages and data ----------------------------------------------

library(ggplot2)
library(tidyverse)

KEGG = read.csv("diffgene_enrichment.csv",header = T)

#
kegg = KEGG[KEGG$pvalue <0.05,]

#
kegg$Group = factor(kegg$Group,levels = c("0.01","0.25"))

# plot
p = ggplot(kegg,aes(x=reorder(Description,Count),y=Count)) + 
  geom_bar(stat='identity',aes(fill=pvalue)) +
  ylab("Gene Counts") +
  xlab("KEGG Pathway")+
  coord_flip()+
  theme(axis.text.y=element_text(color='black', size=12),
        axis.text.x = element_text(size = 10),
        legend.text = element_text(size=10))+
  facet_grid(Group~.,space="free_y",scales="free_y") +
  scale_fill_gradient(high = "#C6C6C6",low ="#DD604B" )

p



