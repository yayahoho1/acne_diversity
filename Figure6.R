# Figure 6A --------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggplot2)
library(ggrepel)

data = read.csv("differential_metabolites.csv",header = T,check.names = F)

# color
data$Regulation = ifelse(data$Pvalue<0.05 & abs(data$log2FC)>= 0,ifelse(data$log2FC >= 0,'Up','Down'),'Normal')
color = c(Down="#CBD289",Normal="#D2D2D2",Up="#C47A9B")

# volcano plot
p = ggplot(data,aes(log2FC, -log10(Pvalue),col=Regulation)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_point(aes(size=VIP),alpha=.6) +
  scale_colour_manual(values = color) +
  theme_bw() +
  xlab("Log2FC") +
  ylab("-Log10(Pvalue)") +
  ggrepel::geom_text_repel(
  aes(label=label),
  show.legend=F,
  max.overlaps = 50,
  box.padding=unit(0.5, "lines"),
  segment.color = "black", segment.size = 0.5,
  arrow = arrow(length=unit(0.01, "npc"))
)
p


# Figure 6B --------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggsignif)
library(ggpubr)
library(ggplot2)

data = read.csv("Carnosine.csv",row.names = 1,header = T,check.names = F)

# color palette
color = c("#90A6BB","#FFD579","#9D4E3F")

#
my_comparisons = list(c("NC","AD"),c("NC","Acne"),c("AD","Acne"))

# log2
data$Carnosine = log2(data$Carnosine)

# plot
p = ggboxplot(data, x="Group", y="Carnosine",
            color="black",
            fill="Group", 
            palette = "npg",
            bxp.errorbar=T,
            bxp.errorbar.width=0.05,
            outlier.shape=NA,
            short.panel.labs = FALSE,
            order = c("NC","AD","Acne"),
            width = 0.7) +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test") +
  stat_compare_means(label.y = 26.1) + 
  ylab("log2(intensity)") + 
  scale_fill_manual(values=color)
p




