# Figure S9A,B -----------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(dplyr)
library(ggplot2)
library(ggbreak)
library(ggplot2)
library(reshape2)

data = read.csv("significant_SNP_metabolome_summary_top50.csv",header = T)
#data = read.csv("significant_SNP_transcriptome_summary_top50.csv",header = T)

#
data$significant_SNP_number=as.numeric(data$significant_SNP_number)

# plot
p <- ggplot(df,aes(x=reorder(Metabolome,significant_SNP_number),y=significant_SNP_number)) +
  geom_bar(stat='identity',fill="#DD604B",width = 0.7) +
  ylab("significant SNP number") + 
  xlab("Metabolome") +
  #theme_bw()+
  coord_flip() +
  mytheme<-theme_classic()+
  theme(text=element_text(family = "sans",colour ="gray30",size = 12),
        legend.text=element_text(colour ="gray30",size = 8),
        legend.title=element_text(colour ="gray30",size = 10),
        legend.key.size=unit(4,units = "mm"),
        legend.position=c(0.10,0.88),
        axis.line = element_line(size = 0.4,colour = "gray30"),
        axis.ticks = element_line(size = 0.4,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"),
        axis.text.y = element_text(size=10,color='black'),
        axis.text.x = element_text(size=10,color='black')
        #plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
        #units="inches")
  ) +
  scale_y_continuous(expand = c(0,5))

p

#
p + scale_y_break(c(15000,30000),scales = .2,space = 0.3)


# Figure S9C -------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(circlize)

data = read.csv("edges.csv",header = T)
groupinfo = read.csv("nodes_order.csv",header = T)

#
data$sig_SNP_number = 1

#
chordDiagram(data, 
             annotationTrack = c('grid','name'),
             annotationTrackHeight = c(0.05, 0.01)
)

#
par(cex = 0.3, mar = c(0, 0, 0, 0))

#
nodes = groupinfo$name
orders = groupinfo$Order
group = structure(orders,names=as.character(nodes))

# edit initialising parameters
#
circos.par(canvas.ylim=c(-1.5,1.5), # edit  canvas size 
           track.margin = c(0.01, 0.01), # adjust bottom and top margin
           track.height = 0.01)

#
chordDiagram(data, group=group, annotationTrack = "grid",preAllocateTracks = 1,annotationTrackHeight = 0.01,big.gap=1,small.gap =0.1 )

#
circos.trackPlotRegion(track.index=1, panel.fun=function(x,y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name=get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name,facing="clockwise",
              niceFacing=TRUE,adj=c(0,0.5))},bg.border=NA)

#
Genome = groupinfo[groupinfo$Group == "Genome",]
Transcriptome = groupinfo[groupinfo$Group == "Transcriptome",]
Metabolome = groupinfo[groupinfo$Group == "Metabolome",]

#
highlight.sector(
  Genome$name, track.index = 1, col = "#FF6347", 
  text = "Genome", cex = 2, text.col = "red", 
  niceFacing = TRUE
  #lwd = par(2)
)

#
highlight.sector(
  Transcriptome$name, track.index = 1, col = "#80b1d3", 
  text = "Transcriptome", cex = 2, text.col = "red", 
  niceFacing = TRUE
)

#
highlight.sector(
  Metabolome$name, track.index = 1, col = "#fdb462", 
  text = "Metabolome", cex = 2, text.col = "red", 
  niceFacing = TRUE
)

#
circos.clear()
















