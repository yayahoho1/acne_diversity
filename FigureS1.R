# Figure S1A -------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggplot2)

alldata = read.table("SLST_site_stat_for_Group.txt",header = T, sep = "\t",check.names = F)

# color palette
color = c("#BCEE68", "#A2CD5A", "#6E8B3D", "#CDC673", "#8B864E", "#FFB6C1", "#CD8C95", "#8B5F65", "#FFA07A", "#CD8162",
         "#8B5742", "#87CEFA", "#A4D3EE", "#607B8B", "#1E90FF", "#3A5FCD", "#C7C7C7FF")

#
alldata$Group = factor(alldata$Group,levels = c("NC","AD","Acne"))

# plot
p = ggplot(alldata,aes(x=Site,y=Strain_Number,fill=SLST)) +
  geom_bar(stat='identity',width = 0.6) +
  ylab("Strain Number") +
  xlab("") +
  labs(fill = "SLST") +
  facet_grid(~Group) +
  scale_y_continuous(expand = c(0,0.1),limits = c(0,170)) +
  scale_fill_manual(values = color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.text = element_text(size=9,color='black'), 
    plot.title = element_text(size=15,color="black"), 
    legend.key = element_blank(), 
    axis.text = element_text(size=10,color='black'), 
    strip.text = element_text(size=12,color="black"), 
    strip.background = element_blank())


p


# Figure S1B -------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggtree)
library(ggplot2)

tree = read.tree("midpointroot.tree")
sampleinfo = read.table("strain_info.txt",header = T)

# color palette
#colors = c("#ff005e","#174656", "#779966", "#992222","#feb2a3", "#447755", "#aa9977",
#         "#cbc0ae", "#e2c6ff", "#9999cc","#ffdd9a", "#fee2cb", "#ddd6c9","#e5ee99", 
#         "#bbbb77","#fff993", "#ffaa66", "#888855","#bbddee", "#99aa88", "#ddcc88",
#         "#4c5e68","#ffa200", "#eaff00","#00ff00","#0092ff", "#7700ff", "#aa00ff",
#         "#c2c2c2","#7a7a7a","#3c3a7a", "#c1437a")

#color = c("#C79757FF", "#8D41C4FF", "#6DC954FF", "#CA5896FF", "#CACF4FFF", "#7374C8FF", 
#        "#D24F3BFF", "#87CFACFF", "#442F56FF", "#4A603BFF", "#A0A9C0FF", "#7C403AFF" )  

color = c("#A6CEE3FF","#1F78B4FF","#B2DF8AFF", "#33A02CFF", "#FB9A99FF", "#E31A1CFF", 
        "#FDBF6FFF","#FF7F00FF","#CAB2D6FF","#6A3D9AFF","#FFFF99FF","#B15928FF")

# plot
p = ggtree(tree, layout="circular", ladderize = T,size=0.01) +  
  geom_tiplab(aes(label=NA),align=F) + # no label
  theme(legend.position = "right",legend.text=element_text(size=rel(1))) +
  geom_treescale()

p
# annotation
p1 = p %<+% sampleinfo + geom_tippoint(aes(color=Individual,shape =Site),size=4) +
  scale_colour_manual(values = color,na.value = "#ffffff",name="Subjects")
p1

