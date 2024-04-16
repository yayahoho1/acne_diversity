# Figure 1B --------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggtree)
library(ggplot2)
library(ggnewscale)

# read tree
tree = read.tree("RAxML.midpointroot.tree")

# read sample infomation
sampleinfo = read.table("all_MLST_SLST_with_groupinfo.txt",header = T,row.names = 1)

# plot tree
p = ggtree(tree, layout="fan", ladderize = T,size=0.1,open.angle=10)+
  geom_tiplab(aes(label=NA),align=F)+ # no label
  theme(legend.position = "right",legend.text=element_text(size=rel(1)))+
  geom_treescale()

# Group
groupinfo = sampleinfo[4]
groupinfo$Group = factor(groupinfo$Group,levels=c("NC","AD","Acne"))

p1 = gheatmap(p, groupinfo, color=NA,offset = .01, width=0.03, font.size=3, colnames_angle=0, hjust=-.001,
            colnames_position="top",colnames_offset_y = .9 ) + 
  theme(legend.position = "right",legend.text=element_text(size=rel(1)),legend.key.size=unit(0.4,'cm'))

# color palette
color1 = c("#9D4E3F","#FFD579","#90A6BB")

#
p2 = p1 + scale_fill_manual(values = color1,na.value = "#ffffff",name="Group")
p2

# Skin Site
siteinfo=sampleinfo[5]

p3 = p2 + new_scale_fill()
p4 = gheatmap(p3,siteinfo, color=NA,offset = 0.03, width=0.03, font.size=3, colnames_angle=0, hjust=-.001,
            colnames_position="top",colnames_offset_y = .9 ) + 
  theme(legend.position = "right",legend.text=element_text(size=rel(1)),legend.key.size=unit(0.4,'cm'))

# color palette
color2 = c("#626C83","#D5A551","#C3735C","#999999")

p5 = p4 + scale_fill_manual(values = color2,na.value = "#ffffff",name="Skin sites")

# MLST
stinfo=sampleinfo[1]
head(stinfo)

p6 = p5 + new_scale_fill()
p7 = gheatmap(p6,stinfo, color=NA,offset = 0.05, width=0.03, font.size=3, colnames_angle=0, hjust=-.001,
            colnames_position="top",colnames_offset_y = .9 ) + 
  theme(legend.position = "right",legend.text=element_text(size=rel(1)),legend.key.size=unit(0.4,'cm'))

# color palette
color3=c("#EEEED1","#FAEBD7", "#CDC0B0", "#8B8378", "#76EEC6", "#66CDAA", "#458B74", "#E0EEEE", "#C1CDCD",
         "#838B8B", "#98F5FF", "#7AC5CD", "#53868B", "#FFB90F", "#CD950C", "#8B6508")

p8=p7 + scale_fill_manual(values = color3,na.value = "#ffffff",name="MLST")
p8

# SLST
slstinfo=sampleinfo[6]
head(slstinfo)


p9 = p8 + new_scale_fill()
p10 = gheatmap(p9,slstinfo, color=NA,offset = 0.07, width=0.03, font.size=3, colnames_angle=0, hjust=-.001,
               colnames_position="top",colnames_offset_y = .9 ) + 
  theme(legend.position = "right",legend.text=element_text(size=rel(1)),legend.key.size=unit(0.4,'cm'))

# color palette
color4 = c("#BCEE68", "#A2CD5A", "#6E8B3D", "#CDC673", "#8B864E", "#FFB6C1", "#CD8C95", "#8B5F65", "#FFA07A", "#CD8162",
         "#8B5742", "#87CEFA", "#A4D3EE", "#607B8B", "#1E90FF", "#3A5FCD", "#27408B", "#8B668B", "#CD5555")

p11 = p10 + scale_fill_manual(values = color4,na.value = "#ffffff",name="SLST")
p11


# Figure 1C --------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggplot2)
library(ggh4x) # Nested axes

# read data
alldata = read.table("SLST_indivadual_site_stat_with_groupinfo.txt",header = T, sep = "\t")
#
alldata$symbol = factor(alldata$symbol, levels = c("NC","AD","Acne"))

# color palette
color4 = c("#BCEE68", "#A2CD5A", "#6E8B3D", "#CDC673", "#8B864E", "#FFB6C1", "#CD8C95", "#8B5F65", "#FFA07A", "#CD8162",
         "#8B5742", "#87CEFA", "#A4D3EE", "#607B8B", "#1E90FF", "#3A5FCD", "grey60","#27408B", "#8B668B", "#CD5555")

#color4 = c("#607B8B","#CD950C","#458B74","#CD5555","#9467BDFF","#8C564BFF","#E377C2FF",
#          "#BCBD22FF","#17BECFFF","#AEC7E8FF","#FFBB78FF", "#98DF8AFF","#C7C7C7FF","#FF9896FF","#C5B0D5FF",
#          "#C49C94FF","#F7B6D2FF","#DBDB8DFF","#9EDAE5FF","#7F7F7FFF","#C7C7C7FF")

# plot
p = ggplot(alldata,aes(x=interaction(Site,Indivadual),y=number,fill=Type))+
  geom_bar(stat='identity',width=0.9)+
  ylab("Strain Number")+
  xlab("")+
  facet_wrap(~symbol,scales = "free_x")+
  labs(fill = "SLST")+ 
  theme_bw()+
  guides(x = "axis_nested")+ # Nested axes
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size=12),
    legend.text = element_text(size=12),
    ggh4x.axis.nesttext.x = element_text(colour = "blue",size=12),
    axis.text.x = element_text(size=10,color='black',angle = 90,vjust = 0.5,hjust = 0.5)
  )+
  scale_y_continuous(expand = c(0,0.1),limits = c(0,21))+
  scale_fill_manual(values = color4)


# Figure 1D --------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ComplexHeatmap)
library(circlize)

# read data
data = read.csv("SLST_matrix.csv",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
data = t(data)
groupinfo = read.table("sample_group_info.txt",row.names = 1,header = 1)
head(groupinfo)

# Determine whether the sample order is consistent
all(rownames(groupinfo) %in% colnames(data))
all(rownames(groupinfo) == colnames(data))

# Reorder the sample
data = data[, rownames(groupinfo)]
all(rownames(groupinfo) == colnames(data))

#
data = as.matrix(data)
prevalence = t(t(data)/colSums(data))

#
groupinfo$Group=factor(groupinfo$Group,levels=c("NC", "AD", "Acne"))

# color palette
col_fun = colorRamp2(c(0, 1), c("white", "red"))
col_fun(c(0, 0.5, 1))

# plot
ComplexHeatmap::Heatmap(t(prevalence), 
                          col = col_fun,
                          row_names_side = "right",
                          show_row_names = T,
                          column_names_gp = gpar(fontsize = 10),
                          show_column_names = T,
                          cluster_rows = T,cluster_columns = T,
                          row_split = as.factor(groupinfo$Group),
                          cluster_row_slices = F,
                          row_names_max_width = unit(10,"cm"),
                          rect_gp = gpar(col = "grey70", lwd = 0.000001),
                          width = ncol(t(prevalence))*unit(6, "mm"),
                          height = nrow(t(prevalence))*unit(5, "mm"),
                          heatmap_legend_param = list(title= "Prevalence",
                                                      title_position = "topcenter", 
                                                      legend_height=unit(5,"cm"))
)






