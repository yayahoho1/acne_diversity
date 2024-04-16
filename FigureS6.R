# Figure S6A,B,C ---------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ggsignif)
library(ggpubr)
library(ggplot2)
library(egg)

data = read.table("HGT_belong_to_accessory_summary.txt",header = T,sep = "\t")

# color palette
color = c("#ab6f26","#d6b474")

#
my_comparisons = list(c("high","low"))

# plot
p = ggboxplot(data, x="Group", y="Accessory.gene.number",color="black",
            fill="Group", 
            palette = "npg",
            #add = "jitter",
            bxp.errorbar=T,
            bxp.errorbar.width=0.05,
            outlier.shape=NA, 
            short.panel.labs = FALSE,
            width = 0.7) +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test") + 
  ylab("The number of HGT belonging to accessory genes") + 
  scale_fill_manual(values=color)
p



# Figure S6D -------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

#install.packages("magick")
library(magick) #
library(circlize)
library(ComplexHeatmap)

data=read.table("potential_donor_taxonID.txt",header = T,row.names = 1,sep = "\t",check.names = F)
groupinfo=read.table("group_info.txt",header = T,row.names = 1,sep = "\t")


# Determine whether the sample order is consistent
all(rownames(groupinfo) %in% colnames(data))
all(rownames(groupinfo) == colnames(data))
# Reorder the sample；
data <- data[, rownames(groupinfo)]
all(rownames(groupinfo) == colnames(data))

#
groupinfo$Group = factor(groupinfo$Group,levels = c("NC", "AD", "Acne"))

# color palette
col_fun = colorRamp2(c(0,1),c("#DCDCDC","#B73508"))

# annotation
hac = HeatmapAnnotation(df = groupinfo,
                        col = list(Group=c("NC" = "#90A6BB", "AD" = "#FFD579","Acne"="#9D4E3F")),
                        which = "col"
)

# plot
df=as.matrix(data)
ComplexHeatmap::Heatmap(df, 
                        col = col_fun,
                        row_names_side = "right",
                        show_row_names = T,
                        column_names_gp = gpar(fontsize = 10),
                        row_names_gp = gpar(fontsize = 9),
                        show_column_names = F,
                        cluster_rows = T,cluster_columns = T,
                        column_split = as.factor(groupinfo$Group),
                        cluster_column_slices = F,
                        row_names_max_width = unit(10,"cm"),
                        #rect_gp = gpar(col = "grey", lwd = 0.000001),
                        top_annotation = hac,
                        width = ncol(data)*unit(0.1, "mm"), #
                        height = nrow(t(df))*unit(0.2, "mm"),#
                        show_heatmap_legend = F,
                        
)




