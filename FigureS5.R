# Figure S5A,B,C ---------------------------------------------------------------
# Load required packages and data ----------------------------------------------
library(pheatmap)
library(ComplexHeatmap)
library(magick)
library(circlize)

data = read.table("antiMash_summary.txt",header = T,row.names = 1,check.names = F)
#data = read.csv("plasmid_summary.csv",header = T,row.names = 1,check.names = F)
#data = read.csv("phage_taxonomy_summary.csv",header = T,row.names = 1,check.names = F)

# Determine whether the sample order is consistent
all(rownames(groupinfo) %in% colnames(data))
all(rownames(groupinfo) == colnames(data))

# Reorder the sample
data = data[, rownames(groupinfo)]
all(rownames(groupinfo) == colnames(data))

#
groupinfo$Group = factor(groupinfo$Group,levels = c("NC", "AD", "Acne","ConwillA_2022","ZhouW_2020","JoglekarP_2023","KashafSS_2023"))

# color palette
col_fun = colorRamp2(c(0,1),c("#DCDCDC","#B73508"))

# annotation
hac = HeatmapAnnotation(df = groupinfo,
                        col = list(Group=c("NC" = "#90A6BB", "AD" = "#FFD579","Acne"="#9D4E3F","ConwillA_2022"="#B2DF8AFF","ZhouW_2020"="#FB9A99FF","JoglekarP_2023"="#FFFF99FF","KashafSS_2023"="#60A5DFFF")),
                        which = "col"
)

# plot
data = as.matrix(data)
p = ComplexHeatmap::Heatmap(df, 
                          col = col_fun,
                          row_names_side = "right",
                          show_row_names = T,
                          column_names_gp = gpar(fontsize = 10),
                          show_column_names = F,
                          cluster_rows = T,cluster_columns = T,
                          column_split = as.factor(groupinfo$Group),
                          cluster_column_slices = F,
                          row_names_max_width = unit(10,"cm"),
                          rect_gp = gpar(col = "grey", lwd = 0.000001),
                          top_annotation = hac,
                          #left_annotation=ha,
                          #width = ncol(data)*unit(6, "mm"), #
                          #height = nrow(t(df))*unit(0.04, "mm"),#
                          show_heatmap_legend = F
                          
)

p

