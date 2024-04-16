# Figure S4A,B -----------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ComplexHeatmap)
library(circlize)

data = read.table("gene_prevalence_without_core_gene.txt",row.names = 1,header = T,check.names = F)
groupinfo = read.table("Sample_group_info.txt",header = T,row.names = 1,sep = "\t")

# Determine whether the sample order is consistent
all(rownames(groupinfo) %in% colnames(data))
all(rownames(groupinfo) == colnames(data))

# Reorder the sample
data <- data[, rownames(groupinfo)]
all(rownames(groupinfo) == colnames(data))

# color palette
col_fun = colorRamp2(c(0,0.5,1),c("#e7e7e5","#bbded7","#D24735"))

# annotation
hac = HeatmapAnnotation(df = groupinfo,
                        col = list(Group=c("NC" = "#90A6BB", "AD" = "#FFD579","Acne"="#9D4E3F")),
                        which = "col"
)

# log
data=log(data+1e-7)

#
groupinfo$Group=factor(groupinfo$Group,levels=c("NC", "AD", "Acne"))

# plot
data=as.matrix(data)
ComplexHeatmap::Heatmap(data, 
                        col = col_fun,
                        row_names_side = "right",
                        show_row_names = F,
                        column_names_gp = gpar(fontsize = 7),
                        show_column_names = T,
                        cluster_rows = T,
                        cluster_columns = T,
                        column_split = as.factor(groupinfo$Group),
                        cluster_column_slices = F,
                        row_names_max_width = unit(10,"cm"),
                        rect_gp = gpar(col = NA, lwd = 0.000001),
                        top_annotation = hac,
                        #left_annotation=ha,
                        heatmap_legend_param = list(title= "log(prevalence+1e-7)",
                                                    title_position = "topcenter", 
                                                    legend_direction= "horizontal",
                                                    legend_height=unit(5,"cm"))
)



# Figure S2C -------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(ComplexHeatmap)
library(circlize)

data=read.table("individual-specific_gene_prevalence.txt",row.names = 1,header = T,check.names = F)

# color palette
col_fun = colorRamp2(c(0,0.5,1),c("#e7e7e5","#bbded7","#D24735"))

#
groupinfo$Group=factor(groupinfo$Group,levels=c("NC", "AD", "Acne"))

# plot
data=as.matrix(data)
ComplexHeatmap::Heatmap(data, 
                        col = col_fun,
                        row_names_side = "right",
                        show_row_names = F,
                        column_names_gp = gpar(fontsize = 7),
                        show_column_names = T,
                        cluster_rows = T,
                        cluster_columns = T,
                        column_split = as.factor(groupinfo$Group),
                        cluster_column_slices = F,
                        row_names_max_width = unit(10,"cm"),
                        rect_gp = gpar(col = NA, lwd = 0.000001),
                        heatmap_legend_param = list(title= "Prevalence",
                                                    title_position = "topcenter", 
                                                    legend_direction= "horizontal",
                                                    legend_height=unit(5,"cm"))
)


































